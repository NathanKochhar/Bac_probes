"""Batch probe design for every genus in KAN2M3-genus.csv (no BLAST).

Equivalent to `bac-probes find --no-blast` per genus (same k-mer filters and
exact off-target counting), but reads SILVA once and scores all genera in
parallel instead of 2 full SILVA passes per genus.

Selection goal: ~10 highly specific probes whose POOLED coverage is ~50%.
  precision = target hits / (target hits + exact off-target hits)
  For precision floors from strict (1.0 = zero off-targets) to loose, run a
  greedy set cover (max new target seqs covered, up to N_PROBES, skipping
  k-mers that share a 16-mer with a picked probe), then top up to N_PROBES with
  redundant non-overlapping probes if coverage saturates. Use the strictest floor
  whose pool reaches POOL_TARGET coverage; otherwise the floor that covers most.
"""
import multiprocessing as mp
from collections import Counter
from pathlib import Path

import pandas as pd

from bac_probes.database import iter_silva
from bac_probes.kmers import VALID_BASES, extract_kmers, gc_content, reverse_complement

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
SILVA = ROOT / "silva_db" / "SILVA_138.1_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz"
INPUT_CSV = OUT.parent / "KAN2M3-genus.csv"

K = 32
GC_RANGE = (0.35, 0.65)
MAX_HOMOPOLYMER = 5
MIN_CONS = 0.02          # candidate floor; low so rare-but-specific k-mers are kept
N_PROBES = 10
POOL_TARGET = 0.50
PRECISION_FLOORS = [1.0, 0.999, 0.995, 0.99, 0.98, 0.95, 0.90]
OVERLAP_WORD = 16
WORKERS = 14

# CSV name -> SILVA 138.1 genus name (only where they differ)
SILVA_NAME = {"Clostridium": "Clostridium sensu stricto 1"}


# ── worker functions ──────────────────────────────────────────────────────────

def conserve(args):
    """Candidate k-mers for one genus with the target seq indices containing each."""
    genus, seqs = args
    per_seq = [extract_kmers(s, K, gc_range=GC_RANGE, max_homopolymer=MAX_HOMOPOLYMER) for s in seqs]
    counts = Counter()
    for ks in per_seq:
        counts.update(ks)
    need = max(1, MIN_CONS * len(seqs))
    cands = {k for k, c in counts.items() if c >= need}
    cover = {k: [] for k in cands}
    for i, ks in enumerate(per_seq):
        for k in ks & cands:
            cover[k].append(i)
    return genus, cover


_OWNERS = None


def _init(owners):
    global _OWNERS
    _OWNERS = owners


def offtarget(chunk):
    """Exact off-target counts {(genus, kmer): n} for a chunk of (genus, seq)."""
    counts = Counter()
    for seq_genus, seq in chunk:
        seen = set()
        for i in range(len(seq) - K + 1):
            kmer = seq[i:i + K]
            if kmer in _OWNERS and kmer not in seen:
                seen.add(kmer)
        for kmer in seen:
            for g in _OWNERS[kmer]:
                if g != seq_genus:
                    counts[(g, kmer)] += 1
    return counts


def top_offtargets(chunk_and_probes):
    """Off-target genus tallies for the selected probes only."""
    chunk, probes = chunk_and_probes
    tallies = {}
    for seq_genus, seq in chunk:
        for g, kmer in probes:
            if seq_genus != g and kmer in seq:
                tallies.setdefault((g, kmer), Counter())[seq_genus or "unassigned"] += 1
    return tallies


# ── selection ─────────────────────────────────────────────────────────────────

def words(seq):
    return {seq[i:i + OVERLAP_WORD] for i in range(len(seq) - OVERLAP_WORD + 1)}


def greedy(df, cover, n_target, floor):
    pool = df[df.precision >= floor]
    sets = {k: set(cover[k]) for k in pool.kmer}
    prec = dict(zip(pool.kmer, pool.precision))
    covered, seen, picked = set(), set(), []
    while len(picked) < N_PROBES and sets:
        best = max(sets, key=lambda k: (len(sets[k] - covered), prec[k]))
        if not sets[best] - covered:
            break
        new = sets.pop(best)
        if words(best) & seen:
            continue
        picked.append((best, len(new - covered)))
        covered |= new
        seen |= words(best)
    # Coverage saturated before N_PROBES: top up with extra non-overlapping probes
    # from the same precision pool (more signal per cell, no specificity cost).
    for k in pool.sort_values(["precision", "conservation_pct"], ascending=False).kmer:
        if len(picked) >= N_PROBES:
            break
        if k in sets and not words(k) & seen:
            picked.append((k, 0))
            seen |= words(k)
    return picked, len(covered) / n_target


def main():
    csv_names = pd.read_csv(INPUT_CSV)["Taxa"].str.strip().tolist()
    targets = {SILVA_NAME.get(n, n): n for n in csv_names}

    print("Reading SILVA …", flush=True)
    records, target_seqs = [], {g: [] for g in targets}
    for _acc, tax, seq in iter_silva(SILVA):
        g = tax.get("genus", "")
        records.append((g, seq))
        if g in target_seqs:
            target_seqs[g].append(seq)
    print(f"  {len(records):,} sequences", flush=True)

    with mp.Pool(WORKERS) as pool:
        print("Scoring conservation per genus …", flush=True)
        jobs = sorted(target_seqs.items(), key=lambda kv: -len(kv[1]))
        covers = dict(pool.imap_unordered(conserve, jobs))

    owners = {}
    for g, cover in covers.items():
        for k in cover:
            owners.setdefault(k, []).append(g)
    print(f"  {len(owners):,} candidate k-mers across {len(covers)} genera", flush=True)

    print("Exact off-target scan …", flush=True)
    step = 5000
    chunks = [records[i:i + step] for i in range(0, len(records), step)]
    off = Counter()
    with mp.Pool(WORKERS, initializer=_init, initargs=(owners,)) as pool:
        for c in pool.imap_unordered(offtarget, chunks):
            off.update(c)

    (OUT / "find_runs").mkdir(exist_ok=True)
    (OUT / "selected").mkdir(exist_ok=True)
    n_total = len(records)
    selections, summary = {}, []
    for g, csv_name in targets.items():
        n_t = len(target_seqs[g])
        n_bg = n_total - n_t
        cover = covers[g]
        rows = []
        for k, idx in cover.items():
            hits, ot = len(idx), off.get((g, k), 0)
            rows.append({
                "kmer": k, "probe": reverse_complement(k),
                "gc_pct": round(gc_content(k) * 100, 1),
                "conservation_pct": round(100 * hits / n_t, 2),
                "target_seqs": n_t, "target_hits": hits,
                "exact_offtarget": ot, "exact_bg_seqs": n_bg,
                "exact_specificity": round(1 - ot / n_bg, 6),
                "precision": round(hits / (hits + ot), 4),
            })
        df = pd.DataFrame(rows).sort_values(["precision", "conservation_pct"], ascending=False)
        safe = g.replace(" ", "_")
        df.to_csv(OUT / "find_runs" / f"{safe}_find.tsv", sep="\t", index=False)

        best = None
        for floor in PRECISION_FLOORS:
            picked, pooled = greedy(df, cover, n_t, floor)
            if best is None or pooled > best[2] + 1e-9:
                best = (floor, picked, pooled)
            if pooled >= POOL_TARGET:
                best = (floor, picked, pooled)
                break
        floor, picked, pooled = best
        sel = df.set_index("kmer").loc[[k for k, _ in picked]].reset_index()
        sel.insert(0, "new_seqs_covered", [n for _, n in picked])
        sel["cumulative_coverage_pct"] = (100 * sel.new_seqs_covered.cumsum() / n_t).round(2)
        sel.insert(0, "probe_rank", range(1, len(sel) + 1))
        sel.insert(0, "silva_genus", g)
        sel.insert(0, "taxa", csv_name)
        selections[g] = sel
        summary.append({
            "taxa": csv_name, "silva_genus": g, "target_seqs": n_t,
            "candidates_found": len(df), "precision_floor_used": floor,
            "probes_selected": len(sel), "pooled_coverage_pct": round(100 * pooled, 2),
            "conservation_min_pct": sel.conservation_pct.min() if len(sel) else None,
            "conservation_max_pct": sel.conservation_pct.max() if len(sel) else None,
            "precision_min": sel.precision.min() if len(sel) else None,
            "precision_mean": round(sel.precision.mean(), 4) if len(sel) else None,
            "exact_offtarget_total": int(sel.exact_offtarget.sum()) if len(sel) else None,
            "exact_offtarget_max": int(sel.exact_offtarget.max()) if len(sel) else None,
        })
        print(f"  {g}: floor {floor} -> {len(sel)} probes, pooled {100 * pooled:.1f}%", flush=True)

    print("Tallying off-target genera for selected probes …", flush=True)
    probes = [(g, k) for g, s in selections.items() for k in s.kmer]
    tallies = {}
    with mp.Pool(WORKERS) as pool:
        for t in pool.imap_unordered(top_offtargets, [(c, probes) for c in chunks]):
            for key, cnt in t.items():
                tallies.setdefault(key, Counter()).update(cnt)

    def fmt(g, k):
        c = tallies.get((g, k))
        return "; ".join(f"{name} ({n})" for name, n in c.most_common(3)) if c else ""

    all_sel = []
    for g, sel in selections.items():
        sel["top_offtarget_genera"] = [fmt(g, k) for k in sel.kmer]
        sel.to_csv(OUT / "selected" / f"{g.replace(' ', '_')}_selected_probes.csv", index=False)
        all_sel.append(sel)
    pd.concat(all_sel).to_csv(OUT / "Genus_all_selected_probes.csv", index=False)

    s = pd.DataFrame(summary)
    s["top_offtarget_genera"] = [
        "; ".join(f"{n} ({c})" for n, c in sum(
            (tallies.get((r.silva_genus, k), Counter()) for k in selections[r.silva_genus].kmer), Counter()
        ).most_common(3))
        for r in s.itertuples()
    ]
    s.to_csv(OUT / "Genus_summary.csv", index=False)
    print(s.drop(columns=["top_offtarget_genera"]).to_string(index=False))


if __name__ == "__main__":
    main()
