"""Batch probe design for every family in KAN2M3-family.csv, with BLAST scoring.

Stage 1 (exact, same as the genus run): read SILVA once, find k-mers present in
>= MIN_CONS of each family's sequences, count exact off-target hits in one
parallel pass. Same filters as `bac-probes find`.

Stage 2 (BLAST): per family, BLAST a pool of up to 2 x POOL_HALF candidates
at distinct loci (no shared 16-mer): half chosen by coverage, half by exact
precision, all with exact precision >= MIN_EXACT_PREC. Uses the package's
run_blast / parse_blast_results (blastn-short, position-weighted near-miss
scoring). Results are cached per family, so a rerun skips finished families.

Selection goal: ~10 highly specific probes with POOLED coverage ~50%.
  Specificity = blast_weighted_specificity (exact off-targets count fully;
  near-misses weighted by mismatch position).
  For floors from strict to loose, greedy set cover (max new target seqs,
  up to N_PROBES, no shared 16-mers), topped up to N_PROBES if coverage
  saturates. Use the strictest floor reaching POOL_TARGET, else best coverage.
"""
import json
import multiprocessing as mp
import pickle
import time
from collections import Counter
from pathlib import Path

import pandas as pd

from bac_probes.database import iter_silva
from bac_probes.kmers import extract_kmers, gc_content, reverse_complement
from bac_probes.specificity import parse_blast_results, run_blast

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
SILVA = ROOT / "silva_db" / "SILVA_138.1_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz"
BLAST_DB = ROOT / "silva_db" / "silva"
INPUT_CSV = OUT.parent / "KAN2M3-family.csv"
CACHE = OUT / "cache"

LEVEL = "family"
K = 32
GC_RANGE = (0.35, 0.65)
MAX_HOMOPOLYMER = 5
MIN_CONS = 0.02
MIN_EXACT_PREC = 0.90
POOL_HALF = 60
N_PROBES = 10
POOL_TARGET = 0.50
SPEC_FLOORS = [0.999, 0.995, 0.99, 0.98, 0.95, 0.90]
OVERLAP_WORD = 16
WORKERS = 14
BLAST_THREADS = 14
BLAST_IDENTITY = 85.0
BLAST_MAX_TARGET_SEQS = 50000
BLAST_CHUNK = 20


# ── stage 1 workers ───────────────────────────────────────────────────────────

def conserve(args):
    taxon, seqs = args
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
    return taxon, cover


_OWNERS = None


def _init(owners):
    global _OWNERS
    _OWNERS = owners


def offtarget(chunk):
    counts = Counter()
    for seq_taxon, seq in chunk:
        seen = set()
        for i in range(len(seq) - K + 1):
            kmer = seq[i:i + K]
            if kmer in _OWNERS:
                seen.add(kmer)
        for kmer in seen:
            for t in _OWNERS[kmer]:
                if t != seq_taxon:
                    counts[(t, kmer)] += 1
    return counts


def stage1(targets):
    ck = CACHE / "stage1.pkl"
    if ck.exists():
        print("Stage 1: loading cached exact results", flush=True)
        return pickle.loads(ck.read_bytes())

    print("Stage 1: reading SILVA …", flush=True)
    records, target_seqs = [], {t: [] for t in targets}
    for _acc, tax, seq in iter_silva(SILVA):
        t = tax.get(LEVEL, "")
        records.append((t, seq))
        if t in target_seqs:
            target_seqs[t].append(seq)
    print(f"  {len(records):,} sequences", flush=True)

    with mp.Pool(WORKERS) as pool:
        print("  scoring conservation …", flush=True)
        jobs = sorted(target_seqs.items(), key=lambda kv: -len(kv[1]))
        covers = dict(pool.imap_unordered(conserve, jobs))

    owners = {}
    for t, cover in covers.items():
        for k in cover:
            owners.setdefault(k, []).append(t)
    print(f"  {len(owners):,} candidate k-mers", flush=True)

    print("  exact off-target scan …", flush=True)
    chunks = [records[i:i + 5000] for i in range(0, len(records), 5000)]
    off = Counter()
    with mp.Pool(WORKERS, initializer=_init, initargs=(owners,)) as pool:
        for c in pool.imap_unordered(offtarget, chunks):
            off.update(c)

    n_total = len(records)
    tables = {}
    for t in targets:
        n_t = len(target_seqs[t])
        n_bg = n_total - n_t
        rows = []
        for k, idx in covers[t].items():
            hits, ot = len(idx), off.get((t, k), 0)
            rows.append({
                "kmer": k, "probe": reverse_complement(k),
                "gc_pct": round(gc_content(k) * 100, 1),
                "conservation_pct": round(100 * hits / n_t, 2),
                "target_seqs": n_t, "target_hits": hits,
                "exact_offtarget": ot, "exact_bg_seqs": n_bg,
                "exact_specificity": round(1 - ot / n_bg, 6),
                "exact_precision": round(hits / (hits + ot), 4),
            })
        tables[t] = pd.DataFrame(rows).sort_values(
            ["exact_precision", "conservation_pct", "kmer"], ascending=[False, False, True])
    result = (tables, covers)
    CACHE.mkdir(exist_ok=True)
    ck.write_bytes(pickle.dumps(result))
    return result


# ── stage 2: BLAST ────────────────────────────────────────────────────────────

def words(seq):
    return {seq[i:i + OVERLAP_WORD] for i in range(len(seq) - OVERLAP_WORD + 1)}


def distinct(kmers, limit, seen=None):
    seen = set() if seen is None else seen
    out = []
    for k in kmers:
        w = words(k)
        if w & seen:
            continue
        out.append(k)
        seen |= w
        if len(out) == limit:
            break
    return out


def blast_pool(df):
    ok = df[df.exact_precision >= MIN_EXACT_PREC]
    # every sort ends on the k-mer itself so ties break the same way on every run
    by_cov = distinct(ok.sort_values(["conservation_pct", "kmer"], ascending=[False, True]).kmer, POOL_HALF)
    by_prec = distinct(ok.sort_values(["exact_precision", "conservation_pct", "kmer"],
                                      ascending=[False, False, True]).kmer,
                       POOL_HALF, seen=set().union(*map(words, by_cov)) if by_cov else set())
    return by_cov + by_prec


def blast_family(taxon, kmers):
    ck = CACHE / f"blast_{taxon}.json"
    scores = json.loads(ck.read_text()) if ck.exists() else {}
    todo = [k for k in kmers if k not in scores]  # only BLAST candidates not cached yet
    if not todo:
        return {k: scores[k] for k in kmers}
    for i in range(0, len(todo), BLAST_CHUNK):
        chunk = todo[i:i + BLAST_CHUNK]
        out = run_blast(chunk, BLAST_DB, threads=BLAST_THREADS,
                        perc_identity=BLAST_IDENTITY, max_target_seqs=BLAST_MAX_TARGET_SEQS)
        scores.update(parse_blast_results(out, chunk, taxon, LEVEL,
                                          max_target_seqs=BLAST_MAX_TARGET_SEQS))
        del out
    for v in scores.values():
        v["blast_capped"] = bool(v["blast_capped"])
    ck.write_text(json.dumps(scores))
    return {k: scores[k] for k in kmers}


# ── selection ─────────────────────────────────────────────────────────────────

def greedy(df, cover, n_target, floor):
    pool = df[df.blast_weighted_specificity >= floor]
    sets = {k: set(cover[k]) for k in pool.kmer}
    spec = dict(zip(pool.kmer, pool.blast_weighted_specificity))
    covered, seen, picked = set(), set(), []
    while len(picked) < N_PROBES and sets:
        best = min(sets, key=lambda k: (-len(sets[k] - covered), -spec[k], k))
        if not sets[best] - covered:
            break
        new = sets.pop(best)
        if words(best) & seen:
            continue
        picked.append((best, len(new - covered)))
        covered |= new
        seen |= words(best)
    for k in pool.sort_values(["blast_weighted_specificity", "conservation_pct", "kmer"],
                              ascending=[False, False, True]).kmer:
        if len(picked) >= N_PROBES:
            break
        if k in sets and not words(k) & seen:
            picked.append((k, 0))
            seen |= words(k)
    return picked, len(covered) / n_target


def main():
    targets = [t for t in pd.read_csv(INPUT_CSV)["Taxa"].dropna().str.strip() if t]
    tables, covers = stage1(targets)

    (OUT / "find_runs").mkdir(exist_ok=True)
    for t in targets:
        tables[t].to_csv(OUT / "find_runs" / f"{t}_find.tsv", sep="\t", index=False)

    pools = {t: blast_pool(tables[t]) for t in targets}
    total = sum(map(len, pools.values()))
    print(f"Stage 2: BLAST {total} candidates across {len(targets)} families "
          f"(identity {BLAST_IDENTITY}%, max_target_seqs {BLAST_MAX_TARGET_SEQS})", flush=True)
    blast = {}
    t0, done = time.time(), 0
    for t in targets:
        blast[t] = blast_family(t, pools[t])
        done += len(pools[t])
        print(f"  {t}: {len(pools[t])} BLASTed  [{done}/{total}, {(time.time() - t0) / 60:.1f} min]", flush=True)

    (OUT / "selected").mkdir(exist_ok=True)
    all_sel, summary = [], []
    print("Selecting …", flush=True)
    for t in targets:
        n_t = int(tables[t].target_seqs.iloc[0])
        df = tables[t][tables[t].kmer.isin(blast[t])].copy()
        for col in ["blast_total_hits", "blast_target_hits", "blast_offtarget_hits",
                    "blast_weighted_offtarget", "blast_specificity",
                    "blast_weighted_specificity", "blast_top_offtarget", "blast_capped"]:
            df[col] = df.kmer.map(lambda k, c=col: blast[t][k][c])
        df.sort_values(["blast_weighted_specificity", "conservation_pct", "kmer"],
                       ascending=[False, False, True]).to_csv(
            OUT / "find_runs" / f"{t}_blast_scored.tsv", sep="\t", index=False)

        best = None
        for floor in SPEC_FLOORS:
            picked, pooled = greedy(df, covers[t], n_t, floor)
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
        sel.insert(0, "taxa", t)
        sel.to_csv(OUT / "selected" / f"{t}_selected_probes.csv", index=False)
        all_sel.append(sel)
        top_ot = Counter(x for x in sel.blast_top_offtarget if isinstance(x, str) and x)
        summary.append({
            "taxa": t, "target_seqs": n_t,
            "candidates_found": len(tables[t]), "candidates_blasted": len(df),
            "blast_spec_floor_used": floor, "probes_selected": len(sel),
            "pooled_coverage_pct": round(100 * pooled, 2),
            "conservation_min_pct": sel.conservation_pct.min() if len(sel) else None,
            "conservation_max_pct": sel.conservation_pct.max() if len(sel) else None,
            "blast_weighted_spec_min": sel.blast_weighted_specificity.min() if len(sel) else None,
            "blast_weighted_spec_mean": round(sel.blast_weighted_specificity.mean(), 4) if len(sel) else None,
            "exact_precision_min": sel.exact_precision.min() if len(sel) else None,
            "exact_offtarget_total": int(sel.exact_offtarget.sum()) if len(sel) else None,
            "blast_offtarget_hits_total": int(sel.blast_offtarget_hits.sum()) if len(sel) else None,
            "probes_blast_capped": int(sel.blast_capped.sum()) if len(sel) else 0,
            "top_offtarget_families": "; ".join(f"{n} ({c})" for n, c in top_ot.most_common(3)),
        })
        print(f"  {t}: floor {floor} -> {len(sel)} probes, pooled {100 * pooled:.1f}%", flush=True)

    pd.concat(all_sel).to_csv(OUT / "Family_all_selected_probes.csv", index=False)
    s = pd.DataFrame(summary)
    s.to_csv(OUT / "Family_summary.csv", index=False)
    print(s.drop(columns=["top_offtarget_families"]).to_string(index=False))


if __name__ == "__main__":
    main()
