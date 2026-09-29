"""Batch probe design for every species in KAN2M3-species.csv, with BLAST scoring.

SILVA species labels carry strain/subspecies suffixes ("Escherichia coli O157:H7",
"Enterobacter hormaechei subsp. oharae"), so every label is normalised to its
first two words ("Genus species").

Each SILVA sequence is classified relative to a target species as:
  target    - labelled with the target species (any strain/subspecies)
  ambiguous - same SILVA genus as the target but no real species name
              ("Fusobacterium sp.", "uncultured bacterium", "uncultured
              Helicobacter sp." ...). Could be the target; excluded from both
              target and off-target counts.
  off       - everything else, including every NAMED sister species.
The same rule is used for exact counting and for BLAST hit classification.

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
from bac_probes.database import parse_taxonomy
from bac_probes.specificity import _hit_mismatch_positions, effective_binding_score, run_blast

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
SILVA = ROOT / "silva_db" / "SILVA_138.1_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz"
BLAST_DB = ROOT / "silva_db" / "silva"
INPUT_CSV = OUT.parent / "KAN2M3-species.csv"
CACHE = OUT / "cache"

LEVEL = "species"
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


def binomial(name):
    """'Escherichia coli O157:H7' -> 'Escherichia coli'."""
    return " ".join(name.split()[:2])


AMBIGUOUS_EPITHETS = {"sp.", "sp", "bacterium", "genomosp.", "cf.", "aff.", "oral"}


def classify(genus, species_raw, target, target_genus):
    """'target', 'ambiguous' or 'off' for one SILVA sequence vs one target species."""
    if binomial(species_raw) == target:
        return "target"
    if genus == target_genus:
        w = species_raw.split()
        named = (len(w) >= 2 and w[0][:1].isupper() and w[1][:1].islower()
                 and w[0].lower() not in ("uncultured", "unidentified", "unclassified", "metagenome")
                 and w[1] not in AMBIGUOUS_EPITHETS)
        if not named:
            return "ambiguous"
    return "off"


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
_TGENUS = None


def _init(owners, tgenus):
    global _OWNERS, _TGENUS
    _OWNERS, _TGENUS = owners, tgenus


def offtarget(chunk):
    counts = Counter()
    for genus, species_raw, seq in chunk:
        seen = set()
        for i in range(len(seq) - K + 1):
            kmer = seq[i:i + K]
            if kmer in _OWNERS:
                seen.add(kmer)
        for kmer in seen:
            for t in _OWNERS[kmer]:
                if classify(genus, species_raw, t, _TGENUS[t]) == "off":
                    counts[(t, kmer)] += 1
    return counts


def stage1(targets):
    ck = CACHE / "stage1.pkl"
    if ck.exists():
        print("Stage 1: loading cached exact results", flush=True)
        return pickle.loads(ck.read_bytes())

    print("Stage 1: reading SILVA …", flush=True)
    records, target_seqs, target_genera = [], {t: [] for t in targets}, {}
    for _acc, tax, seq in iter_silva(SILVA):
        genus, species_raw = tax.get("genus", ""), tax.get(LEVEL, "")
        records.append((genus, species_raw, seq))
        t = binomial(species_raw)
        if t in target_seqs:
            target_seqs[t].append(seq)
            target_genera.setdefault(t, Counter())[genus] += 1
    print(f"  {len(records):,} sequences", flush=True)
    # SILVA genus of each target = most common genus among its sequences
    tgenus = {t: target_genera[t].most_common(1)[0][0] if t in target_genera else "" for t in targets}
    n_ambig = {t: sum(1 for g, sp, _ in records if classify(g, sp, t, tgenus[t]) == "ambiguous")
               for t in targets}
    for t in targets:
        print(f"  {t}: {len(target_seqs[t])} target, SILVA genus '{tgenus[t]}', "
              f"{n_ambig[t]} ambiguous same-genus seqs excluded", flush=True)

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
    with mp.Pool(WORKERS, initializer=_init, initargs=(owners, tgenus)) as pool:
        for c in pool.imap_unordered(offtarget, chunks):
            off.update(c)

    n_total = len(records)
    tables = {}
    for t in targets:
        n_t = len(target_seqs[t])
        n_bg = n_total - n_t - n_ambig[t]
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
            ["exact_precision", "conservation_pct"], ascending=False)
    result = (tables, covers, tgenus, n_ambig)
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
    by_cov = distinct(ok.sort_values("conservation_pct", ascending=False).kmer, POOL_HALF)
    by_prec = distinct(ok.sort_values(["exact_precision", "conservation_pct"], ascending=False).kmer,
                       POOL_HALF, seen=set().union(*map(words, by_cov)) if by_cov else set())
    return by_cov + by_prec


def parse_blast_species(blast_output, kmers, target, tgenus):
    """parse_blast_results with the 3-way target / ambiguous / off classification."""
    idx = {f"kmer_{i}": k for i, k in enumerate(kmers)}
    data = {k: {"target": 0, "ambig": 0, "off": 0, "offw": 0.0, "taxa": Counter()} for k in kmers}
    for line in blast_output.splitlines():
        parts = line.split("\t")
        if len(parts) < 8 or parts[0] not in idx:
            continue
        kmer = idx[parts[0]]
        tax = parse_taxonomy(parts[1])
        cls = classify(tax.get("genus", ""), tax.get("species", ""), target, tgenus)
        d = data[kmer]
        if cls == "target":
            d["target"] += 1
        elif cls == "ambiguous":
            d["ambig"] += 1
        else:
            mm = _hit_mismatch_positions(parts[5], int(parts[6]), int(parts[7]), len(kmer))
            d["off"] += 1
            d["offw"] += 1.0 - effective_binding_score(mm, len(kmer))
            d["taxa"][binomial(tax.get("species", "")) or tax.get("genus", "")] += 1
    out = {}
    for k, d in data.items():
        total = d["target"] + d["off"]
        wtotal = d["target"] + d["offw"]
        top = d["taxa"].most_common(1)
        out[k] = {
            "blast_total_hits": total + d["ambig"],
            "blast_target_hits": d["target"],
            "blast_ambiguous_hits": d["ambig"],
            "blast_offtarget_hits": d["off"],
            "blast_weighted_offtarget": round(d["offw"], 3),
            "blast_specificity": round(d["target"] / total, 6) if total else 1.0,
            "blast_weighted_specificity": round(d["target"] / wtotal, 6) if wtotal else 1.0,
            "blast_top_offtarget": top[0][0] if top else "",
            "blast_capped": total + d["ambig"] >= BLAST_MAX_TARGET_SEQS,
        }
    return out


def blast_family(taxon, kmers, tgenus):
    ck = CACHE / f"blast_{taxon.replace(' ', '_')}.json"
    if ck.exists():
        cached = json.loads(ck.read_text())
        if set(cached) >= set(kmers):
            return cached
    scores = {}
    for i in range(0, len(kmers), BLAST_CHUNK):
        chunk = kmers[i:i + BLAST_CHUNK]
        out = run_blast(chunk, BLAST_DB, threads=BLAST_THREADS,
                        perc_identity=BLAST_IDENTITY, max_target_seqs=BLAST_MAX_TARGET_SEQS)
        scores.update(parse_blast_species(out, chunk, taxon, tgenus))
        del out
    for v in scores.values():
        v["blast_capped"] = bool(v["blast_capped"])
    ck.write_text(json.dumps(scores))
    return scores


# ── selection ─────────────────────────────────────────────────────────────────

def greedy(df, cover, n_target, floor):
    pool = df[df.blast_weighted_specificity >= floor]
    sets = {k: set(cover[k]) for k in pool.kmer}
    spec = dict(zip(pool.kmer, pool.blast_weighted_specificity))
    covered, seen, picked = set(), set(), []
    while len(picked) < N_PROBES and sets:
        best = max(sets, key=lambda k: (len(sets[k] - covered), spec[k]))
        if not sets[best] - covered:
            break
        new = sets.pop(best)
        if words(best) & seen:
            continue
        picked.append((best, len(new - covered)))
        covered |= new
        seen |= words(best)
    for k in pool.sort_values(["blast_weighted_specificity", "conservation_pct"], ascending=False).kmer:
        if len(picked) >= N_PROBES:
            break
        if k in sets and not words(k) & seen:
            picked.append((k, 0))
            seen |= words(k)
    return picked, len(covered) / n_target


def main():
    targets = [t for t in pd.read_csv(INPUT_CSV)["Taxa"].dropna().str.strip() if t]
    tables, covers, tgenus, n_ambig = stage1(targets)

    (OUT / "find_runs").mkdir(exist_ok=True)
    for t in targets:
        tables[t].to_csv(OUT / "find_runs" / f"{t.replace(' ', '_')}_find.tsv", sep="\t", index=False)

    pools = {t: blast_pool(tables[t]) for t in targets}
    total = sum(map(len, pools.values()))
    print(f"Stage 2: BLAST {total} candidates across {len(targets)} species "
          f"(identity {BLAST_IDENTITY}%, max_target_seqs {BLAST_MAX_TARGET_SEQS})", flush=True)
    blast = {}
    t0, done = time.time(), 0
    for t in targets:
        blast[t] = blast_family(t, pools[t], tgenus[t])
        done += len(pools[t])
        print(f"  {t}: {len(pools[t])} BLASTed  [{done}/{total}, {(time.time() - t0) / 60:.1f} min]", flush=True)

    (OUT / "selected").mkdir(exist_ok=True)
    all_sel, summary = [], []
    print("Selecting …", flush=True)
    for t in targets:
        n_t = int(tables[t].target_seqs.iloc[0])
        df = tables[t][tables[t].kmer.isin(blast[t])].copy()
        for col in ["blast_total_hits", "blast_target_hits", "blast_ambiguous_hits", "blast_offtarget_hits",
                    "blast_weighted_offtarget", "blast_specificity",
                    "blast_weighted_specificity", "blast_top_offtarget", "blast_capped"]:
            df[col] = df.kmer.map(lambda k, c=col: blast[t][k][c])
        df.sort_values(["blast_weighted_specificity", "conservation_pct"], ascending=False).to_csv(
            OUT / "find_runs" / f"{t.replace(' ', '_')}_blast_scored.tsv", sep="\t", index=False)

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
        sel.to_csv(OUT / "selected" / f"{t.replace(' ', '_')}_selected_probes.csv", index=False)
        all_sel.append(sel)
        top_ot = Counter(x for x in sel.blast_top_offtarget if isinstance(x, str) and x)
        summary.append({
            "taxa": t, "silva_genus": tgenus[t], "target_seqs": n_t,
            "ambiguous_seqs_excluded": n_ambig[t],
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
            "top_offtarget_species": "; ".join(f"{n} ({c})" for n, c in top_ot.most_common(3)),
        })
        print(f"  {t}: floor {floor} -> {len(sel)} probes, pooled {100 * pooled:.1f}%", flush=True)

    pd.concat(all_sel).to_csv(OUT / "Species_all_selected_probes.csv", index=False)
    s = pd.DataFrame(summary)
    s.to_csv(OUT / "Species_summary.csv", index=False)
    print(s.drop(columns=["top_offtarget_species"]).to_string(index=False))


if __name__ == "__main__":
    main()
