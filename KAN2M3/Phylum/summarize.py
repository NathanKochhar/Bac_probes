"""Select ~10 highly specific, non-overlapping probes per phylum from the
`bac-probes find --no-blast` outputs and write per-phylum + summary CSVs.

Selection per phylum:
  1. Compute on-target precision = target hits / (target hits + exact off-target
     hits), i.e. the fraction of SILVA sequences carrying the k-mer that belong
     to the phylum. Raw off-target counts are not comparable between a phylum of
     1.8k and one of 136k sequences; precision is. Keep precision >= MIN_PREC.
  2. Try conservation floors from 50% downward; use the highest floor that
     yields N_PROBES non-overlapping probes (falls back to the lowest floor).
  3. Within the floor, rank by precision, then highest conservation.
  4. Skip any k-mer sharing a 16-mer with an already-picked probe (same locus).
Pooled coverage (target seqs hit by >=1 selected probe) is computed with one
extra SILVA pass.
"""
import sys
from pathlib import Path

import pandas as pd

from bac_probes.database import iter_silva

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
SILVA = ROOT / "silva_db" / "SILVA_138.1_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz"
INPUT_CSV = OUT.parent / "KAN2M3-phylum.csv"

N_PROBES = 10
MIN_PREC = float(sys.argv[1]) if len(sys.argv) > 1 else 0.95
FLOORS = [50, 45, 40, 35, 30, 25]
OVERLAP_WORD = 16

# CSV name -> SILVA 138.1 phylum name
SILVA_NAME = {
    "Bacillota-Firmicutes": "Firmicutes",
    "Pseudomonadota-Proteobacteria": "Proteobacteria",
    "Bacteroidota-Bacteroidota": "Bacteroidota",
    "Actinomycetota-Actinobacteriota": "Actinobacteriota",
    "Fusobacteriota-Fusobacteriota": "Fusobacteriota",
    "Campylobacterota-Campilobacterota": "Campylobacterota",
}


def words(seq: str) -> set[str]:
    return {seq[i:i + OVERLAP_WORD] for i in range(len(seq) - OVERLAP_WORD + 1)}


def pick(df: pd.DataFrame, floor: float) -> pd.DataFrame:
    pool = df[(df.precision >= MIN_PREC) & (df.conservation_pct >= floor)]
    pool = pool.sort_values(["precision", "conservation_pct"], ascending=[False, False])
    chosen, seen = [], set()
    for row in pool.itertuples():
        w = words(row.kmer)
        if w & seen:
            continue
        chosen.append(row.Index)
        seen |= w
        if len(chosen) == N_PROBES:
            break
    return df.loc[chosen]


def main() -> None:
    taxa = pd.read_csv(INPUT_CSV)["Taxa"].str.strip().tolist()
    selections, meta = {}, {}
    for taxon in taxa:
        name = SILVA_NAME[taxon]
        df = pd.read_csv(OUT / "find_runs" / f"{name}_find.tsv", sep="\t")
        target_hits = df.conservation_pct / 100 * df.target_seqs
        df["target_hits"] = target_hits.round().astype(int)
        df["precision"] = (target_hits / (target_hits + df.exact_offtarget)).round(4)
        for floor in FLOORS:
            sel = pick(df, floor)
            if len(sel) >= N_PROBES:
                break
        selections[name] = sel
        meta[name] = {"taxon": taxon, "floor": floor, "n_candidates": len(df),
                      "n_specific": int((df.precision >= MIN_PREC).sum())}
        print(f"{name}: floor {floor}% -> {len(sel)} probes")

    # Pooled coverage: one SILVA pass over all target phyla
    kmers = {n: list(s.kmer) for n, s in selections.items()}
    covered = {n: 0 for n in kmers}
    totals = {n: 0 for n in kmers}
    for _acc, tax, seq in iter_silva(SILVA):
        n = tax.get("phylum", "")
        if n in kmers:
            totals[n] += 1
            if any(k in seq for k in kmers[n]):
                covered[n] += 1

    all_rows, summary = [], []
    for name, sel in selections.items():
        m = meta[name]
        sel = sel.copy()
        sel.insert(0, "probe_rank", range(1, len(sel) + 1))
        sel.insert(0, "silva_phylum", name)
        sel.insert(0, "taxa", m["taxon"])
        sel.to_csv(OUT / "selected" / f"{name}_selected_probes.csv", index=False)
        all_rows.append(sel)
        summary.append({
            "taxa": m["taxon"],
            "silva_phylum": name,
            "target_seqs": totals[name],
            "candidates_found": m["n_candidates"],
            "candidates_precision_ge_%s" % MIN_PREC: m["n_specific"],
            "conservation_floor_used_pct": m["floor"],
            "probes_selected": len(sel),
            "conservation_min_pct": sel.conservation_pct.min(),
            "conservation_mean_pct": round(sel.conservation_pct.mean(), 2),
            "conservation_max_pct": sel.conservation_pct.max(),
            "pooled_coverage_pct": round(100 * covered[name] / totals[name], 2) if totals[name] else 0,
            "exact_offtarget_min": sel.exact_offtarget.min(),
            "exact_offtarget_median": sel.exact_offtarget.median(),
            "exact_offtarget_max": sel.exact_offtarget.max(),
            "precision_min": sel.precision.min(),
            "precision_mean": round(sel.precision.mean(), 4),
            "exact_specificity_min": sel.exact_specificity.min(),
            "best_probe": sel.probe.iloc[0] if len(sel) else "",
        })

    pd.concat(all_rows).to_csv(OUT / "Phylum_all_selected_probes.csv", index=False)
    s = pd.DataFrame(summary)
    s.to_csv(OUT / "Phylum_summary.csv", index=False)
    print(s.to_string(index=False))


if __name__ == "__main__":
    main()
