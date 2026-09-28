"""BLAST-based near-match off-target specificity scoring for candidate k-mers."""

import math
import subprocess
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

from .database import find_blast, parse_taxonomy


def _kmers_to_fasta(kmers: list[str]) -> str:
    return "\n".join(f">kmer_{i}\n{kmer}" for i, kmer in enumerate(kmers))


# ── Position-weighted mismatch scoring ────────────────────────────────────────

def mismatch_weight(pos: int, k: int) -> float:
    """
    Destabilisation weight for a mismatch at 0-indexed query position `pos`
    in a probe of length `k`.

    Uses a sine-squared curve so weight → 0 at both ends and peaks at 1.0 in
    the exact centre. This reflects the empirical observation that central
    mismatches are far more disruptive to probe hybridisation than terminal
    ones — relevant for Xenium in-situ probes where ligation/extension
    requires a well-formed central duplex.

    Examples for k=32:
        pos  0 (end)    → weight ≈ 0.002
        pos  8 (¼)      → weight ≈ 0.50
        pos 15 (centre) → weight ≈ 1.000
        pos 31 (end)    → weight ≈ 0.002
    """
    return math.sin(math.pi * (pos + 0.5) / k) ** 2


def _parse_btop(btop: str, qstart: int = 1) -> list[int]:
    """
    Parse a BLAST BTOP string into 0-indexed query positions of every mismatch
    or gap in the alignment.

    BTOP positions are relative to the start of the alignment, so `qstart`
    (BLAST's 1-based query start) is used to map them back onto the probe.

    BTOP format:
        <int>   – N consecutive matches (advance both query and subject)
        XY      – mismatch: query has X, subject has Y (advance both by 1)
        -Y      – gap in query / insertion in subject (subject advances, query does
                  not); recorded at the query base following the insertion
        X-      – gap in subject / deletion in subject (query advances by 1)
    """
    positions: list[int] = []
    query_pos = qstart - 1
    i = 0
    while i < len(btop):
        c = btop[i]
        if c.isdigit():
            j = i + 1
            while j < len(btop) and btop[j].isdigit():
                j += 1
            query_pos += int(btop[i:j])
            i = j
        elif c == "-":
            # Gap in query: subject has an extra base, query position doesn't advance.
            # The bulge still destabilises the duplex, so penalise it at this position.
            positions.append(query_pos)
            i += 2
        else:
            # Mismatch (XY) or gap in subject (X-): query advances by 1.
            positions.append(query_pos)
            query_pos += 1
            i += 2
    return positions


def _hit_mismatch_positions(btop: str, qstart: int, qend: int, qlen: int) -> list[int]:
    """
    0-indexed probe positions that do not pair with the subject in a BLAST hit:
    mismatches and gaps inside the alignment, plus any probe bases BLAST left
    unaligned at either end (positions before qstart or after qend).
    """
    overhang = list(range(qstart - 1)) + list(range(qend, qlen))
    return overhang + _parse_btop(btop, qstart)


def effective_binding_score(mismatch_positions: list[int], k: int) -> float:
    """
    Probe specificity score against a single off-target hit, given where its
    mismatches fall relative to the probe centre.

    Returns a value in [0, 1] where **1 = specific** (off-target won't bind)
    and **0 = non-specific** (off-target will bind):
        0.0  – perfect match  → certain off-target binding (no protection)
        ~1.0 – one central mismatch → off-target unlikely to hybridise (well protected)
        ~0.002 – one terminal mismatch → off-target may still bind (minimal protection)

    The score is the sum of position weights clamped to [0, 1]. A single
    perfectly-central mismatch (weight = 1.0) fully protects the probe, while
    a terminal mismatch (weight ≈ 0.002) provides almost no protection.
    """
    if not mismatch_positions:
        return 0.0
    penalty = sum(mismatch_weight(p, k) for p in mismatch_positions)
    return min(1.0, penalty)


# ── BLAST interface ───────────────────────────────────────────────────────────

def run_blast(
    kmers: list[str],
    blast_db: str | Path,
    threads: int = 4,
    perc_identity: float = 85.0,
    max_target_seqs: int = 1000,
) -> str:
    """
    BLAST a list of k-mers against a SILVA BLAST database.

    Uses blastn-short (optimised for sequences < 50 bp). The BTOP field is
    requested so that mismatch positions can be extracted for position-weighted
    off-target scoring.

    Returns raw tabular output
    (format 6: qseqid stitle pident length qlen btop qstart qend).
    """
    query_fasta = _kmers_to_fasta(kmers)

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as qf:
        qf.write(query_fasta)
        query_path = qf.name

    try:
        result = subprocess.run(
            [
                find_blast("blastn"),
                "-task", "blastn-short",
                "-query", query_path,
                "-db", str(blast_db),
                "-outfmt", "6 qseqid stitle pident length qlen btop qstart qend",
                "-perc_identity", str(perc_identity),
                "-qcov_hsp_perc", "80",
                "-dust", "no",
                "-num_threads", str(threads),
                "-max_target_seqs", str(max_target_seqs),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError:
        raise RuntimeError(find_blast("blastn")) from None  # surfaces helpful install message
    finally:
        Path(query_path).unlink(missing_ok=True)

    return result.stdout


# ── Result parsing ────────────────────────────────────────────────────────────

def parse_blast_results(
    blast_output: str,
    kmers: list[str],
    target_name: str,
    level: str,
    max_target_seqs: int = 1000,
) -> dict[str, dict]:
    """
    Parse tabular BLAST output and classify every hit as on-target or off-target.

    Off-target hits are scored in two ways:
      • Raw count  – every off-target hit counts as 1 (standard specificity).
      • Weighted   – each off-target hit is multiplied by its effective_binding_score,
                     which is low when mismatches fall near the probe centre.  This
                     makes the weighted score more lenient toward near-matches that are
                     unlikely to actually hybridise (important for Xenium probes).

    Returns per-kmer dict with keys:
        blast_total_hits          – raw number of BLAST hits (target + off-target)
        blast_target_hits         – hits within the target taxon
        blast_offtarget_hits      – raw count of off-target hits
        blast_weighted_offtarget  – sum of effective binding scores for off-target hits
        blast_specificity         – target_hits / total_hits  (raw)
        blast_weighted_specificity – target_hits / (target_hits + weighted_offtarget)
        blast_top_offtarget       – most frequent off-target taxon name
    """
    kmer_index = {f"kmer_{i}": kmer for i, kmer in enumerate(kmers)}
    target_name_lower = target_name.lower()

    data: dict[str, dict] = defaultdict(
        lambda: {
            "target": 0,
            "offtarget_raw": 0,
            "offtarget_weighted": 0.0,
            "offtarget_taxa": Counter(),
        }
    )

    for line in blast_output.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 6:
            continue

        qseqid, stitle, pident, length, qlen, btop = parts[:6]
        kmer = kmer_index.get(qseqid)
        if kmer is None:
            continue

        # stitle from a -parse_seqids BLAST db is the sequence description
        # *after* the accession (i.e. the taxonomy string directly).
        # Do not split on spaces — species names contain spaces.
        hit_taxon = parse_taxonomy(stitle).get(level, "")

        if hit_taxon.lower() == target_name_lower:
            data[kmer]["target"] += 1
        else:
            k = len(kmer)
            if len(parts) >= 8:
                mismatch_pos = _hit_mismatch_positions(btop, int(parts[6]), int(parts[7]), k)
            else:
                mismatch_pos = _parse_btop(btop)
            specificity = effective_binding_score(mismatch_pos, k)
            data[kmer]["offtarget_raw"] += 1
            data[kmer]["offtarget_weighted"] += (1.0 - specificity)  # binding concern = 1 - specificity
            if hit_taxon:
                data[kmer]["offtarget_taxa"][hit_taxon] += 1

    results: dict[str, dict] = {}
    for kmer in kmers:
        entry = data.get(
            kmer,
            {"target": 0, "offtarget_raw": 0, "offtarget_weighted": 0.0, "offtarget_taxa": Counter()},
        )
        total_raw = entry["target"] + entry["offtarget_raw"]
        total_weighted = entry["target"] + entry["offtarget_weighted"]

        spec_raw = entry["target"] / total_raw if total_raw > 0 else 1.0
        spec_weighted = entry["target"] / total_weighted if total_weighted > 0 else 1.0

        top_ot_list = entry["offtarget_taxa"].most_common(1)

        results[kmer] = {
            "blast_total_hits": total_raw,
            "blast_target_hits": entry["target"],
            "blast_offtarget_hits": entry["offtarget_raw"],
            "blast_weighted_offtarget": round(entry["offtarget_weighted"], 3),
            "blast_specificity": round(spec_raw, 6),
            "blast_weighted_specificity": round(spec_weighted, 6),
            "blast_top_offtarget": top_ot_list[0][0] if top_ot_list else "",
            # True when BLAST returned exactly max_target_seqs hits — specificity
            # estimates are unreliable because the result set is truncated.
            "blast_capped": total_raw >= max_target_seqs,
        }

    return results
