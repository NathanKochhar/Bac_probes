"""BLAST-based near-match off-target specificity scoring for candidate k-mers."""

import subprocess
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

from .database import parse_taxonomy


def _kmers_to_fasta(kmers: list[str]) -> str:
    return "\n".join(f">kmer_{i}\n{kmer}" for i, kmer in enumerate(kmers))


def run_blast(
    kmers: list[str],
    blast_db: str | Path,
    threads: int = 4,
    perc_identity: float = 85.0,
    max_target_seqs: int = 1000,
) -> str:
    """
    BLAST a list of k-mers against a SILVA BLAST database.

    Uses blastn-short (tuned for sequences < 50 bp) and tabular output format.
    perc_identity=85 allows ~5 mismatches in a 32-mer, catching near-cross-binders.

    Returns raw tabular output (format 6 with custom columns).
    """
    query_fasta = _kmers_to_fasta(kmers)
    blast_db = str(blast_db)

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as qf:
        qf.write(query_fasta)
        query_path = qf.name

    try:
        result = subprocess.run(
            [
                "blastn",
                "-task", "blastn-short",
                "-query", query_path,
                "-db", blast_db,
                "-outfmt", "6 qseqid stitle pident length qlen",
                "-perc_identity", str(perc_identity),
                "-qcov_hsp_perc", "80",   # >=80% of the k-mer must align
                "-dust", "no",            # disable low-complexity masking for short seqs
                "-num_threads", str(threads),
                "-max_target_seqs", str(max_target_seqs),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError:
        raise RuntimeError(
            "blastn not found. Install BLAST+:\n"
            "  conda install -c bioconda blast"
        ) from None
    finally:
        Path(query_path).unlink(missing_ok=True)

    return result.stdout


def parse_blast_results(
    blast_output: str,
    kmers: list[str],
    target_name: str,
    level: str,
) -> dict[str, dict]:
    """
    Parse tabular BLAST output and classify every hit as on-target or off-target.

    SILVA stitle format: "<accession> Domain;Phylum;Class;..."
    We parse taxonomy directly from the hit title so no external taxonomy map is needed.

    Returns a dict:
        {
          kmer: {
            "blast_total_hits": int,
            "blast_target_hits": int,
            "blast_offtarget_hits": int,
            "blast_specificity": float,   # target_hits / total_hits
            "blast_top_offtarget": str,   # most common off-target taxon at this level
          }
        }
    """
    kmer_index = {f"kmer_{i}": kmer for i, kmer in enumerate(kmers)}
    target_name_lower = target_name.lower()

    # Accumulate counts per kmer
    data: dict[str, dict] = defaultdict(
        lambda: {"target": 0, "offtarget": 0, "offtarget_taxa": Counter()}
    )

    for line in blast_output.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 5:
            continue

        qseqid, stitle, pident, length, qlen = parts[:5]
        kmer = kmer_index.get(qseqid)
        if kmer is None:
            continue

        # Extract taxonomy from the hit title
        # stitle: "AB016480.1.1455 Bacteria;Firmicutes;Clostridia;..."
        if " " in stitle:
            tax_str = stitle.split(" ", 1)[1]
        else:
            tax_str = stitle

        taxonomy = parse_taxonomy(tax_str)
        hit_taxon = taxonomy.get(level, "")

        if hit_taxon.lower() == target_name_lower:
            data[kmer]["target"] += 1
        else:
            data[kmer]["offtarget"] += 1
            if hit_taxon:
                data[kmer]["offtarget_taxa"][hit_taxon] += 1

    results: dict[str, dict] = {}
    for kmer in kmers:
        entry = data.get(kmer, {"target": 0, "offtarget": 0, "offtarget_taxa": Counter()})
        total = entry["target"] + entry["offtarget"]
        specificity = entry["target"] / total if total > 0 else 1.0
        top_ot_list = entry["offtarget_taxa"].most_common(1)
        top_ot = top_ot_list[0][0] if top_ot_list else ""

        results[kmer] = {
            "blast_total_hits": total,
            "blast_target_hits": entry["target"],
            "blast_offtarget_hits": entry["offtarget"],
            "blast_specificity": round(specificity, 6),
            "blast_top_offtarget": top_ot,
        }

    return results
