"""SILVA 16S rRNA database download, parsing, and BLAST database construction."""

import gzip
import shutil
import subprocess
import urllib.request
from pathlib import Path
from typing import Iterator

from Bio import SeqIO

SILVA_URL = (
    "https://ftp.arb-silva.de/release_138_1/Exports/"
    "SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz"
)
SILVA_FILENAME = "SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz"

# Canonical SILVA taxonomy depth; entries with fewer levels get empty strings.
TAXONOMY_LEVELS = ["domain", "phylum", "class", "order", "family", "genus", "species"]


def parse_taxonomy(taxonomy_str: str) -> dict[str, str]:
    """Parse a SILVA semicolon-delimited taxonomy string into level→name dict."""
    taxa = [t.strip() for t in taxonomy_str.split(";")]
    return {level: (taxa[i] if i < len(taxa) else "") for i, level in enumerate(TAXONOMY_LEVELS)}


def iter_silva(fasta_path: str | Path) -> Iterator[tuple[str, dict[str, str], str]]:
    """
    Yield (accession, taxonomy_dict, sequence) for every entry in a SILVA FASTA.

    Handles gzipped (.gz) and plain FASTA files.
    RNA U→T substitution is applied so all sequences are DNA.

    SILVA header format:
        >AB016480.1.1455 Bacteria;Firmicutes;Clostridia;...
    """
    fasta_path = Path(fasta_path)
    opener = gzip.open if fasta_path.suffix == ".gz" else open

    with opener(fasta_path, "rt") as fh:
        for record in SeqIO.parse(fh, "fasta"):
            # record.description is the full header (id + rest); strip the id
            desc = record.description
            taxonomy_str = desc.split(" ", 1)[1] if " " in desc else ""
            taxonomy = parse_taxonomy(taxonomy_str)
            seq = str(record.seq).upper().replace("U", "T")
            yield record.id, taxonomy, seq


def download_silva(output_dir: str | Path = ".") -> Path:
    """Download the SILVA 138.1 NR99 SSU FASTA (~1.5 GB compressed)."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    dest = output_dir / SILVA_FILENAME

    if dest.exists():
        print(f"SILVA already present: {dest}")
        return dest

    print(f"Downloading SILVA 138.1 NR99 from:\n  {SILVA_URL}")
    print("This download is ~1.5 GB and may take several minutes.")

    def _progress(count: int, block: int, total: int) -> None:
        mb_done = count * block / 1e6
        mb_total = total / 1e6
        print(f"\r  {mb_done:.0f} / {mb_total:.0f} MB", end="", flush=True)

    urllib.request.urlretrieve(SILVA_URL, dest, reporthook=_progress)
    print(f"\nSaved to {dest}")
    return dest


def decompress_silva(gz_path: str | Path, out_path: str | Path | None = None) -> Path:
    """Decompress a gzipped SILVA FASTA; returns path to plain FASTA."""
    gz_path = Path(gz_path)
    if out_path is None:
        out_path = gz_path.with_suffix("")  # strip .gz
    out_path = Path(out_path)

    if out_path.exists():
        return out_path

    print(f"Decompressing {gz_path.name} …")
    with gzip.open(gz_path, "rb") as src, open(out_path, "wb") as dst:
        shutil.copyfileobj(src, dst)
    print(f"Decompressed to {out_path}")
    return out_path


def build_blast_db(fasta_path: str | Path, db_dir: str | Path) -> Path:
    """
    Build a BLAST nucleotide database from a SILVA FASTA.

    Decompresses the file if gzipped (makeblastdb needs a seekable file for
    -parse_seqids). The uncompressed FASTA is kept in db_dir for reuse.

    Returns the BLAST database path prefix (e.g. db_dir/silva).
    """
    fasta_path = Path(fasta_path)
    db_dir = Path(db_dir)
    db_dir.mkdir(parents=True, exist_ok=True)
    db_prefix = db_dir / "silva"

    if (db_dir / "silva.nhr").exists() or (db_dir / "silva.00.nhr").exists():
        print(f"BLAST database already exists at {db_prefix}")
        return db_prefix

    # makeblastdb needs a plain (non-gzipped) seekable file when using -parse_seqids
    if fasta_path.suffix == ".gz":
        plain = db_dir / fasta_path.stem  # e.g. silva.fasta
        fasta_path = decompress_silva(fasta_path, plain)

    print("Building BLAST database (this takes a few minutes) …")
    try:
        subprocess.run(
            [
                "makeblastdb",
                "-in", str(fasta_path),
                "-dbtype", "nucl",
                "-out", str(db_prefix),
                "-parse_seqids",
                "-title", "SILVA_138.1_NR99",
            ],
            check=True,
        )
    except FileNotFoundError:
        raise RuntimeError(
            "makeblastdb not found. Install BLAST+:\n"
            "  conda install -c bioconda blast"
        )
    print(f"BLAST database built at {db_prefix}")
    return db_prefix


def list_taxa(fasta_path: str | Path, level: str) -> list[tuple[str, int]]:
    """
    Return (taxon_name, sequence_count) pairs for all taxa at the given level,
    sorted by descending sequence count.
    """
    if level not in TAXONOMY_LEVELS:
        raise ValueError(f"Unknown level '{level}'. Choose from: {TAXONOMY_LEVELS}")

    counts: dict[str, int] = {}
    for _, taxonomy, _ in iter_silva(fasta_path):
        name = taxonomy.get(level, "")
        if name:
            counts[name] = counts.get(name, 0) + 1

    return sorted(counts.items(), key=lambda x: x[1], reverse=True)
