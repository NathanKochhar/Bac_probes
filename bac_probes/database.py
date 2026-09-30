"""SILVA 16S rRNA database download, parsing, and BLAST database construction."""

import glob
import gzip
import shutil
import subprocess
import urllib.request
from pathlib import Path
from typing import Iterator

from Bio import SeqIO

# Full-length alignment, truncated version — same taxonomy headers as the
# unaligned release but sequences contain '.' / '-' gap characters that are
# stripped on read.  Using this file keeps the tool aligned with the standard
# SILVA release used by most 16S pipelines.
SILVA_URL = (
    "https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/"
    "SILVA_138.1_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz"
)
SILVA_FILENAME = "SILVA_138.1_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz"

# Canonical SILVA taxonomy depth; entries with fewer levels get empty strings.
TAXONOMY_LEVELS = ["domain", "phylum", "class", "order", "family", "genus", "species"]


# ── BLAST binary discovery ────────────────────────────────────────────────────

def find_blast(name: str) -> str:
    """
    Locate a BLAST+ binary (e.g. 'blastn', 'makeblastdb') without requiring
    the user to manually update PATH.

    Search order:
      1. Anything already on the current PATH  (shutil.which)
      2. /tmp/ncbi-blast-*/bin/  (manual NCBI tarballs, newest version first)
      3. /usr/local/bin, /opt/homebrew/bin, /opt/local/bin  (package managers)
    """
    import shutil as _shutil
    found = _shutil.which(name)
    if found:
        return found

    candidates: list[str] = sorted(
        glob.glob(f"/tmp/ncbi-blast-*/bin/{name}"), reverse=True
    )
    for p in candidates + [
        f"/usr/local/bin/{name}",
        f"/opt/homebrew/bin/{name}",
        f"/opt/local/bin/{name}",
    ]:
        if Path(p).exists():
            return p

    raise RuntimeError(
        f"{name} not found. Install BLAST+:\n"
        "  conda install -c bioconda blast\n"
        "  brew install blast\n"
        "  # or download: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
    )


# ── Taxonomy helpers ──────────────────────────────────────────────────────────

def parse_taxonomy(taxonomy_str: str) -> dict[str, str]:
    """Parse a SILVA semicolon-delimited taxonomy string into level→name dict."""
    taxa = [t.strip() for t in taxonomy_str.split(";")]
    return {level: (taxa[i] if i < len(taxa) else "") for i, level in enumerate(TAXONOMY_LEVELS)}


def taxon_matches(name: str, target: str, level: str) -> bool:
    """
    True if a SILVA taxon name belongs to the target taxon (case-insensitive).

    At species level SILVA labels many sequences with a strain or subspecies
    ("Helicobacter pylori 26695", "Campylobacter jejuni subsp. doylei"), so a
    species target also matches any name that starts with it followed by a
    space. Other levels require an exact match, so e.g. the genus "Clostridium"
    does not swallow "Clostridium sensu stricto 1".
    """
    name, target = name.strip().lower(), target.strip().lower()
    if not target:
        return False
    if name == target:
        return True
    return level == "species" and name.startswith(target + " ")


_UNNAMED_FIRST_WORDS = {"uncultured", "unidentified", "unclassified", "metagenome", "bacterium"}
_UNNAMED_EPITHETS = {"sp.", "sp", "bacterium", "cf.", "aff.", "oral", "genomosp."}


def is_unnamed_species(name: str) -> bool:
    """
    True if a SILVA species field gives no real species name, e.g.
    "uncultured bacterium", "Fusobacterium sp.", "uncultured Helicobacter sp.",
    "Fusobacterium sp. oral taxon 203". Such sequences may belong to any
    species of their SILVA genus, including the target.
    """
    w = name.split()
    if len(w) < 2:
        return True
    return (
        w[0].lower() in _UNNAMED_FIRST_WORDS
        or not w[0][:1].isupper()
        or not w[1][:1].islower()
        or w[1] in _UNNAMED_EPITHETS
    )


# ── FASTA iteration ───────────────────────────────────────────────────────────

def iter_silva(fasta_path: str | Path) -> Iterator[tuple[str, dict[str, str], str]]:
    """
    Yield (accession, taxonomy_dict, sequence) for every entry in a SILVA FASTA.

    Works with both the unaligned release and the full-alignment/truncated
    release (SILVA_*_full_align_trunc.fasta.gz).  Gap characters ('.' terminal
    gaps, '-' internal gaps) are stripped so the returned sequence is always
    plain DNA.  U→T substitution is also applied.

    SILVA header format:
        >AB016480.1.1455 Bacteria;Firmicutes;Clostridia;...
    """
    fasta_path = Path(fasta_path)
    opener = gzip.open if fasta_path.suffix == ".gz" else open

    with opener(fasta_path, "rt") as fh:
        for record in SeqIO.parse(fh, "fasta"):
            desc = record.description
            taxonomy_str = desc.split(" ", 1)[1] if " " in desc else ""
            taxonomy = parse_taxonomy(taxonomy_str)
            seq = str(record.seq).upper().replace("U", "T").replace(".", "").replace("-", "")
            if seq:
                yield record.id, taxonomy, seq


# ── Download ──────────────────────────────────────────────────────────────────

def download_silva(output_dir: str | Path = ".") -> Path:
    """Download the SILVA 138.1 NR99 full-alignment truncated FASTA (~1 GB compressed)."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    dest = output_dir / SILVA_FILENAME

    if dest.exists():
        print(f"SILVA already present: {dest}")
        return dest

    print(f"Downloading SILVA 138.1 NR99 (full alignment) from:\n  {SILVA_URL}")
    print("This download is ~1 GB and may take several minutes.")

    def _progress(count: int, block: int, total: int) -> None:
        mb_done = count * block / 1e6
        mb_total = total / 1e6
        print(f"\r  {mb_done:.0f} / {mb_total:.0f} MB", end="", flush=True)

    urllib.request.urlretrieve(SILVA_URL, dest, reporthook=_progress)
    print(f"\nSaved to {dest}")
    return dest


# ── BLAST database construction ───────────────────────────────────────────────

def _write_degapped_fasta(src_path: Path, out_path: Path) -> None:
    """
    Stream a (possibly gzipped) SILVA FASTA and write a gap-stripped plain
    FASTA suitable for makeblastdb.

    The aligned SILVA release contains '.' and '-' characters that would cause
    makeblastdb to reject or mishandle sequences.  Gaps are removed here so the
    BLAST database contains only biological sequence.
    """
    opener = gzip.open if src_path.suffix == ".gz" else open
    written = 0
    print(f"Writing gap-stripped FASTA to {out_path.name} …")
    with opener(src_path, "rt") as src, open(out_path, "w") as dst:
        for record in SeqIO.parse(src, "fasta"):
            seq = str(record.seq).upper().replace("U", "T").replace(".", "").replace("-", "")
            if seq:
                dst.write(f">{record.description}\n{seq}\n")
                written += 1
    print(f"  Wrote {written:,} sequences")


def build_blast_db(fasta_path: str | Path, db_dir: str | Path) -> Path:
    """
    Build a BLAST nucleotide database from a SILVA FASTA.

    For the aligned release, gaps are stripped before indexing so that BLAST
    operates on plain biological sequence.  The gap-stripped FASTA is kept in
    db_dir as 'silva_nogap.fasta' for reuse on subsequent runs.

    Returns the BLAST database path prefix (db_dir/silva).
    """
    fasta_path = Path(fasta_path)
    db_dir = Path(db_dir)
    db_dir.mkdir(parents=True, exist_ok=True)
    db_prefix = db_dir / "silva"

    if (db_dir / "silva.nhr").exists() or (db_dir / "silva.00.nhr").exists():
        print(f"BLAST database already exists at {db_prefix}")
        return db_prefix

    # Always write a fresh gap-stripped FASTA (handles both aligned and
    # unaligned releases; stripping from an unaligned file is a no-op).
    nogap = db_dir / "silva_nogap.fasta"
    if not nogap.exists():
        _write_degapped_fasta(fasta_path, nogap)

    makeblastdb = find_blast("makeblastdb")
    print("Building BLAST database (this takes a few minutes) …")
    subprocess.run(
        [
            makeblastdb,
            "-in", str(nogap),
            "-dbtype", "nucl",
            "-out", str(db_prefix),
            "-parse_seqids",
            "-title", "SILVA_138.1_NR99",
        ],
        check=True,
    )
    print(f"BLAST database built at {db_prefix}")
    return db_prefix


def list_taxa(fasta_path: str | Path, level: str) -> list[tuple[str, int]]:
    """
    Return (taxon_name, sequence_count) pairs for all taxa at the given level,
    sorted by descending sequence count.
    """
    from tqdm import tqdm

    if level not in TAXONOMY_LEVELS:
        raise ValueError(f"Unknown level '{level}'. Choose from: {TAXONOMY_LEVELS}")

    counts: dict[str, int] = {}
    for _, taxonomy, _ in tqdm(iter_silva(fasta_path), desc="  scanning", unit=" seq"):
        name = taxonomy.get(level, "")
        if name:
            counts[name] = counts.get(name, 0) + 1

    return sorted(counts.items(), key=lambda x: x[1], reverse=True)
