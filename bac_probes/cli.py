"""Command-line interface for bac-probes."""

import sys
from pathlib import Path

import click
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

from .database import (
    TAXONOMY_LEVELS,
    build_blast_db,
    download_silva,
    iter_silva,
    list_taxa,
)
from .kmers import gc_content, reverse_complement, score_conservation, score_offtarget_exact
from .specificity import parse_blast_results, run_blast


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option()
def cli() -> None:
    """bac-probes: find taxon-specific 16S rRNA k-mers from the SILVA database.

    \b
    Typical workflow:
      1. bac-probes download --output-dir ./silva_db
      2. bac-probes build-db ./silva_db/SILVA_*.fasta.gz --db-dir ./silva_db
      3. bac-probes find  ./silva_db/SILVA_*.fasta.gz  Streptococcus  --level genus
    """


# ─── download ──────────────────────────────────────────────────────────────────

@cli.command()
@click.option(
    "--output-dir", "-o", default="silva_db", show_default=True,
    help="Directory to save the SILVA FASTA file.",
)
def download(output_dir: str) -> None:
    """Download the SILVA 138.1 NR99 16S SSU reference database (~1.5 GB)."""
    download_silva(output_dir)


# ─── build-db ──────────────────────────────────────────────────────────────────

@cli.command("build-db")
@click.argument("silva_fasta")
@click.option(
    "--db-dir", default="silva_db", show_default=True,
    help="Directory where the BLAST database will be created.",
)
def build_db(silva_fasta: str, db_dir: str) -> None:
    """Build a BLAST database from SILVA_FASTA for off-target scoring.

    \b
    SILVA_FASTA  Path to the SILVA NR99 FASTA (gzipped or plain).
    """
    build_blast_db(silva_fasta, db_dir)
    click.echo("Done. Use '--blast-db silva_db/silva' with 'bac-probes find'.")


# ─── list-taxa ─────────────────────────────────────────────────────────────────

@cli.command("list-taxa")
@click.argument("silva_fasta")
@click.option(
    "--level", "-l",
    type=click.Choice(TAXONOMY_LEVELS, case_sensitive=False),
    default="genus", show_default=True,
    help="Taxonomic level to list.",
)
@click.option("--top", default=30, show_default=True, help="Number of top taxa to show.")
def list_taxa_cmd(silva_fasta: str, level: str, top: int) -> None:
    """List taxon names and sequence counts at a given taxonomic level.

    Useful for finding the exact spelling required by 'bac-probes find'.

    \b
    SILVA_FASTA  Path to the SILVA NR99 FASTA (gzipped or plain).
    """
    click.echo(f"Counting taxa at level '{level}' …")
    taxa = list_taxa(silva_fasta, level)
    click.echo(f"\n{'Taxon':<50} {'Sequences':>10}")
    click.echo("-" * 62)
    for name, count in taxa[:top]:
        click.echo(f"{name:<50} {count:>10,}")
    if len(taxa) > top:
        click.echo(f"… and {len(taxa) - top} more (pipe through 'head' or increase --top)")


# ─── find ──────────────────────────────────────────────────────────────────────

@cli.command()
@click.argument("silva_fasta")
@click.argument("target_name")
@click.option(
    "--level", "-l",
    type=click.Choice(TAXONOMY_LEVELS, case_sensitive=False),
    default="genus", show_default=True,
    help="Taxonomic level for the target name.",
)
@click.option(
    "--target-fasta", default=None,
    help="Plain FASTA of target sequences (e.g. downloaded from NCBI). "
         "When provided, target sequences come from this file instead of SILVA; "
         "SILVA is still used for background off-target scoring.",
)
@click.option("--kmer-size", "-k", default=32, show_default=True, help="K-mer length (bp).")
@click.option(
    "--min-conservation", default=0.80, show_default=True,
    help="Minimum fraction of target sequences containing the k-mer [0–1].",
)
@click.option("--gc-min", default=0.35, show_default=True, help="Minimum k-mer GC content.")
@click.option("--gc-max", default=0.65, show_default=True, help="Maximum k-mer GC content.")
@click.option(
    "--max-homopolymer", default=5, show_default=True,
    help="Reject k-mers with a homopolymer run >= this length.",
)
@click.option(
    "--blast/--no-blast", "use_blast", default=True, show_default=True,
    help="Run BLAST to score near-match off-target binding (requires build-db first).",
)
@click.option(
    "--blast-db", default=None,
    help="BLAST database prefix (default: silva_db/silva).",
)
@click.option(
    "--blast-identity", default=85.0, show_default=True,
    help="Minimum BLAST percent identity to count as an off-target hit.",
)
@click.option(
    "--top-n-blast", default=500, show_default=True,
    help="Submit at most this many top candidates to BLAST (ranked by conservation).",
)
@click.option("--threads", default=4, show_default=True, help="CPU threads for BLAST.")
@click.option(
    "--output", "-o", default="probes.tsv", show_default=True,
    help="Output TSV file path.",
)
def find(
    silva_fasta: str,
    target_name: str,
    level: str,
    target_fasta: str | None,
    kmer_size: int,
    min_conservation: float,
    gc_min: float,
    gc_max: float,
    max_homopolymer: int,
    use_blast: bool,
    blast_db: str | None,
    blast_identity: float,
    top_n_blast: int,
    threads: int,
    output: str,
) -> None:
    """Find taxon-specific k-mers for a bacterial clade.

    \b
    SILVA_FASTA  Path to the SILVA NR99 FASTA (gzipped or plain).
    TARGET_NAME  Taxonomic name to target (e.g. 'Firmicutes', 'Streptococcus').
                 Also used to label the target in BLAST output columns.

    \b
    Output columns:
      kmer                       – k-mer sequence (same strand as reference)
      probe                      – reverse complement (order this as the probe)
      gc_pct                     – GC content (%)
      conservation_pct           – % of target sequences containing this k-mer exactly
      target_seqs                – number of target sequences used
      exact_offtarget            – number of non-target SILVA seqs containing it exactly
      exact_bg_seqs              – total non-target sequences scanned
      exact_specificity          – 1 − (exact_offtarget / exact_bg_seqs)
      blast_total_hits           – total BLAST hits (if --blast)
      blast_target_hits          – BLAST hits within target taxon
      blast_offtarget_hits       – raw count of off-target BLAST hits
      blast_weighted_offtarget   – off-target hits weighted by mismatch position;
                                   central mismatches score near 0 (unlikely to bind),
                                   terminal mismatches score near 1 (may still bind)
      blast_specificity          – target_hits / total_hits (raw)
      blast_weighted_specificity – target_hits / (target_hits + weighted_offtarget);
                                   primary ranking metric for Xenium probe design
      blast_top_offtarget        – most frequent off-target taxon
    """
    silva_path = Path(silva_fasta)
    blast_db_path = Path(blast_db) if blast_db else Path("silva_db/silva")
    target_lower = target_name.lower()

    # ── Pass 1: collect target sequences ─────────────────────────────────────
    if target_fasta:
        # Load target sequences from a plain FASTA (e.g. NCBI download).
        # SILVA is used only for background off-target scoring.
        click.echo(f"\nLoading target sequences from {target_fasta} …")
        target_seqs = [
            str(r.seq).upper().replace("U", "T")
            for r in SeqIO.parse(target_fasta, "fasta")
        ]
        if not target_seqs:
            click.echo(f"ERROR: No sequences found in {target_fasta}.", err=True)
            sys.exit(1)
        click.echo(f"  Target sequences : {len(target_seqs):,}  (from {target_fasta})")
        click.echo(f"\nCounting background sequences in SILVA …")
        n_bg_pass1 = sum(
            1 for _acc, taxonomy, _seq in tqdm(iter_silva(silva_path), desc="  reading", unit=" seq")
            if taxonomy.get(level, "").lower() != target_lower
        )
        click.echo(f"  Background seqs  : {n_bg_pass1:,}")
    else:
        click.echo(f"\nPass 1/2 — scanning SILVA for '{target_name}' at level '{level}' …")
        target_seqs: list[str] = []
        n_bg_pass1 = 0
        for _acc, taxonomy, seq in tqdm(iter_silva(silva_path), desc="  reading", unit=" seq"):
            if taxonomy.get(level, "").lower() == target_lower:
                target_seqs.append(seq)
            else:
                n_bg_pass1 += 1

        if not target_seqs:
            click.echo(
                f"\nERROR: No sequences found for '{target_name}' at level '{level}'.\n"
                "  • Check spelling (run 'bac-probes list-taxa' to browse names).\n"
                "  • Try a higher level (phylum → class → order → family → genus → species).",
                err=True,
            )
            sys.exit(1)

        click.echo(
            f"  Target sequences : {len(target_seqs):,}\n"
            f"  Background seqs  : {n_bg_pass1:,}"
        )

    # ── Conservation scoring ──────────────────────────────────────────────────
    click.echo(f"\nExtracting conserved {kmer_size}-mers (conservation ≥ {min_conservation:.0%}) …")
    conservation, n_target = score_conservation(
        tqdm(target_seqs, desc="  k-mer scan", unit=" seq"),
        k=kmer_size,
        min_conservation=min_conservation,
        gc_range=(gc_min, gc_max),
        max_homopolymer=max_homopolymer,
    )

    if not conservation:
        click.echo(
            "\nNo k-mers passed the conservation threshold.\n"
            "  • Lower --min-conservation (e.g. 0.60).\n"
            "  • Widen --gc-min / --gc-max.\n"
            "  • Increase --max-homopolymer.",
            err=True,
        )
        sys.exit(1)

    click.echo(f"  Candidates after conservation filter: {len(conservation):,}")

    # ── Pass 2: exact off-target scoring (streaming background) ──────────────
    click.echo("\nPass 2/2 — exact off-target scoring across background sequences …")

    candidate_set = set(conservation.keys())

    def _bg_stream():
        for _acc, taxonomy, seq in iter_silva(silva_path):
            if taxonomy.get(level, "").lower() != target_lower:
                yield seq

    offtarget_counts, n_bg = score_offtarget_exact(
        candidate_set,
        tqdm(_bg_stream(), desc="  scanning", total=n_bg_pass1, unit=" seq"),
        k=kmer_size,
    )

    # ── Build results table ───────────────────────────────────────────────────
    rows = []
    for kmer, cons in conservation.items():
        ot_count = offtarget_counts.get(kmer, 0)
        ot_frac = ot_count / n_bg if n_bg > 0 else 0.0
        rows.append(
            {
                "kmer": kmer,
                "probe": reverse_complement(kmer),
                "gc_pct": round(gc_content(kmer) * 100, 1),
                "conservation_pct": round(cons * 100, 2),
                "target_seqs": n_target,
                "exact_offtarget": ot_count,
                "exact_bg_seqs": n_bg,
                "exact_specificity": round(1.0 - ot_frac, 6),
            }
        )

    df = pd.DataFrame(rows).sort_values(
        ["exact_specificity", "conservation_pct"], ascending=[False, False]
    )

    # ── BLAST near-match scoring ──────────────────────────────────────────────
    if use_blast:
        db_exists = (
            blast_db_path.with_suffix(".nhr").exists()
            or (blast_db_path.parent / (blast_db_path.name + ".00.nhr")).exists()
        )
        if not db_exists:
            click.echo(
                f"\nWARNING: BLAST database not found at '{blast_db_path}'.\n"
                "  Run 'bac-probes build-db <silva.fasta.gz>' first, or use --no-blast.",
                err=True,
            )
        else:
            n_blast = min(top_n_blast, len(df))
            click.echo(f"\nBLASTing top {n_blast} candidates (identity ≥ {blast_identity:.0f}%) …")
            top_kmers: list[str] = df.head(n_blast)["kmer"].tolist()

            blast_out = run_blast(
                top_kmers,
                blast_db_path,
                threads=threads,
                perc_identity=blast_identity,
            )
            blast_scores = parse_blast_results(blast_out, top_kmers, target_name, level)

            for col, key in [
                ("blast_total_hits",           "blast_total_hits"),
                ("blast_target_hits",          "blast_target_hits"),
                ("blast_offtarget_hits",       "blast_offtarget_hits"),
                ("blast_weighted_offtarget",   "blast_weighted_offtarget"),
                ("blast_specificity",          "blast_specificity"),
                ("blast_weighted_specificity", "blast_weighted_specificity"),
                ("blast_top_offtarget",        "blast_top_offtarget"),
            ]:
                default = "" if col == "blast_top_offtarget" else pd.NA
                df[col] = df["kmer"].map(lambda k, _k=key, _d=default: blast_scores.get(k, {}).get(_k, _d))

            # Re-sort: weighted specificity is the primary metric for Xenium probes
            blast_mask = df["blast_weighted_specificity"].notna()
            df_blasted = df[blast_mask].sort_values(
                ["blast_weighted_specificity", "conservation_pct"], ascending=[False, False]
            )
            df_rest = df[~blast_mask]
            df = pd.concat([df_blasted, df_rest], ignore_index=True)

    # ── Write output ──────────────────────────────────────────────────────────
    df.to_csv(output, sep="\t", index=False)

    click.echo(f"\nResults written to: {output}")
    click.echo(f"Total k-mer candidates: {len(df):,}")

    display_cols = ["kmer", "gc_pct", "conservation_pct", "exact_specificity"]
    if "blast_weighted_specificity" in df.columns:
        display_cols += ["blast_specificity", "blast_weighted_specificity"]

    click.echo("\nTop 10 candidates:")
    click.echo(df[display_cols].head(10).to_string(index=False))
