"""Command-line interface for bac-probes."""

import sys
from collections import Counter
from pathlib import Path

import click
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

from .database import (
    TAXONOMY_LEVELS,
    is_unnamed_species,
    taxon_matches,
    build_blast_db,
    download_silva,
    iter_silva,
    list_taxa,
)
from .kmers import gc_content, greedy_cocktail, merge_overlapping_kmers, rescore_merged_conservation, reverse_complement, score_conservation, score_conservation_with_coverage, score_offtarget_exact
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
    # Output uses non-ASCII characters (≥, …, —). On Windows, piped/redirected
    # output defaults to a legacy code page that can't encode them.
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure") and (stream.encoding or "").lower() != "utf-8":
            stream.reconfigure(encoding="utf-8")


# ─── download ──────────────────────────────────────────────────────────────────

@cli.command()
@click.option(
    "--output-dir", "-o", default="silva_db", show_default=True,
    help="Directory to save the SILVA FASTA file.",
)
def download(output_dir: str) -> None:
    """Download the SILVA 138.1 NR99 16S SSU reference database (~1 GB compressed)."""
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
@click.option(
    "--blast-max-target-seqs", default=1000, show_default=True,
    help="BLAST max_target_seqs — increase for large target taxa where cap inflates specificity.",
)
@click.option(
    "--merge-overlapping/--no-merge-overlapping", default=False, show_default=True,
    help="Merge overlapping k-mers into longer contiguous probe sequences.",
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
    blast_max_target_seqs: int,
    merge_overlapping: bool,
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
      blast_weighted_offtarget   – off-target hits weighted by binding concern;
                                   central mismatches contribute ~0 (unlikely to bind),
                                   terminal mismatches contribute ~1 (may still bind)
      blast_specificity          – target_hits / total_hits (raw)
      blast_weighted_specificity – target_hits / (target_hits + weighted_offtarget);
                                   primary ranking metric for Xenium probe design
      blast_top_offtarget        – most frequent off-target taxon
    """
    silva_path = Path(silva_fasta)
    blast_db_path = Path(blast_db) if blast_db else Path("silva_db/silva")

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
            if not taxon_matches(taxonomy.get(level, ""), target_name, level)
        )
        click.echo(f"  Background seqs  : {n_bg_pass1:,}")
    else:
        click.echo(f"\nPass 1/2 — scanning SILVA for '{target_name}' at level '{level}' …")
        target_seqs: list[str] = []
        n_bg_pass1 = 0
        for _acc, taxonomy, seq in tqdm(iter_silva(silva_path), desc="  reading", unit=" seq"):
            if taxon_matches(taxonomy.get(level, ""), target_name, level):
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

    # ── Merge overlapping k-mers into longer probes ───────────────────────────
    if merge_overlapping:
        click.echo(f"\nMerging overlapping {kmer_size}-mers …")
        merged_seqs = merge_overlapping_kmers(set(conservation.keys()), kmer_size)
        conservation = rescore_merged_conservation(
            target_seqs, merged_seqs, min_conservation,
            gc_range=(gc_min, gc_max), max_homopolymer=max_homopolymer,
        )
        click.echo(f"  {len(merged_seqs)} merged sequences → {len(conservation):,} passed conservation filter")
        if not conservation:
            click.echo(
                "\nNo merged sequences passed the conservation threshold.\n"
                "  • Lower --min-conservation (merged probes are longer and harder to conserve).\n"
                "  • Try --no-merge-overlapping to see unmerged candidates.",
                err=True,
            )
            sys.exit(1)

    # ── Pass 2: exact off-target scoring (streaming background) ──────────────
    click.echo("\nPass 2/2 — exact off-target scoring across background sequences …")

    candidate_set = set(conservation.keys())

    def _bg_stream():
        for _acc, taxonomy, seq in iter_silva(silva_path):
            if not taxon_matches(taxonomy.get(level, ""), target_name, level):
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

            max_tseqs = blast_max_target_seqs
            blast_out = run_blast(
                top_kmers,
                blast_db_path,
                threads=threads,
                perc_identity=blast_identity,
                max_target_seqs=max_tseqs,
            )
            blast_scores = parse_blast_results(
                blast_out, top_kmers, target_name, level, max_target_seqs=max_tseqs
            )

            for col, key in [
                ("blast_total_hits",           "blast_total_hits"),
                ("blast_target_hits",          "blast_target_hits"),
                ("blast_offtarget_hits",       "blast_offtarget_hits"),
                ("blast_weighted_offtarget",   "blast_weighted_offtarget"),
                ("blast_specificity",          "blast_specificity"),
                ("blast_weighted_specificity", "blast_weighted_specificity"),
                ("blast_top_offtarget",        "blast_top_offtarget"),
                ("blast_capped",               "blast_capped"),
            ]:
                default = "" if col == "blast_top_offtarget" else pd.NA
                df[col] = df["kmer"].map(lambda k, _k=key, _d=default: blast_scores.get(k, {}).get(_k, _d))

            n_capped = df["blast_capped"].sum() if "blast_capped" in df.columns else 0
            if n_capped:
                click.echo(
                    f"  WARNING: {n_capped} probe(s) hit the BLAST max_target_seqs={max_tseqs} cap — "
                    "blast_specificity may be underestimated for large target taxa.",
                    err=True,
                )

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

    display_cols = ["kmer", "gc_pct", "conservation_pct", "exact_offtarget", "exact_specificity"]
    if "blast_weighted_specificity" in df.columns:
        display_cols += ["blast_offtarget_hits", "blast_specificity", "blast_weighted_specificity", "blast_top_offtarget"]

    click.echo("\nTop 10 candidates:")
    click.echo(df[display_cols].head(10).to_string(index=False))


# ─── coverage ──────────────────────────────────────────────────────────────────

@cli.command()
@click.argument("silva_fasta")
@click.argument("target_name")
@click.option(
    "--level", "-l",
    type=click.Choice(TAXONOMY_LEVELS, case_sensitive=False),
    default="genus", show_default=True,
    help="Taxonomic level for TARGET_NAME.",
)
@click.option(
    "--breakdown", "-b",
    type=click.Choice(TAXONOMY_LEVELS, case_sensitive=False),
    default="species", show_default=True,
    help="Sub-level to break coverage down by.",
)
@click.option(
    "--kmer", "kmers", multiple=True, metavar="SEQ",
    help="K-mer sequence(s) to test (repeatable).",
)
@click.option(
    "--from-tsv", default=None,
    help="Auto-pick top k-mers from a 'bac-probes find' results TSV.",
)
@click.option(
    "--top-n", default=3, show_default=True,
    help="Number of top k-mers to pick from --from-tsv.",
)
@click.option(
    "--min-seqs", default=5, show_default=True,
    help="Minimum sequence count to show a breakdown row.",
)
@click.option(
    "--output", "-o", default=None,
    help="Write coverage table to this TSV path.",
)
def coverage(
    silva_fasta: str,
    target_name: str,
    level: str,
    breakdown: str,
    kmers: tuple[str, ...],
    from_tsv: str | None,
    top_n: int,
    min_seqs: int,
    output: str | None,
) -> None:
    """Show per-taxon k-mer coverage for probes within a target group.

    Scans SILVA sequences that belong to TARGET_NAME (at --level) and
    reports what fraction contain each k-mer, broken down by --breakdown
    sub-taxon (e.g., species within a genus).  Also reports pool coverage
    (any k-mer hit).

    \b
    SILVA_FASTA  Path to the SILVA NR99 FASTA (gzipped or plain).
    TARGET_NAME  Taxon to analyse (e.g. 'Fusobacterium').

    Supply k-mers with --kmer or auto-pick top probes with --from-tsv:

    \b
    Examples:
      bac-probes coverage silva_db/SILVA_*.fasta.gz Fusobacterium \\
          --level genus --breakdown species \\
          --kmer GATGGGGAAGCCAGCTTACTGGACAGATACTG \\
          --kmer ATGCAGGGCTCAACTCTGTATTGCGTTGGAAA

      bac-probes coverage silva_db/SILVA_*.fasta.gz Bacteroides \\
          --from-tsv bacteroides_genus_probes.tsv --top-n 2
    """
    # ── Resolve k-mers ─────────────────────────────────────────────────────────
    kmer_list: list[str] = list(kmers)
    if from_tsv:
        df_probes = pd.read_csv(from_tsv, sep="\t")
        sort_col = (
            "blast_weighted_specificity"
            if "blast_weighted_specificity" in df_probes.columns
            else "exact_specificity"
        )
        kmer_list = df_probes.nlargest(top_n, sort_col)["kmer"].tolist()
        click.echo(f"Using top {len(kmer_list)} k-mers from {from_tsv}:")
        for k in kmer_list:
            click.echo(f"  {k}")
        click.echo()

    if not kmer_list:
        click.echo("ERROR: supply at least one --kmer or use --from-tsv.", err=True)
        sys.exit(1)

    kmer_upper = [k.upper() for k in kmer_list]

    # ── Scan SILVA ─────────────────────────────────────────────────────────────
    total: dict[str, int] = {}
    hits: dict[str, dict[str, int]] = {k: {} for k in kmer_upper}
    pool_hits: dict[str, int] = {}

    click.echo(f"Scanning SILVA for '{target_name}' at level '{level}' …")
    for _acc, taxonomy, seq in tqdm(iter_silva(Path(silva_fasta)), desc="  reading", unit=" seq"):
        if not taxon_matches(taxonomy.get(level, ""), target_name, level):
            continue
        sub = taxonomy.get(breakdown, "") or "unknown"
        total[sub] = total.get(sub, 0) + 1
        any_hit = False
        for k in kmer_upper:
            if k in seq:
                hits[k][sub] = hits[k].get(sub, 0) + 1
                any_hit = True
        if any_hit:
            pool_hits[sub] = pool_hits.get(sub, 0) + 1

    if not total:
        click.echo(
            f"ERROR: no sequences found for '{target_name}' at level '{level}'.", err=True
        )
        sys.exit(1)

    # ── Build table ─────────────────────────────────────────────────────────────
    short_labels = [f"k{i+1}_{k[:8]}" for i, k in enumerate(kmer_upper)]

    rows = []
    for sub, n in sorted(total.items(), key=lambda x: -x[1]):
        if n < min_seqs:
            continue
        row: dict = {"subtaxon": sub, "total_seqs": n}
        for label, k in zip(short_labels, kmer_upper):
            h = hits[k].get(sub, 0)
            row[label] = f"{h}/{n} ({h/n*100:.0f}%)"
        ph = pool_hits.get(sub, 0)
        row["pool_coverage"] = f"{ph}/{n} ({ph/n*100:.0f}%)"
        rows.append(row)

    # Summary row (all sequences)
    n_all = sum(total.values())
    summary: dict = {"subtaxon": "ALL", "total_seqs": n_all}
    for label, k in zip(short_labels, kmer_upper):
        h_all = sum(hits[k].values())
        summary[label] = f"{h_all}/{n_all} ({h_all/n_all*100:.1f}%)"
    ph_all = sum(pool_hits.values())
    summary["pool_coverage"] = f"{ph_all}/{n_all} ({ph_all/n_all*100:.1f}%)"

    df_out = pd.concat([pd.DataFrame([summary]), pd.DataFrame(rows)], ignore_index=True)

    # ── Output ─────────────────────────────────────────────────────────────────
    click.echo(f"\nK-mer key:")
    for label, k in zip(short_labels, kmer_upper):
        click.echo(f"  {label} = {k}")
    click.echo()
    click.echo(df_out.to_string(index=False))

    if output:
        df_out.to_csv(output, sep="\t", index=False)
        click.echo(f"\nCoverage table written to {output}")


# ─── cocktail ──────────────────────────────────────────────────────────────────

@cli.command()
@click.argument("silva_fasta")
@click.argument("target_name")
@click.option(
    "--level", "-l",
    type=click.Choice(TAXONOMY_LEVELS, case_sensitive=False),
    default="genus", show_default=True,
    help="Taxonomic level for TARGET_NAME.",
)
@click.option(
    "--n-probes", default=8, show_default=True,
    help="Maximum number of probes in the cocktail.",
)
@click.option(
    "--coverage-target", default=0.90, show_default=True,
    help="Stop adding probes once this fraction of target sequences is covered.",
)
@click.option(
    "--min-conservation", default=0.20, show_default=True,
    help="Minimum fraction of target sequences a k-mer must appear in to be a candidate.",
)
@click.option(
    "--min-specificity", default=0.99, show_default=True,
    help="Minimum exact_specificity a candidate probe must have to be considered.",
)
@click.option("--gc-min", default=0.35, show_default=True, help="Minimum GC fraction.")
@click.option("--gc-max", default=0.65, show_default=True, help="Maximum GC fraction.")
@click.option("--max-homopolymer", default=5, show_default=True, help="Reject k-mers with homopolymer runs >= this length.")
@click.option("-k", "--kmer-size", default=32, show_default=True, help="K-mer length (bp).")
@click.option("--blast/--no-blast", "use_blast", default=True, show_default=True, help="Run BLAST on probes.")
@click.option(
    "--blast-select/--no-blast-select", default=True, show_default=True,
    help="BLAST top candidates before selection and use blast_weighted_specificity for greedy ranking. "
         "Requires --blast. When off, greedy uses exact_specificity and BLAST only scores the final selection.",
)
@click.option(
    "--top-n-blast", default=500, show_default=True,
    help="Number of top candidates (by exact specificity) to pre-BLAST when --blast-select is on.",
)
@click.option("--blast-db", default=None, help="BLAST database prefix (created by build-db).")
@click.option("--blast-identity", default=85.0, show_default=True, help="Minimum BLAST percent identity.")
@click.option("--blast-max-target-seqs", default=10000, show_default=True, help="BLAST max_target_seqs.")
@click.option("--threads", default=4, show_default=True, help="CPU threads for BLAST.")
@click.option("--output", "-o", default="cocktail.tsv", show_default=True, help="Output TSV path.")
def cocktail(
    silva_fasta: str,
    target_name: str,
    level: str,
    n_probes: int,
    coverage_target: float,
    min_conservation: float,
    min_specificity: float,
    gc_min: float,
    gc_max: float,
    max_homopolymer: int,
    kmer_size: int,
    use_blast: bool,
    blast_select: bool,
    top_n_blast: int,
    blast_db: str | None,
    blast_identity: float,
    blast_max_target_seqs: int,
    threads: int,
    output: str,
) -> None:
    """Design a probe cocktail covering TARGET_NAME using greedy set cover.

    Finds up to N_PROBES k-mers that together cover >= COVERAGE_TARGET of the
    target taxon's sequences. With --blast-select (default), BLASTs the top
    candidates first and ranks by blast_weighted_specificity so near-miss
    off-targets are penalised during probe selection.
    """
    silva_path = Path(silva_fasta)
    blast_db_path = Path(blast_db) if blast_db else silva_path.parent / "silva"

    # ── Pass 1: collect target sequences ────────────────────────────────────
    click.echo(f"\nPass 1/2 — scanning SILVA for '{target_name}' at level '{level}' …")
    target_seqs: list[str] = []
    n_bg_pass1 = 0
    for _acc, taxonomy, seq in tqdm(iter_silva(silva_path), desc="  reading", unit=" seq"):
        if taxon_matches(taxonomy.get(level, ""), target_name, level):
            target_seqs.append(seq)
        else:
            n_bg_pass1 += 1

    if not target_seqs:
        click.echo(
            f"\nERROR: No sequences found for '{target_name}' at level '{level}'.\n"
            "  • Check spelling (run 'bac-probes list-taxa' to browse names).",
            err=True,
        )
        sys.exit(1)

    click.echo(f"  Target sequences : {len(target_seqs):,}\n  Background seqs  : {n_bg_pass1:,}")

    # ── Conservation scoring with per-sequence coverage ──────────────────────
    click.echo(f"\nScoring {kmer_size}-mers (conservation ≥ {min_conservation:.0%}) …")
    coverage_sets, n_target = score_conservation_with_coverage(
        target_seqs, k=kmer_size, min_conservation=min_conservation,
        gc_range=(gc_min, gc_max), max_homopolymer=max_homopolymer,
    )

    if not coverage_sets:
        click.echo("\nNo k-mers passed the conservation threshold.", err=True)
        sys.exit(1)

    click.echo(f"  Candidates: {len(coverage_sets):,}")

    # ── Pass 2: exact off-target scoring (streaming) ─────────────────────────
    click.echo("\nPass 2/2 — exact off-target scoring …")

    def _bg_stream():
        for _acc, taxonomy, seq in iter_silva(silva_path):
            if not taxon_matches(taxonomy.get(level, ""), target_name, level):
                yield seq

    offtarget_counts, n_bg = score_offtarget_exact(
        set(coverage_sets.keys()),
        tqdm(_bg_stream(), desc="  scanning", total=n_bg_pass1, unit=" seq"),
        k=kmer_size,
    )

    exact_spec = {
        kmer: round(1.0 - offtarget_counts.get(kmer, 0) / n_bg, 6)
        for kmer in coverage_sets
    }

    n_pass = sum(1 for s in exact_spec.values() if s >= min_specificity)
    click.echo(f"  Candidates passing exact specificity ≥ {min_specificity}: {n_pass:,}")

    # ── Optional: BLAST top candidates before selection ──────────────────────
    blast_scores_all: dict = {}
    selection_spec = exact_spec  # default: use exact specificity for greedy
    db_exists = False

    if use_blast:
        db_exists = (
            blast_db_path.with_suffix(".nhr").exists()
            or (blast_db_path.parent / (blast_db_path.name + ".00.nhr")).exists()
        )
        if not db_exists:
            click.echo(
                f"\nWARNING: BLAST database not found at '{blast_db_path}'. Skipping BLAST.",
                err=True,
            )

    if use_blast and blast_select and db_exists:
        # Sort by exact specificity, take top N to BLAST
        top_candidates = sorted(
            [k for k in coverage_sets if exact_spec.get(k, 0) >= min_specificity],
            key=lambda k: (exact_spec[k], len(coverage_sets[k])),
            reverse=True,
        )[:top_n_blast]
        click.echo(f"\nPre-BLASTing top {len(top_candidates)} candidates for weighted specificity …")
        blast_out = run_blast(
            top_candidates, blast_db_path, threads=threads,
            perc_identity=blast_identity, max_target_seqs=blast_max_target_seqs,
        )
        blast_scores_all = parse_blast_results(
            blast_out, top_candidates, target_name, level,
            max_target_seqs=blast_max_target_seqs,
        )
        # Use blast_weighted_specificity for greedy; fall back to exact if missing
        selection_spec = {
            k: blast_scores_all[k].get("blast_weighted_specificity", exact_spec[k])
            if k in blast_scores_all else exact_spec[k]
            for k in top_candidates
        }
        # Restrict coverage_sets to only BLASTed candidates
        coverage_sets_sel = {k: coverage_sets[k] for k in top_candidates}
        n_blast_pass = sum(1 for s in selection_spec.values() if s >= min_specificity)
        click.echo(f"  Candidates passing blast_weighted_specificity ≥ {min_specificity}: {n_blast_pass:,}")
    else:
        coverage_sets_sel = coverage_sets

    # ── Greedy set cover ─────────────────────────────────────────────────────
    click.echo(f"\nRunning greedy set cover (target ≥ {coverage_target:.0%}, max {n_probes} probes) …")
    selected = greedy_cocktail(
        coverage_sets_sel, selection_spec, n_target,
        n_probes=n_probes, coverage_target=coverage_target,
        min_specificity=min_specificity,
    )

    if not selected:
        click.echo(
            "\nNo probes selected — no candidates met the specificity threshold.\n"
            "  • Lower --min-specificity.\n"
            "  • Lower --min-conservation to widen the candidate pool.\n"
            "  • Increase --top-n-blast to BLAST more candidates.",
            err=True,
        )
        sys.exit(1)

    # ── BLAST selected probes (if not already BLASTed during selection) ──────
    if use_blast and db_exists and not blast_select:
        click.echo(f"\nBLASTing {len(selected)} selected probes …")
        blast_out = run_blast(
            selected, blast_db_path, threads=threads,
            perc_identity=blast_identity, max_target_seqs=blast_max_target_seqs,
        )
        blast_scores_all = parse_blast_results(
            blast_out, selected, target_name, level,
            max_target_seqs=blast_max_target_seqs,
        )

    # ── Build results table ──────────────────────────────────────────────────
    cumulative: set[int] = set()
    rows = []
    for rank, kmer in enumerate(selected, 1):
        newly = coverage_sets[kmer] - cumulative
        cumulative |= coverage_sets[kmer]
        rows.append({
            "probe_rank":              rank,
            "kmer":                    kmer,
            "probe":                   reverse_complement(kmer),
            "gc_pct":                  round(gc_content(kmer) * 100, 1),
            "conservation_pct":        round(100 * len(coverage_sets[kmer]) / n_target, 2),
            "seqs_covered":            len(coverage_sets[kmer]),
            "new_seqs_covered":        len(newly),
            "cumulative_seqs":         len(cumulative),
            "cumulative_coverage_pct": round(100 * len(cumulative) / n_target, 2),
            "exact_offtarget":         offtarget_counts.get(kmer, 0),
            "exact_bg_seqs":           n_bg,
            "exact_specificity":       exact_spec[kmer],
        })

    df = pd.DataFrame(rows)
    final_coverage = round(100 * len(cumulative) / n_target, 1)
    click.echo(f"  Selected {len(selected)} probe(s) covering {final_coverage}% of {n_target:,} target sequences")

    if blast_scores_all:
        for col, key in [
            ("blast_total_hits",           "blast_total_hits"),
            ("blast_target_hits",          "blast_target_hits"),
            ("blast_offtarget_hits",       "blast_offtarget_hits"),
            ("blast_weighted_offtarget",   "offtarget_weighted"),
            ("blast_specificity",          "blast_specificity"),
            ("blast_weighted_specificity", "blast_weighted_specificity"),
            ("blast_capped",               "blast_capped"),
            ("blast_top_offtarget",        "blast_top_offtarget"),
        ]:
            df[col] = df["kmer"].map(lambda k, key=key: blast_scores_all.get(k, {}).get(key))

    # ── Output ───────────────────────────────────────────────────────────────
    click.echo(f"\n{'─'*70}")
    click.echo(f"Cocktail for {target_name}  ({len(selected)} probes, {final_coverage}% coverage)")
    click.echo(f"{'─'*70}")
    display_cols = ["probe_rank", "probe", "gc_pct", "conservation_pct",
                    "new_seqs_covered", "cumulative_coverage_pct", "exact_specificity"]
    if "blast_weighted_specificity" in df.columns:
        display_cols.append("blast_weighted_specificity")
    click.echo(df[display_cols].to_string(index=False))

    df.to_csv(output, sep="\t", index=False)
    click.echo(f"\nCocktail written to {output}")


# ─── validate ──────────────────────────────────────────────────────────────────

@cli.command()
@click.argument("silva_fasta")
@click.argument("input_csv")
@click.option(
    "--top-offtargets", default=5, show_default=True,
    help="Number of top off-target organism names to report per probe.",
)
@click.option(
    "--no-header", "no_header", is_flag=True, default=False,
    help="CSV has no header row; columns are assumed to be taxa, level, sequences.",
)
@click.option(
    "--blast/--no-blast", "use_blast", default=False, show_default=True,
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
    "--blast-max-target-seqs", default=10000, show_default=True,
    help="BLAST max_target_seqs cap.",
)
@click.option(
    "--threads", default=4, show_default=True,
    help="CPU threads for BLAST.",
)
@click.option(
    "--exclude-unnamed-congeners", is_flag=True, default=False,
    help="Species level only: also report specificity excluding off-target hits in the "
         "target's SILVA genus that have no species name ('uncultured bacterium', "
         "'Genus sp.'), which may be the target itself. Adds *_excl_unnamed columns.",
)
@click.option(
    "--output", "-o", default=None,
    help="Write results to this TSV path (default: print to stdout only).",
)
def validate(
    silva_fasta: str,
    input_csv: str,
    top_offtargets: int,
    no_header: bool,
    use_blast: bool,
    blast_db: str | None,
    blast_identity: float,
    blast_max_target_seqs: int,
    threads: int,
    exclude_unnamed_congeners: bool,
    output: str | None,
) -> None:
    """Validate probe sequences against SILVA: report coverage and off-target specificity.

    \b
    INPUT_CSV  CSV with columns: taxa, level, sequences
               One probe per row. 'sequences' is the probe sequence (either strand —
               both orientations are checked against SILVA).

    \b
    Output columns (exact scoring, always):
      taxa                – as supplied
      level               – as supplied
      sequence            – probe sequence as supplied
      probe_len           – length in bp
      target_seqs         – SILVA sequences belonging to this taxon
      coverage_count      – target seqs containing the probe (either strand)
      coverage_pct        – coverage_count / target_seqs × 100
      exact_offtarget     – background seqs with an exact hit
      exact_bg_seqs       – total background seqs scanned
      exact_specificity   – 1 − (exact_offtarget / exact_bg_seqs)
      top_offtargets      – most frequent off-target organism names (exact hits)

    \b
    Additional columns with --blast:
      blast_total_hits           – total BLAST hits
      blast_target_hits          – BLAST hits within the target taxon
      blast_offtarget_hits       – raw off-target BLAST hit count
      blast_weighted_offtarget   – off-target hits weighted by binding concern
      blast_specificity          – blast_target_hits / blast_total_hits
      blast_weighted_specificity – blast_target_hits / (blast_target_hits + blast_weighted_offtarget)
      blast_top_offtarget        – most frequent off-target taxon in BLAST results
      blast_capped               – True if BLAST hit the max_target_seqs limit

    \b
    Additional columns with --exclude-unnamed-congeners (species-level rows):
      target_genus                             – SILVA genus of the target (most common
                                                 genus among its sequences)
      exact_unnamed_congener_hits              – exact off-target hits in that genus
                                                 with no species name
      exact_offtarget_excl_unnamed             – exact_offtarget without them
      exact_specificity_excl_unnamed           – exact_specificity without them
      blast_unnamed_congener_hits              – same idea for BLAST hits (with --blast)
      blast_specificity_excl_unnamed           – blast_specificity without them
      blast_weighted_specificity_excl_unnamed  – blast_weighted_specificity without them
    """
    silva_path = Path(silva_fasta)

    # ── Read and validate input CSV ───────────────────────────────────────────
    df_in = (
        pd.read_csv(input_csv, header=None, names=["taxa", "level", "sequences"])
        if no_header
        else pd.read_csv(input_csv)
    )
    required = {"taxa", "level", "sequences"}
    missing = required - set(df_in.columns)
    if missing:
        click.echo(f"ERROR: CSV is missing column(s): {', '.join(sorted(missing))}", err=True)
        sys.exit(1)

    df_in["sequences"] = df_in["sequences"].str.strip().str.upper()
    bad_level = ~df_in["level"].isin(TAXONOMY_LEVELS)
    if bad_level.any():
        bad = df_in.loc[bad_level, "level"].unique().tolist()
        click.echo(
            f"ERROR: unknown level(s): {bad}\n"
            f"  Valid levels: {', '.join(TAXONOMY_LEVELS)}",
            err=True,
        )
        sys.exit(1)

    # ── Build per-group data structures for a single SILVA pass ──────────────
    # Each entry: probes, RCs, counters indexed parallel to df_in rows.
    # key = (taxa_lower, level) → group state dict
    groups: dict[tuple[str, str], dict] = {}
    group_order: list[tuple[str, str]] = []  # preserves input order

    for _, row in df_in.iterrows():
        key = (row["taxa"].lower(), row["level"])
        if key not in groups:
            groups[key] = {
                "taxa":              row["taxa"],
                "level":             row["level"],
                "probes":            [],
                "probes_rc":         [],
                "orig_seqs":         [],
                "target_hit":        [],
                "bg_hit":            [],
                "bg_offtarget_taxa": [],
                "target_total":      0,
                "bg_total":          0,
                # --exclude-unnamed-congeners bookkeeping (genus → count)
                "target_genera":     Counter(),
                "bg_unnamed_genera": Counter(),
                "bg_unnamed_hit":    [],
            }
            group_order.append(key)
        g = groups[key]
        seq_up = row["sequences"]
        g["probes"].append(seq_up)
        g["probes_rc"].append(reverse_complement(seq_up))
        g["orig_seqs"].append(seq_up)
        g["target_hit"].append(0)
        g["bg_hit"].append(0)
        g["bg_offtarget_taxa"].append(Counter())
        g["bg_unnamed_hit"].append(Counter())

    track_unnamed = exclude_unnamed_congeners and any(lv == "species" for _, lv in groups)
    if exclude_unnamed_congeners and not track_unnamed:
        click.echo(
            "WARNING: --exclude-unnamed-congeners only applies to species-level rows; ignoring.",
            err=True,
        )

    n_groups = len(groups)
    n_probes_total = len(df_in)
    click.echo(
        f"\nValidating {n_probes_total} probe(s) across {n_groups} taxon group(s) "
        f"in a single SILVA pass …"
    )

    # ── Single pass through SILVA ─────────────────────────────────────────────
    for _acc, taxonomy, seq in tqdm(iter_silva(silva_path), desc="  scanning", unit=" seq"):
        # Pre-compute the finest available taxon name for off-target labelling
        offtarget_name = ""
        for lv in reversed(TAXONOMY_LEVELS):
            name = taxonomy.get(lv, "")
            if name:
                offtarget_name = name
                break

        for key, g in groups.items():
            level = key[1]
            is_target = taxon_matches(taxonomy.get(level, ""), g["taxa"], level)
            unnamed_genus = (
                taxonomy.get("genus", "")
                if track_unnamed and level == "species" and not is_target
                and is_unnamed_species(taxonomy.get("species", ""))
                else None
            )
            if is_target:
                g["target_total"] += 1
                if track_unnamed:
                    g["target_genera"][taxonomy.get("genus", "")] += 1
                for i, (p, prc) in enumerate(zip(g["probes"], g["probes_rc"])):
                    if p in seq or prc in seq:
                        g["target_hit"][i] += 1
            else:
                g["bg_total"] += 1
                if unnamed_genus is not None:
                    g["bg_unnamed_genera"][unnamed_genus] += 1
                for i, (p, prc) in enumerate(zip(g["probes"], g["probes_rc"])):
                    if p in seq or prc in seq:
                        g["bg_hit"][i] += 1
                        if offtarget_name:
                            g["bg_offtarget_taxa"][i][offtarget_name] += 1
                        if unnamed_genus is not None:
                            g["bg_unnamed_hit"][i][unnamed_genus] += 1

    # ── Build output rows (in input order) ───────────────────────────────────
    rows_out: list[dict] = []
    for key in group_order:
        g = groups[key]
        taxa, level = g["taxa"], g["level"]
        bg_total = g["bg_total"]
        target_total = g["target_total"]

        if target_total == 0:
            click.echo(
                f"\nWARNING: no SILVA sequences found for '{taxa}' at level '{level}'.\n"
                "  Check spelling (run 'bac-probes list-taxa').",
                err=True,
            )

        for i, orig_seq in enumerate(g["orig_seqs"]):
            ot = g["bg_hit"][i]
            exact_spec = round(1.0 - ot / bg_total, 6) if bg_total > 0 else 1.0
            cov_pct = round(100 * g["target_hit"][i] / target_total, 2) if target_total > 0 else 0.0
            top_ot = ", ".join(
                f"{name} ({cnt})"
                for name, cnt in g["bg_offtarget_taxa"][i].most_common(top_offtargets)
            )
            row_out = {
                "taxa":              taxa,
                "level":             level,
                "sequence":          orig_seq,
                "probe_len":         len(orig_seq),
                "target_seqs":       target_total,
                "coverage_count":    g["target_hit"][i],
                "coverage_pct":      cov_pct,
                "exact_offtarget":   ot,
                "exact_bg_seqs":     bg_total,
                "exact_specificity": exact_spec,
                "top_offtargets":    top_ot,
            }
            if track_unnamed and level == "species":
                tg = g["target_genera"].most_common(1)[0][0] if g["target_genera"] else ""
                g["target_genus"] = tg
                unnamed_hits = g["bg_unnamed_hit"][i][tg] if tg else 0
                bg_ex = bg_total - (g["bg_unnamed_genera"][tg] if tg else 0)
                row_out.update({
                    "target_genus":                   tg,
                    "exact_unnamed_congener_hits":    unnamed_hits,
                    "exact_offtarget_excl_unnamed":   ot - unnamed_hits,
                    "exact_specificity_excl_unnamed":
                        round(1.0 - (ot - unnamed_hits) / bg_ex, 6) if bg_ex > 0 else 1.0,
                })
            rows_out.append(row_out)

    if not rows_out:
        click.echo("No results to report.", err=True)
        sys.exit(1)

    df_out = pd.DataFrame(rows_out)

    # ── Optional BLAST scoring (per taxon group) ──────────────────────────────
    if use_blast:
        blast_db_path = Path(blast_db) if blast_db else Path("silva_db/silva")
        db_exists = (
            blast_db_path.with_suffix(".nhr").exists()
            or (blast_db_path.parent / (blast_db_path.name + ".00.nhr")).exists()
        )
        if not db_exists:
            click.echo(
                f"\nWARNING: BLAST database not found at '{blast_db_path}'. Skipping BLAST.\n"
                "  Run 'bac-probes build-db' first or pass --blast-db.",
                err=True,
            )
        else:
            all_blast: dict[str, dict] = {}  # seq → blast score dict
            for key in group_order:
                g = groups[key]
                taxa, level = g["taxa"], g["level"]
                probes = g["orig_seqs"]
                click.echo(f"\nBLASTing {len(probes)} probe(s) for '{taxa}' …")
                blast_out = run_blast(
                    probes, blast_db_path, threads=threads,
                    perc_identity=blast_identity, max_target_seqs=blast_max_target_seqs,
                )
                blast_scores = parse_blast_results(
                    blast_out, probes, taxa, level,
                    max_target_seqs=blast_max_target_seqs,
                    target_genus=g.get("target_genus") if track_unnamed else None,
                )
                all_blast.update(blast_scores)

            for col, key in [
                ("blast_total_hits",           "blast_total_hits"),
                ("blast_target_hits",          "blast_target_hits"),
                ("blast_offtarget_hits",       "blast_offtarget_hits"),
                ("blast_weighted_offtarget",   "blast_weighted_offtarget"),
                ("blast_specificity",          "blast_specificity"),
                ("blast_weighted_specificity", "blast_weighted_specificity"),
                ("blast_top_offtarget",        "blast_top_offtarget"),
                ("blast_capped",               "blast_capped"),
            ] + ([
                ("blast_unnamed_congener_hits",             "blast_unnamed_congener_hits"),
                ("blast_specificity_excl_unnamed",          "blast_specificity_excl_unnamed"),
                ("blast_weighted_specificity_excl_unnamed", "blast_weighted_specificity_excl_unnamed"),
            ] if track_unnamed else []):
                default = "" if col == "blast_top_offtarget" else pd.NA
                df_out[col] = df_out["sequence"].map(
                    lambda s, k=key, d=default: all_blast.get(s, {}).get(k, d)
                )

    # ── Display ───────────────────────────────────────────────────────────────
    click.echo()
    display_cols = [
        "taxa", "level", "sequence", "probe_len",
        "target_seqs", "coverage_count", "coverage_pct",
        "exact_offtarget", "exact_specificity", "top_offtargets",
    ]
    if use_blast and "blast_weighted_specificity" in df_out.columns:
        display_cols += ["blast_weighted_specificity", "blast_top_offtarget"]
    click.echo(df_out[display_cols].to_string(index=False))

    if output:
        df_out.to_csv(output, sep="\t", index=False)
        click.echo(f"\nResults written to {output}")
