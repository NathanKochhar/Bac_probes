"""Command-line interface for bac-probes."""

import importlib.metadata
import shlex
import sys
from collections import Counter
from datetime import datetime
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
from .design import (
    DEFAULT_SPEC_FLOORS,
    Target,
    blast_pool,
    candidate_table,
    chunk_records,
    conservation_job,
    coverage_masks,
    init_offtarget_scan,
    offtarget_job,
    read_targets,
    run_jobs,
    safe_filename,
    select_probes,
    selected_table,
    selection_pool,
)
from .kmers import reverse_complement
from .specificity import parse_blast_results, run_blast


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option()
def cli() -> None:
    """bac-probes: design and validate taxon-specific 16S rRNA probes from SILVA.

    \b
    Typical workflow:
      1. bac-probes download --output-dir ./silva_db
      2. bac-probes build-db ./silva_db/SILVA_*.fasta.gz --db-dir ./silva_db
      3. bac-probes design   ./silva_db/SILVA_*.fasta.gz  Streptococcus  --level genus
      4. bac-probes validate ./silva_db/SILVA_*.fasta.gz  my_probes.csv  --blast
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


# ─── design ────────────────────────────────────────────────────────────────────

@cli.command()
@click.argument("silva_fasta")
@click.argument("target")
@click.option(
    "--level", "-l",
    type=click.Choice(TAXONOMY_LEVELS, case_sensitive=False),
    default="genus", show_default=True,
    help="Taxonomic level for a single TARGET, or for CSV rows without a level.",
)
@click.option(
    "--target-fasta", default=None,
    help="Single-taxon mode only: take target sequences from this FASTA (e.g. an NCBI "
         "download) instead of SILVA. SILVA is still the background for off-targets.",
)
@click.option("--kmer-size", "-k", default=32, show_default=True, help="Probe length (bp).")
@click.option(
    "--min-conservation", default=0.02, show_default=True,
    help="Keep candidates present in at least this fraction of target sequences.",
)
@click.option("--gc-min", default=0.35, show_default=True, help="Minimum probe GC fraction.")
@click.option("--gc-max", default=0.65, show_default=True, help="Maximum probe GC fraction.")
@click.option(
    "--max-homopolymer", default=5, show_default=True,
    help="Reject probes with a homopolymer run >= this length.",
)
@click.option(
    "--merge-overlapping/--no-merge-overlapping", default=False, show_default=True,
    help="Merge overlapping k-mers into longer probes before scoring.",
)
@click.option(
    "--exclude-unnamed-congeners", is_flag=True, default=False,
    help="Species-level targets: don't count off-target hits in the target's SILVA genus "
         "that have no species name ('uncultured bacterium', 'Genus sp.') when selecting. "
         "Adds *_excl_unnamed columns.",
)
@click.option(
    "--blast/--no-blast", "use_blast", default=False, show_default=True,
    help="Score near-match off-targets with BLAST and select on blast_weighted_specificity.",
)
@click.option("--blast-db", default=None, help="BLAST database prefix (default: silva_db/silva).")
@click.option("--blast-identity", default=85.0, show_default=True, help="Minimum BLAST percent identity.")
@click.option(
    "--blast-max-target-seqs", default=50000, show_default=True,
    help="BLAST max_target_seqs; must exceed the largest target's sequence count.",
)
@click.option(
    "--blast-pool", "blast_pool_size", default=120, show_default=True,
    help="Candidates BLASTed per taxon (distinct sites; half by coverage, half by specificity).",
)
@click.option(
    "--select/--no-select", "do_select", default=True, show_default=True,
    help="Select a probe set per taxon. --no-select only writes candidate tables.",
)
@click.option("--n-probes", default=10, show_default=True, help="Maximum probes selected per taxon.")
@click.option(
    "--pool-target", default=0.50, show_default=True,
    help="Pooled coverage to aim for: the strictest specificity floor reaching it is used.",
)
@click.option(
    "--spec-floors", default=",".join(f"{f:g}" for f in DEFAULT_SPEC_FLOORS), show_default=True,
    help="Specificity floors to try, comma-separated (strictest first is not required).",
)
@click.option(
    "--candidates/--no-candidates", "write_candidates", default=True, show_default=True,
    help="Write the full candidate table per taxon to candidates/.",
)
@click.option("--threads", default=4, show_default=True, help="Worker processes / BLAST threads.")
@click.option(
    "--output-dir", "-o", default="bac_probes_design", show_default=True,
    help="Directory for all outputs (created if missing).",
)
def design(
    silva_fasta: str,
    target: str,
    level: str,
    target_fasta: str | None,
    kmer_size: int,
    min_conservation: float,
    gc_min: float,
    gc_max: float,
    max_homopolymer: int,
    merge_overlapping: bool,
    exclude_unnamed_congeners: bool,
    use_blast: bool,
    blast_db: str | None,
    blast_identity: float,
    blast_max_target_seqs: int,
    blast_pool_size: int,
    do_select: bool,
    n_probes: int,
    pool_target: float,
    spec_floors: str,
    write_candidates: bool,
    threads: int,
    output_dir: str,
) -> None:
    """Design probes for one taxon or a CSV list of taxa.

    \b
    SILVA_FASTA  SILVA NR99 FASTA (gzipped or plain).
    TARGET       A taxon name (with --level), or a CSV with a 'taxa' column and
                 optional 'level' column, one taxon per row.

    \b
    For every taxon it finds candidate k-mers present in >= --min-conservation of
    the taxon's SILVA sequences, counts exact off-target matches across the rest
    of SILVA (one pass for all taxa), optionally BLASTs the best candidates, and
    selects up to --n-probes probes at distinct sites that together cover
    --pool-target of the taxon, using the strictest --spec-floors value that gets
    there. Specificity is exact_precision (target hits / all SILVA hits), or
    blast_weighted_specificity with --blast.

    \b
    Outputs in --output-dir:
      design_summary.csv           one row per taxon
      all_selected_probes.csv      every selected probe
      selected/<taxon>_selected_probes.csv
      candidates/<taxon>_candidates.tsv   full ranked candidate table
      command.txt                  the command, version and settings used

    \b
    Examples:
      bac-probes design silva.fasta.gz Fusobacterium --level genus
      bac-probes design silva.fasta.gz taxa.csv --blast --threads 14 -o results/
    """
    silva_path = Path(silva_fasta)
    out_dir = Path(output_dir)
    try:
        floors = [float(f) for f in spec_floors.split(",") if f.strip()]
        pairs = read_targets(target, level.lower())
    except ValueError as exc:
        click.echo(f"ERROR: {exc}", err=True)
        sys.exit(1)
    if not floors:
        click.echo("ERROR: --spec-floors is empty.", err=True)
        sys.exit(1)
    if target_fasta and len(pairs) > 1:
        click.echo("ERROR: --target-fasta works with a single taxon only.", err=True)
        sys.exit(1)

    targets = [
        Target(name, lv, exclude_unnamed=exclude_unnamed_congeners and lv == "species")
        for name, lv in pairs
    ]
    if exclude_unnamed_congeners and not any(t.exclude_unnamed for t in targets):
        click.echo("WARNING: --exclude-unnamed-congeners only applies to species-level taxa.", err=True)
    gc_range = (gc_min, gc_max)

    # ── Pass 1: target sequences ──────────────────────────────────────────────
    click.echo(f"\nDesigning probes for {len(targets)} taxon/taxa. Pass 1/2 — reading SILVA …")
    n_total = 0
    unnamed_by_genus: Counter = Counter()
    for _acc, taxonomy, seq in tqdm(iter_silva(silva_path), desc="  reading", unit=" seq"):
        n_total += 1
        if is_unnamed_species(taxonomy.get("species", "")):
            unnamed_by_genus[taxonomy.get("genus", "")] += 1
        for t in targets:
            if taxon_matches(taxonomy.get(t.level, ""), t.name, t.level):
                t.genera[taxonomy.get("genus", "")] += 1
                if not target_fasta:
                    t.seqs.append(seq)
    if target_fasta:
        t = targets[0]
        t.seqs = [str(r.seq).upper().replace("U", "T") for r in SeqIO.parse(target_fasta, "fasta")]
        click.echo(f"  {len(t.seqs):,} target sequences from {target_fasta}")
    for t in targets:
        if t.exclude_unnamed:
            t.n_ambiguous = unnamed_by_genus[t.genus]

    missing = [t for t in targets if not t.seqs]
    for t in missing:
        click.echo(
            f"WARNING: no sequences for '{t.name}' at level '{t.level}' — skipped. "
            "Check spelling with 'bac-probes list-taxa'.",
            err=True,
        )
    targets = [t for t in targets if t.seqs]
    if not targets:
        click.echo("ERROR: none of the taxa were found in SILVA.", err=True)
        sys.exit(1)

    # ── Candidates per target ─────────────────────────────────────────────────
    click.echo(f"\nScoring {kmer_size}-mers (conservation ≥ {min_conservation:.0%}) …")
    jobs = [(i, t.seqs, kmer_size, min_conservation, gc_range, max_homopolymer, merge_overlapping)
            for i, t in sorted(enumerate(targets), key=lambda it: -len(it[1].seqs))]
    hits: dict[int, dict[str, int]] = dict(run_jobs(conservation_job, jobs, threads))
    owners: dict[str, list[int]] = {}
    for i, h in hits.items():
        for probe_seq in h:
            owners.setdefault(probe_seq, []).append(i)
    for i in sorted(hits):
        t = targets[i]
        click.echo(f"  {t.name}: {len(t.seqs):,} seqs, {len(hits[i]):,} candidates")

    # ── Pass 2: exact off-targets for all candidates ──────────────────────────
    click.echo("\nPass 2/2 — exact off-target scan …")
    light = [Target(t.name, t.level, t.exclude_unnamed, [], t.genera, t.n_ambiguous) for t in targets]
    records = (
        (taxonomy, seq)
        for _acc, taxonomy, seq in tqdm(iter_silva(silva_path), desc="  scanning", total=n_total, unit=" seq")
    )
    off: Counter = Counter()
    unnamed: Counter = Counter()
    for o, u in run_jobs(offtarget_job, chunk_records(records), threads,
                         initializer=init_offtarget_scan, initargs=(owners, light, kmer_size)):
        off.update(o)
        unnamed.update(u)

    tables = {i: candidate_table(t, hits[i], off, unnamed, i, n_total) for i, t in enumerate(targets)}
    metric_exact = {i: ("exact_precision_excl_unnamed" if t.exclude_unnamed else "exact_precision")
                    for i, t in enumerate(targets)}
    min_floor = min(floors)

    # ── Optional BLAST ────────────────────────────────────────────────────────
    blast_cols = ["blast_total_hits", "blast_target_hits", "blast_offtarget_hits",
                  "blast_weighted_offtarget", "blast_specificity", "blast_weighted_specificity",
                  "blast_top_offtarget", "blast_capped"]
    if use_blast:
        blast_db_path = Path(blast_db) if blast_db else Path("silva_db/silva")
        if not (blast_db_path.with_suffix(".nhr").exists()
                or (blast_db_path.parent / (blast_db_path.name + ".00.nhr")).exists()):
            click.echo(f"ERROR: BLAST database not found at '{blast_db_path}'. Run 'bac-probes build-db' "
                       "or drop --blast.", err=True)
            sys.exit(1)
        click.echo(f"\nBLAST (identity ≥ {blast_identity:g}%, max_target_seqs {blast_max_target_seqs}) …")
        for i, t in enumerate(targets):
            df = tables[i]
            if df.empty:
                continue
            pool = blast_pool(df, metric_exact[i], min_floor, blast_pool_size)
            scores: dict[str, dict] = {}
            for j in range(0, len(pool), 20):
                chunk = pool[j:j + 20]
                raw = run_blast(chunk, blast_db_path, threads=threads,
                                perc_identity=blast_identity, max_target_seqs=blast_max_target_seqs)
                scores.update(parse_blast_results(
                    raw, chunk, t.name, t.level, max_target_seqs=blast_max_target_seqs,
                    target_genus=t.genus if t.exclude_unnamed else None,
                ))
            cols = blast_cols + (["blast_unnamed_congener_hits", "blast_specificity_excl_unnamed",
                                  "blast_weighted_specificity_excl_unnamed"] if t.exclude_unnamed else [])
            for col in cols:
                df[col] = df.kmer.map(lambda p, c=col: scores.get(p, {}).get(c, pd.NA))
            n_capped = sum(1 for s in scores.values() if s.get("blast_capped"))
            click.echo(
                f"  {t.name}: {len(pool)} BLASTed"
                + (f" ({n_capped} hit the max_target_seqs cap; their BLAST scores are "
                   "incomplete — see blast_capped)" if n_capped else "")
            )

    # ── Outputs ───────────────────────────────────────────────────────────────
    out_dir.mkdir(parents=True, exist_ok=True)
    if write_candidates or not do_select:
        (out_dir / "candidates").mkdir(exist_ok=True)
        for i, t in enumerate(targets):
            tables[i].to_csv(out_dir / "candidates" / f"{safe_filename(t.name)}_candidates.tsv",
                             sep="\t", index=False)

    summary_rows, all_selected = [], []
    if do_select:
        click.echo("\nSelecting probes …")
        (out_dir / "selected").mkdir(exist_ok=True)
        for i, t in enumerate(targets):
            df = tables[i]
            if use_blast:
                metric = ("blast_weighted_specificity_excl_unnamed" if t.exclude_unnamed
                          else "blast_weighted_specificity")
                df_sel = df[df[metric].notna()].copy() if metric in df else df.iloc[0:0]
                df_sel[metric] = df_sel[metric].astype(float)
            else:
                metric = metric_exact[i]
                df_sel = df
            pool = list(df_sel.kmer) if use_blast else selection_pool(df_sel, metric, min_floor)
            masks = coverage_masks(t.seqs, pool)
            floor, picked, pooled = select_probes(df_sel, masks, len(t.seqs), metric, floors,
                                                  n_probes, pool_target)
            sel = selected_table(df_sel, picked, len(t.seqs))
            sel.insert(0, "level", t.level)
            sel.insert(0, "taxa", t.name)
            sel.to_csv(out_dir / "selected" / f"{safe_filename(t.name)}_selected_probes.csv", index=False)
            all_selected.append(sel)
            row = {
                "taxa": t.name, "level": t.level, "target_seqs": len(t.seqs),
                "candidates_found": len(df), "candidates_considered": len(pool),
                "spec_metric": metric, "spec_floor_used": floor if picked else None,
                "probes_selected": len(sel),
                "pooled_coverage_pct": round(100 * pooled, 2),
                "goal_reached": pooled >= pool_target,
                "conservation_min_pct": sel.conservation_pct.min() if len(sel) else None,
                "conservation_max_pct": sel.conservation_pct.max() if len(sel) else None,
                "spec_min": sel[metric].min() if len(sel) else None,
                "spec_mean": round(sel[metric].mean(), 4) if len(sel) else None,
                "exact_offtarget_total": int(sel.exact_offtarget.sum()) if len(sel) else None,
            }
            if t.exclude_unnamed:
                row["silva_genus"] = t.genus
                row["ambiguous_seqs_excluded"] = t.n_ambiguous
            if use_blast and len(sel):
                row["probes_blast_capped"] = int(sel.blast_capped.fillna(False).astype(bool).sum())
                top = Counter(x for x in sel.blast_top_offtarget if isinstance(x, str) and x)
                row["top_offtargets"] = "; ".join(f"{n} ({c})" for n, c in top.most_common(3))
            summary_rows.append(row)
            click.echo(f"  {t.name}: floor {floor} → {len(sel)} probes, pooled {100 * pooled:.1f}%")

        pd.concat(all_selected, ignore_index=True).to_csv(out_dir / "all_selected_probes.csv", index=False)
        summary = pd.DataFrame(summary_rows)
        summary.to_csv(out_dir / "design_summary.csv", index=False)
        click.echo("\n" + summary[["taxa", "target_seqs", "spec_floor_used", "probes_selected",
                                   "pooled_coverage_pct", "goal_reached"]].to_string(index=False))

    _write_command_file(out_dir, click.get_current_context())
    click.echo(f"\nResults written to {out_dir}/")


def _write_command_file(out_dir: Path, ctx: click.Context) -> None:
    """Record how the outputs were produced: command line, version, date and all settings."""
    try:
        version = importlib.metadata.version("bac-probes")
    except importlib.metadata.PackageNotFoundError:
        version = "unknown"
    argv = sys.argv[1:] if Path(sys.argv[0]).stem == "bac-probes" else None
    lines = [
        f"# bac-probes {version}, run {datetime.now().isoformat(timespec='seconds')}",
        "bac-probes " + (shlex.join(argv) if argv else f"{ctx.info_name} ..."),
        "",
        "# resolved settings",
    ] + [f"{k} = {v}" for k, v in ctx.params.items()]
    (out_dir / "command.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


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
    "--breakdown", "-b",
    type=click.Choice(TAXONOMY_LEVELS, case_sensitive=False), default=None,
    help="Also report coverage of each taxon's probes broken down by this sub-level "
         "(e.g. species within a genus), with pooled coverage. Written to "
         "<output>_breakdown.tsv when --output is given.",
)
@click.option(
    "--min-seqs", default=5, show_default=True,
    help="With --breakdown: hide sub-taxa with fewer sequences than this.",
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
    breakdown: str | None,
    min_seqs: int,
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

    \b
    With --breakdown LEVEL, a second table shows, for each taxon, how many of its
    sequences in each sub-taxon (e.g. each species of a genus) each probe hits,
    plus pooled coverage (hit by any of the taxon's probes). Probes are numbered
    k1, k2, … in input order within each taxon; the first row per taxon is ALL.
    """
    silva_path = Path(silva_fasta)
    if breakdown:
        breakdown = breakdown.lower()

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
                # --breakdown bookkeeping (sub-taxon → count)
                "bd_total":          Counter(),
                "bd_pool":           Counter(),
                "bd_hit":            [],
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
        g["bd_hit"].append(Counter())

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
                sub = (taxonomy.get(breakdown, "") or "unknown") if breakdown else None
                any_hit = False
                for i, (p, prc) in enumerate(zip(g["probes"], g["probes_rc"])):
                    if p in seq or prc in seq:
                        g["target_hit"][i] += 1
                        any_hit = True
                        if sub is not None:
                            g["bd_hit"][i][sub] += 1
                if sub is not None:
                    g["bd_total"][sub] += 1
                    if any_hit:
                        g["bd_pool"][sub] += 1
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

    if breakdown:
        df_bd = _breakdown_table(groups, group_order, min_seqs)
        click.echo(f"\nCoverage by {breakdown} (k1, k2, … = probes in input order per taxon):")
        for key in group_order:
            g = groups[key]
            for i, p in enumerate(g["orig_seqs"], 1):
                click.echo(f"  {g['taxa']}  k{i} = {p}")
        click.echo()
        click.echo(df_bd.to_string(index=False))
        if output:
            out_path = Path(output)
            bd_path = out_path.with_name(f"{out_path.stem}_breakdown{out_path.suffix or '.tsv'}")
            df_bd.to_csv(bd_path, sep="\t", index=False)
            click.echo(f"\nBreakdown written to {bd_path}")


def _fmt_cov(hits: int, total: int, decimals: int = 0) -> str:
    return f"{hits}/{total} ({hits / total * 100:.{decimals}f}%)" if total else "0/0"


def _breakdown_table(groups: dict, group_order: list, min_seqs: int) -> pd.DataFrame:
    """Per-taxon coverage table broken down by sub-taxon (ALL row first)."""
    rows = []
    for key in group_order:
        g = groups[key]
        n_probes = len(g["orig_seqs"])
        n_all = sum(g["bd_total"].values())
        all_row = {"taxa": g["taxa"], "level": g["level"], "subtaxon": "ALL", "total_seqs": n_all}
        for i in range(n_probes):
            all_row[f"k{i + 1}"] = _fmt_cov(sum(g["bd_hit"][i].values()), n_all, 1)
        all_row["pool_coverage"] = _fmt_cov(sum(g["bd_pool"].values()), n_all, 1)
        rows.append(all_row)
        for sub, n in sorted(g["bd_total"].items(), key=lambda x: (-x[1], x[0])):
            if n < min_seqs:
                continue
            row = {"taxa": g["taxa"], "level": g["level"], "subtaxon": sub, "total_seqs": n}
            for i in range(n_probes):
                row[f"k{i + 1}"] = _fmt_cov(g["bd_hit"][i][sub], n)
            row["pool_coverage"] = _fmt_cov(g["bd_pool"][sub], n)
            rows.append(row)
    df = pd.DataFrame(rows)
    k_cols = sorted((c for c in df.columns if c.startswith("k") and c[1:].isdigit()), key=lambda c: int(c[1:]))
    return df[["taxa", "level", "subtaxon", "total_seqs"] + k_cols + ["pool_coverage"]].fillna("")
