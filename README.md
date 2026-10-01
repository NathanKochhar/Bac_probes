# bac-probes

Design and validate bacterial 16S rRNA probes that are **conserved within a taxonomic group** (phylum to species) and **specific** — i.e., unlikely to bind outside of it.

Designed for **Xenium in-situ probe design**, with position-weighted off-target scoring that down-penalises near-matches whose mismatches fall near the probe centre (where mismatches are most destabilising).

Uses the [SILVA 138.1 NR99](https://www.arb-silva.de/) SSU reference database (~510,000 curated 16S sequences) as the sequence source.

## How it works

1. **Conservation scan** — all k-mers (default 32-mers) from the target taxon's 16S sequences are counted. A k-mer's conservation score is the fraction of target sequences that contain it exactly.
2. **Exact off-target scan** — candidate k-mers are checked against every non-target SILVA sequence for exact matches. This is fast and memory-efficient because only the candidates (not all background k-mers) are held in memory.
3. **BLAST near-match scoring** *(optional)* — top candidates are BLASTed against the full SILVA database using `blastn-short`. Each hit is classified as on-target or off-target using SILVA taxonomy. Off-target hits are scored by **effective binding probability**: mismatches near the probe centre (which destabilise hybridisation most) contribute less to the off-target score than terminal mismatches.
4. **Probe selection** — `design` picks up to `--n-probes` probes at distinct sites that together cover `--pool-target` of the taxon's sequences (pooled coverage), using the strictest specificity floor that gets there. See [Probe selection](#probe-selection).

### Position-weighted scoring (Xenium-specific)

Each off-target BLAST hit is assigned a **specificity score** in [0, 1], where **1 = specific** (off-target won't bind) and **0 = non-specific** (off-target will bind):

```
weight(pos) = sin(π × (pos + 0.5) / k)²
specificity_score = min(1, Σ weight(mismatch_positions))
```

- A perfect-match hit scores **0.0** (certain off-target binding — probe is not protected).
- A single central mismatch (pos ≈ k/2) has weight ≈ 1.0 → specificity score ≈ **1.0** (probe won't hybridise — well protected).
- A single terminal mismatch (pos = 0 or k−1) has weight ≈ 0.002 → specificity score ≈ **0.002** (probe may still bind — minimal protection).

The `blast_weighted_specificity` column uses these scores instead of raw hit counts, making it a better predictor of real-world probe performance than `blast_specificity`.

### Candidate columns (`design`)

| Column | Description |
|---|---|
| `kmer` | K-mer sequence (same strand as SILVA reference) |
| `probe` | Reverse complement — order this as the in-situ probe sequence |
| `length` | Probe length (bp); longer than `-k` only with `--merge-overlapping` |
| `gc_pct` | GC content (%) |
| `conservation_pct` | % of target sequences containing this k-mer exactly |
| `target_seqs` | Number of target sequences used |
| `target_hits` | Target sequences containing the k-mer |
| `exact_offtarget` | Non-target SILVA sequences containing the k-mer exactly |
| `exact_bg_seqs` | Total non-target sequences scanned |
| `exact_specificity` | `1 − (exact_offtarget / exact_bg_seqs)` |
| `exact_precision` | `target_hits / (target_hits + exact_offtarget)` — share of all SILVA sequences carrying the k-mer that belong to the taxon; the selection metric without `--blast` |
| `blast_total_hits` | Total BLAST hits (target + off-target) |
| `blast_target_hits` | BLAST hits within the target taxon |
| `blast_offtarget_hits` | Raw count of off-target BLAST hits |
| `blast_weighted_offtarget` | Sum of effective binding scores for off-target hits; down-weights near-misses with central mismatches |
| `blast_specificity` | `blast_target_hits / blast_total_hits` (raw) |
| `blast_weighted_specificity` | `blast_target_hits / (blast_target_hits + blast_weighted_offtarget)` — **primary ranking metric** |
| `blast_capped` | `True` when BLAST hit the `max_target_seqs` limit; specificity estimates may be unreliable for large target taxa |
| `blast_top_offtarget` | Most frequent off-target taxon in BLAST results |

`blast_*` columns are filled for the candidates that were BLASTed (`--blast-pool` per taxon). Candidate tables are sorted by `exact_precision`, then `conservation_pct`, then k-mer.

## Installation

```bash
pip install -e .
```

BLAST+ is required for `--blast` scoring:

```bash
# Conda (Linux/macOS):
conda install -c bioconda blast

# Homebrew (macOS):
brew install blast

# Or download directly from NCBI:
# https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
```

`bac-probes` auto-detects BLAST+ in common locations (`/opt/homebrew`, `/usr/local`, `/tmp/ncbi-blast-*/`) so `export PATH=...` is not required.

## Workflow

### 1. Download SILVA

```bash
bac-probes download --output-dir ./silva_db
```

Downloads `SILVA_138.1_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz` (~1 GB compressed) into `./silva_db/`. Gap characters (`.` / `-`) from the multiple sequence alignment are stripped automatically on read.

### 2. Build BLAST database *(optional but recommended)*

```bash
bac-probes build-db \
    ./silva_db/SILVA_138.1_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz \
    --db-dir ./silva_db
```

Strips alignment gaps and runs `makeblastdb`. Takes a few minutes. Required for `--blast` scoring.

### 3. Browse taxon names

```bash
bac-probes list-taxa \
    ./silva_db/SILVA_138.1_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz \
    --level genus --top 50
```

Useful for finding the exact spelling `design` and `validate` need.

### 4. Design probes

`bac-probes design` takes one taxon, or a CSV listing many, and writes candidate tables plus a selected probe set per taxon. All taxa share one pass over SILVA, and the work is spread over `--threads` processes.

```bash
# One genus, exact scoring only
bac-probes design \
    ./silva_db/SILVA_138.1_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz \
    Bacteroides --level genus \
    --output-dir bacteroides/

# A list of taxa, with BLAST near-match scoring, 14 worker processes
bac-probes design \
    ./silva_db/SILVA_138.1_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz \
    my_families.csv --blast --threads 14 \
    --output-dir families/

# Species, treating unnamed same-genus sequences as possibly the target
bac-probes design \
    ./silva_db/SILVA_138.1_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz \
    "Fusobacterium nucleatum" --level species \
    --blast --exclude-unnamed-congeners \
    --output-dir f_nucleatum/

# Candidates only (no selection), target sequences from an NCBI download
bac-probes design \
    ./silva_db/SILVA_138.1_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz \
    "Fusobacterium nucleatum" --level species \
    --target-fasta my_ncbi_sequences.fasta --no-select \
    --output-dir f_nucleatum_candidates/
```

A taxa CSV needs a `taxa` column and may have a `level` column (any case); rows without a level use `--level`:

```csv
Taxa,Level
Staphylococcaceae,Family
Escherichia coli,Species
```

Outputs in `--output-dir`:

| File | Contents |
|---|---|
| `design_summary.csv` | One row per taxon: target sequences, floor used, probes selected, pooled coverage, whether `--pool-target` was reached, specificity range, top off-targets (with `--blast`) |
| `all_selected_probes.csv` | Every selected probe, all taxa |
| `selected/<taxon>_selected_probes.csv` | Selected probes in pick order, with `probe_rank`, `new_seqs_covered` and `cumulative_coverage_pct` |
| `candidates/<taxon>_candidates.tsv` | Full ranked candidate table (skip with `--no-candidates`) |
| `command.txt` | The command line, bac-probes version, date and every resolved setting |

#### Probe selection

For each taxon, `design` tries the `--spec-floors` from strictest to loosest (default `1, 0.999, 0.995, 0.99, 0.98, 0.95, 0.9`). At each floor it runs a greedy set cover over the candidates at or above it: every pick is the probe that adds the most not-yet-covered target sequences, skipping probes that share a 16-mer with one already picked (same site). It stops at the first floor whose set reaches `--pool-target`; if none does, it keeps the floor with the highest pooled coverage. If coverage stops growing before `--n-probes`, remaining slots are filled with further specific probes at distinct sites, which adds signal without losing specificity.

The specificity used is `exact_precision`, or `blast_weighted_specificity` with `--blast` (the `_excl_unnamed` versions with `--exclude-unnamed-congeners`). With `--blast`, only the `--blast-pool` candidates per taxon are BLASTed and considered: half chosen for coverage and half for exact specificity, at distinct sites. Ties are always broken on the k-mer sequence, so repeated runs give identical results.

### 5. Validate existing probes

`bac-probes validate` scores probe sequences you already have — from a previous run, a paper, or a vendor — against SILVA. It reports coverage (how many target sequences contain the probe) and exact off-target specificity in a single SILVA pass, regardless of how many taxa are in the input file.

```bash
bac-probes validate \
    ./silva_db/SILVA_138.1_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz \
    my_probes.csv \
    --output validation_results.tsv
```

The input CSV requires three columns: `taxa`, `level`, `sequences` (one probe per row). If your file has no header row, add `--no-header`.

```csv
taxa,level,sequences
Fusobacterium nucleatum,species,TCGCGATTACTAGCGATTCCAACTTCATGTAC
Bacteroides,genus,ATCGATCGATCGATCGATCGATCGATCGATCG
Firmicutes,phylum,GCGTATCCGGATTTATGGGCGTAAAGCGCGTC
```

Both orientations of each probe are checked (forward and reverse complement), so it does not matter whether you supply the k-mer or its reverse complement.

At species level, SILVA labels many sequences with a strain or subspecies (e.g. `Helicobacter pylori 26695`, `Campylobacter jejuni subsp. doylei`). A species target therefore also matches any SILVA name that starts with it followed by a space, so strains count as target rather than off-target. This applies to `design`, `validate` and BLAST hit classification. Other levels still require an exact name match.

Add `--blast` to also score near-match off-target binding via BLAST (same position-weighted scoring as `bac-probes design`). When using `--blast` with `validate`, the probe sequence you supply is treated as the full k-mer — its length determines the position-weighting window. Mismatches in BLAST hits are scored by their distance from the centre of your probe sequence: a mismatch at the centre contributes close to 1.0 (high protection against off-target binding), while terminal mismatches contribute close to 0.0 (minimal protection). This means supplying a truncated or padded sequence will shift the weighting incorrectly, so pass the exact probe length you intend to use.

```bash
bac-probes validate \
    ./silva_db/SILVA_138.1_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz \
    my_probes.csv \
    --blast \
    --blast-max-target-seqs 10000 \
    --output validation_results.tsv
```

Output columns (exact scoring, always present):

| Column | Description |
|---|---|
| `taxa` | As supplied |
| `level` | As supplied |
| `sequence` | Probe sequence as supplied |
| `probe_len` | Length in bp |
| `target_seqs` | SILVA sequences belonging to this taxon |
| `coverage_count` | Target seqs containing the probe (either strand) |
| `coverage_pct` | `coverage_count / target_seqs × 100` |
| `exact_offtarget` | Background seqs with an exact hit |
| `exact_bg_seqs` | Total background seqs scanned |
| `exact_specificity` | `1 − (exact_offtarget / exact_bg_seqs)` |
| `top_offtargets` | Most frequent off-target organism names (exact hits) |

Additional columns added by `--blast` (same definitions as `bac-probes design`): `blast_total_hits`, `blast_target_hits`, `blast_offtarget_hits`, `blast_weighted_offtarget`, `blast_specificity`, `blast_weighted_specificity`, `blast_top_offtarget`, `blast_capped`.

#### Unnamed congeners (`--exclude-unnamed-congeners`)

Many SILVA sequences have no species name — `uncultured bacterium`, `Fusobacterium sp.`, `uncultured Helicobacter sp.` — but are still placed in a genus. For a species-level probe, such sequences in the **target's own genus** may well be the target species, so counting them as off-target can make a probe look much less specific than it is. They often show up as `uncultured bacterium` in `blast_top_offtarget`.

With `--exclude-unnamed-congeners`, species-level rows keep all standard columns and add scores with those sequences set aside. The target genus is the most common SILVA genus among the target's sequences (so *E. coli* → `Escherichia-Shigella`). Unnamed sequences in other genera and named sister species still count as off-target. Treat the standard and `_excl_unnamed` values as worst and best case.

| Column | Description |
|---|---|
| `target_genus` | SILVA genus of the target species |
| `exact_unnamed_congener_hits` | Exact off-target hits that are unnamed and in `target_genus` |
| `exact_offtarget_excl_unnamed` | `exact_offtarget` without them |
| `exact_specificity_excl_unnamed` | `exact_specificity` without them (background also excludes them) |
| `blast_unnamed_congener_hits` | Same for BLAST hits (with `--blast`) |
| `blast_specificity_excl_unnamed` | `blast_specificity` without them |
| `blast_weighted_specificity_excl_unnamed` | `blast_weighted_specificity` without them |

`design --exclude-unnamed-congeners` uses the same rule and adds the matching `exact_*_excl_unnamed` and `exact_precision_excl_unnamed` columns, selecting on the `_excl_unnamed` values.

#### Coverage by sub-taxon (`--breakdown`)

`--breakdown LEVEL` adds a table showing, for each taxon in the input, how many of its sequences in each sub-taxon (for example each species of a genus) each probe hits, and the pooled coverage of all that taxon's probes together:

```bash
bac-probes validate \
    ./silva_db/SILVA_138.1_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz \
    fusobacterium_probes.csv \
    --breakdown species --min-seqs 5 \
    --output fusobacterium_validation.tsv
```

The table is printed and, with `--output`, written next to it as `<output>_breakdown.tsv` (here `fusobacterium_validation_breakdown.tsv`). Columns: `taxa`, `level`, `subtaxon`, `total_seqs`, one column per probe (`k1`, `k2`, … in input order within each taxon, each as `hits/total (pct%)`), and `pool_coverage`. The first row per taxon is `ALL`; sub-taxa with fewer than `--min-seqs` sequences are hidden.

## Key options

### `bac-probes design`

| Option | Default | Description |
|---|---|---|
| `--level` | `genus` | Level for a single taxon, or for CSV rows without a level |
| `--target-fasta` | — | Single taxon only: target sequences from a FASTA (e.g. NCBI) instead of SILVA |
| `-k / --kmer-size` | `32` | Probe length (bp) |
| `--min-conservation` | `0.02` | Keep candidates present in at least this fraction of target sequences |
| `--gc-min / --gc-max` | `0.35 / 0.65` | GC content bounds |
| `--max-homopolymer` | `5` | Reject probes with homopolymer runs ≥ this length |
| `--merge-overlapping` | off | Merge overlapping k-mers into longer probes |
| `--exclude-unnamed-congeners` | off | Species level: don't count unnamed same-genus sequences as off-target when selecting |
| `--blast / --no-blast` | `--no-blast` | BLAST near-match scoring; select on `blast_weighted_specificity` |
| `--blast-db` | `silva_db/silva` | BLAST database prefix (created by `build-db`) |
| `--blast-identity` | `85.0` | Min % identity for a BLAST hit to count |
| `--blast-max-target-seqs` | `50000` | BLAST cap; must exceed the largest taxon's sequence count |
| `--blast-pool` | `120` | Candidates BLASTed per taxon |
| `--select / --no-select` | `--select` | Select a probe set; `--no-select` writes candidate tables only |
| `--n-probes` | `10` | Maximum probes selected per taxon |
| `--pool-target` | `0.50` | Pooled coverage to aim for |
| `--spec-floors` | `1,0.999,0.995,0.99,0.98,0.95,0.9` | Specificity floors to try |
| `--candidates / --no-candidates` | `--candidates` | Write full candidate tables |
| `--threads` | `4` | Worker processes and BLAST threads |
| `-o / --output-dir` | `bac_probes_design` | Output directory |

### `bac-probes validate`

| Option | Default | Description |
|---|---|---|
| `--no-header` | off | CSV has no header row; assumes columns are `taxa`, `level`, `sequences` |
| `--top-offtargets` | `5` | Number of top off-target organism names to report per probe |
| `--blast / --no-blast` | `--no-blast` | Run BLAST to score near-match off-target binding |
| `--blast-db` | `silva_db/silva` | BLAST database prefix (created by `build-db`) |
| `--blast-identity` | `85.0` | Min % identity for a BLAST hit to count |
| `--blast-max-target-seqs` | `10000` | BLAST `max_target_seqs` cap |
| `--threads` | `4` | CPU threads for BLAST |
| `--exclude-unnamed-congeners` | off | Species level: also report specificity without unnamed sequences in the target's genus (see above) |
| `-b / --breakdown` | — | Also report coverage per sub-taxon at this level (see above) |
| `--min-seqs` | `5` | With `--breakdown`: hide sub-taxa with fewer sequences |
| `--output` | — | Write results to TSV (default: print to stdout only) |

## Interpreting results

### Choosing a conservation threshold

`--min-conservation` only decides which candidates are kept; coverage comes from selection pooling several probes. The default `0.02` keeps rare but very specific probes in play, which is what lets selection reach a pooled target with highly specific probes. Raising it (for example to `0.25` for large phyla) shrinks candidate tables and speeds up the run, at the cost of fewer options. A pool of several probes at different sites is usually better than one perfect probe.

### BLAST `blast_capped` warning

When `blast_capped = True`, BLAST returned exactly `--blast-max-target-seqs` hits. This happens when the taxon has more matching sequences than the cap, because the hit list is truncated before off-target near-misses are reached. Raise the cap above the taxon's sequence count (the `design` default, 50,000, covers every family in SILVA). In practice:

- `blast_weighted_specificity` is still a useful **relative** ranking — probes with fewer off-targets consistently score higher even under the cap.
- For absolute specificity estimates on large taxa, rely on `exact_specificity` (which is exhaustive) and confirm near-miss results with a targeted BLAST at `--blast-identity 97`.

### BLAST identity and k-mer length

For 32-mers, BLAST identity thresholds have a discrete effect:

| Identity | Max mismatches in 32 bp | Note |
|---|---|---|
| 97% | 0 | Equivalent to exact match; all off-target hits have `specificity_score = 0.0` (no mismatch protection) |
| 90% | 3 | Catches near-misses; position weighting meaningful |
| 85% (default) | 4 | Broader near-miss coverage |

The default 85% provides the widest near-miss coverage. For a focused check on 1–2 mismatch near-matches, re-run BLAST at `--blast-identity 90`.

### Taxon-level specificity expectations

Specificity varies greatly by how phylogenetically isolated the target genus is:

| Target | Exact off-targets | Notes |
|---|---|---|
| *Fusobacterium* (genus) | **0** | Phylogenetically isolated (Fusobacteriota phylum); excellent probes at genus level |
| *Bacteroides* (genus) | ~10–15 | Bacteroidota; very specific, top off-target is Lachnospiraceae |
| *Porphyromonas* (genus) | ~80–200 | Cross-reacts with *Tannerella* (same family); position weighting helps |
| *Prevotella* (genus) | 600–900 | Indistinguishable from *Alloprevotella* at 16S k-mer level; genus-specific probes not feasible |

**General guidance**: Genera with few close relatives in SILVA (isolated phyla or families) yield much cleaner probes than genera embedded in large, diverse orders like Bacteroidales or Lachnospirales. When `exact_offtarget` is high, check `blast_top_offtarget` — if the main cross-reactor is a genus you need to distinguish from, a different marker gene may be required.

## Memory and runtime

- **Target sequences** for every taxon are held in memory; **background sequences** are streamed from disk.
- `design` reads SILVA twice per run regardless of how many taxa it covers: once for target sequences, once for the exact off-target scan (spread over `--threads` processes). Expect ~10–15 min for a list of genera or families on a multi-core machine.
- BLAST takes ~1–2 s per candidate with 14 threads; `--blast-pool 120` over 25 families is ~45–60 min.
- k-mer counting for large phyla (>100,000 sequences) can use several GB per worker. Lower `--threads` if memory is tight.

## Development

```bash
pip install -e ".[dev]"
pytest tests/
```

Tests cover k-mer extraction, conservation scoring, BTOP parsing, position-weighted binding scores, taxonomy and species matching, BLAST result classification, probe selection, and the `design` and `validate` commands on small synthetic databases — all without requiring SILVA.
