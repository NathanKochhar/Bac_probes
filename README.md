# bac-probes

Find bacterial 16S rRNA k-mers that are **highly conserved within a taxonomic group** (phylum, genus, or species) and **specific** — i.e., unlikely to bind outside of it.

Designed for **Xenium in-situ probe design**, with position-weighted off-target scoring that down-penalises near-matches whose mismatches fall near the probe centre (where mismatches are most destabilising).

Uses the [SILVA 138.1 NR99](https://www.arb-silva.de/) SSU reference database (~510,000 curated 16S sequences) as the sequence source.

## How it works

1. **Conservation scan** — all k-mers (default 32-mers) from the target taxon's 16S sequences are counted. A k-mer's conservation score is the fraction of target sequences that contain it exactly.
2. **Exact off-target scan** — candidate k-mers are checked against every non-target SILVA sequence for exact matches. This is fast and memory-efficient because only the candidates (not all background k-mers) are held in memory.
3. **BLAST near-match scoring** *(optional)* — top candidates are BLASTed against the full SILVA database using `blastn-short`. Each hit is classified as on-target or off-target using SILVA taxonomy. Off-target hits are scored by **effective binding probability**: mismatches near the probe centre (which destabilise hybridisation most) contribute less to the off-target score than terminal mismatches.

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

### Output columns

| Column | Description |
|---|---|
| `kmer` | K-mer sequence (same strand as SILVA reference) |
| `probe` | Reverse complement — order this as the in-situ probe sequence |
| `gc_pct` | GC content (%) |
| `conservation_pct` | % of target sequences containing this k-mer exactly |
| `target_seqs` | Number of target sequences used |
| `exact_offtarget` | Non-target SILVA sequences containing the k-mer exactly |
| `exact_bg_seqs` | Total non-target sequences scanned |
| `exact_specificity` | `1 − (exact_offtarget / exact_bg_seqs)` |
| `blast_total_hits` | Total BLAST hits (target + off-target) |
| `blast_target_hits` | BLAST hits within the target taxon |
| `blast_offtarget_hits` | Raw count of off-target BLAST hits |
| `blast_weighted_offtarget` | Sum of effective binding scores for off-target hits; down-weights near-misses with central mismatches |
| `blast_specificity` | `blast_target_hits / blast_total_hits` (raw) |
| `blast_weighted_specificity` | `blast_target_hits / (blast_target_hits + blast_weighted_offtarget)` — **primary ranking metric** |
| `blast_capped` | `True` when BLAST hit the `max_target_seqs` limit; specificity estimates may be unreliable for large target taxa |
| `blast_top_offtarget` | Most frequent off-target taxon in BLAST results |

**Ranking:** results are sorted by `blast_weighted_specificity` (desc) then `conservation_pct` (desc).

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

Useful for finding the exact spelling required by the `find` command.

### 4. Find specific k-mers

```bash
# Genus-level (default)
bac-probes find \
    ./silva_db/SILVA_138.1_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz \
    Bacteroides \
    --level genus \
    --min-conservation 0.60 \
    --output bacteroides_probes.tsv

# Phylum-level, higher conservation threshold, no BLAST
bac-probes find \
    ./silva_db/SILVA_138.1_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz \
    Bacteroidota \
    --level phylum \
    --min-conservation 0.80 \
    --no-blast \
    --output bacteroidota_probes.tsv

# Species-level, custom k-mer size, custom target sequences from NCBI
bac-probes find \
    ./silva_db/SILVA_138.1_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz \
    "Fusobacterium nucleatum" \
    --level species \
    --kmer-size 28 \
    --target-fasta my_ncbi_sequences.fasta \
    --output fusobacterium_probes.tsv
```

### 5. Inspect per-species coverage

After finding genus-level probes, `bac-probes coverage` shows which species within the genus are covered by each k-mer — and what pool coverage looks like if you use multiple probes.

```bash
# Auto-pick top 2 probes from a find results TSV
bac-probes coverage \
    ./silva_db/SILVA_138.1_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz \
    Fusobacterium \
    --level genus --breakdown species \
    --from-tsv fusobacterium_probes.tsv --top-n 2

# Or specify k-mers explicitly
bac-probes coverage \
    ./silva_db/SILVA_138.1_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz \
    Fusobacterium \
    --level genus --breakdown species \
    --kmer GATGGGGAAGCCAGCTTACTGGACAGATACTG \
    --kmer ATGCAGGGCTCAACTCTGTATTGCGTTGGAAA \
    --output fusobacterium_coverage.tsv
```

Output columns: `subtaxon`, `total_seqs`, one column per k-mer (`k1_…`, `k2_…`), and `pool_coverage` (sequences hit by any k-mer).

### 6. Validate existing probes

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

Add `--blast` to also score near-match off-target binding via BLAST (same position-weighted scoring as `bac-probes find`). When using `--blast` with `validate`, the probe sequence you supply is treated as the full k-mer — its length determines the position-weighting window. Mismatches in BLAST hits are scored by their distance from the centre of your probe sequence: a mismatch at the centre contributes close to 1.0 (high protection against off-target binding), while terminal mismatches contribute close to 0.0 (minimal protection). This means supplying a truncated or padded sequence will shift the weighting incorrectly, so pass the exact probe length you intend to use.

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

Additional columns added by `--blast` (same definitions as `bac-probes find`): `blast_total_hits`, `blast_target_hits`, `blast_offtarget_hits`, `blast_weighted_offtarget`, `blast_specificity`, `blast_weighted_specificity`, `blast_top_offtarget`, `blast_capped`.

## Key options

### `bac-probes find`

| Option | Default | Description |
|---|---|---|
| `--level` | `genus` | Taxonomic level: `domain` / `phylum` / `class` / `order` / `family` / `genus` / `species` |
| `-k / --kmer-size` | `32` | K-mer length (bp) |
| `--min-conservation` | `0.80` | Minimum fraction of target sequences that must contain the k-mer |
| `--gc-min / --gc-max` | `0.35 / 0.65` | GC content bounds for valid k-mers |
| `--max-homopolymer` | `5` | Reject k-mers with homopolymer runs ≥ this length |
| `--blast / --no-blast` | `--blast` | Enable/disable BLAST near-match scoring |
| `--blast-db` | `silva_db/silva` | BLAST database prefix (created by `build-db`) |
| `--blast-identity` | `85.0` | Min % identity for a BLAST hit to count |
| `--top-n-blast` | `500` | Number of top candidates to submit to BLAST (ranked by exact specificity) |
| `--threads` | `4` | CPU threads for BLAST |
| `--target-fasta` | — | Supply target sequences from a plain FASTA (e.g. NCBI download) instead of SILVA |
| `--output` | `probes.tsv` | Output TSV path |

### `bac-probes coverage`

| Option | Default | Description |
|---|---|---|
| `--level` | `genus` | Taxonomic level for the target name |
| `--breakdown` | `species` | Sub-level to break coverage down by |
| `--kmer SEQ` | — | K-mer(s) to test (repeatable) |
| `--from-tsv FILE` | — | Auto-pick top k-mers from a `find` results TSV |
| `--top-n` | `3` | Number of k-mers to pick from `--from-tsv` |
| `--min-seqs` | `5` | Minimum sequences to include a breakdown row |
| `--output` | — | Write table to TSV (default: stdout only) |

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
| `--output` | — | Write results to TSV (default: print to stdout only) |

## Interpreting results

### Choosing a conservation threshold

| Level | Recommended threshold | Rationale |
|---|---|---|
| Phylum | `0.80` (default) | Phyla are diverse; only universal regions pass |
| Genus | `0.60–0.75` | Variable-region probes with partial coverage |
| Species | `0.50–0.70` | Within-species diversity can be high |

Lower conservation finds more candidates at the cost of some sequences being missed. A pool of 2–4 probes covering different loci is often better than one perfect probe.

### BLAST `blast_capped` warning

When `blast_capped = True`, BLAST returned exactly `max_target_seqs` (1000) hits. This happens for large target taxa (e.g. genus with >1000 matching sequences) because the hit list is truncated before including all target sequences. In practice:

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

- **Target sequences** are held in memory (~10–100 MB for genus/species-level targets).
- **Background sequences** are streamed from disk in pass 2 — no full-database RAM required.
- SILVA NR99 has ~510,000 sequences. A typical `find` run takes **10–20 min** on a laptop (most time is the two FASTA scan passes over 1 GB compressed). BLAST on 500 candidates adds ~2 min.
- Phylum-level targets with >10,000 sequences may use 1–4 GB RAM for k-mer counting.

## Development

```bash
pip install -e ".[dev]"
pytest tests/
```

Tests cover k-mer extraction, conservation scoring, BTOP parsing, position-weighted binding scores, taxonomy parsing, and BLAST result classification — all without requiring the SILVA database.
