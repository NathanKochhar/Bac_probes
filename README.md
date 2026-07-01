# bac-probes

Find bacterial 16S rRNA k-mers that are **highly conserved within a taxonomic group** (phylum, genus, or species) and **specific** — i.e., unlikely to bind outside of it.

Designed for **Xenium in-situ probe design**, with position-weighted off-target scoring that down-penalises near-matches whose mismatches fall near the probe centre (where mismatches are most destabilising).

Uses the [SILVA 138.1 NR99](https://www.arb-silva.de/) SSU reference database (~510,000 curated 16S sequences) as the sequence source.

## How it works

1. **Conservation scan** — all k-mers (default 32-mers) from the target taxon's 16S sequences are counted. A k-mer's conservation score is the fraction of target sequences that contain it exactly.
2. **Exact off-target scan** — candidate k-mers are checked against every non-target SILVA sequence for exact matches. This is fast and memory-efficient because only the candidates (not all background k-mers) are held in memory.
3. **BLAST near-match scoring** *(optional)* — top candidates are BLASTed against the full SILVA database using `blastn-short`. Each hit is classified as on-target or off-target using SILVA taxonomy. Off-target hits are scored by **effective binding probability**: mismatches near the probe centre (which destabilise hybridisation most) contribute less to the off-target score than terminal mismatches.

### Position-weighted scoring (Xenium-specific)

Each off-target BLAST hit is assigned an **effective binding score** in [0, 1]:

```
weight(pos) = sin(π × (pos + 0.5) / k)²
binding_score = max(0, 1 − Σ weight(mismatch_positions))
```

- A perfect-match hit scores **1.0** (certain off-target binding).
- A single central mismatch (pos ≈ k/2) has weight ≈ 1.0 → binding score ≈ **0.0** (probe unlikely to hybridise).
- A single terminal mismatch (pos = 0 or k−1) has weight ≈ 0.002 → binding score ≈ **0.998** (probe may still bind).

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

## Key options

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
| 97% | 0 | Equivalent to exact match; all off-target hits have `binding_score = 1.0` |
| 90% | 3 | Catches near-misses; position weighting meaningful |
| 85% (default) | 4 | Broader near-miss coverage |

The default 85% provides the widest near-miss coverage. For a focused check on 1–2 mismatch near-matches, re-run BLAST at `--blast-identity 90`.

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
