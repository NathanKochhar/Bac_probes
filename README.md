# bac-probes

Find bacterial 16S rRNA k-mers that are **highly conserved within a taxonomic group** (phylum, genus, or species) and **specific** — i.e., unlikely to bind outside of it.

Uses the [SILVA](https://www.arb-silva.de/) NR99 SSU reference database as the sequence source and BLAST+ for near-match off-target scoring.

## How it works

1. **Conservation scan** — all k-mers (default 32-mers) from the target taxon's 16S sequences are counted. A k-mer's conservation score is the fraction of target sequences that contain it.
2. **Exact off-target scan** — candidate k-mers (those above the conservation threshold) are checked against every non-target SILVA sequence for exact matches. This is fast and memory-efficient because only the candidates (not all background k-mers) are held in memory.
3. **BLAST near-match scoring** *(optional)* — top candidates are BLASTed against the full SILVA BLAST database using `blastn-short`. This catches sequences that might cross-hybridize even with 1–4 mismatches. Each hit is classified as on-target or off-target using SILVA taxonomy.

### Output columns

| Column | Description |
|---|---|
| `kmer` | K-mer sequence (same strand as SILVA reference) |
| `probe` | Reverse complement — use this as the actual probe/primer sequence |
| `gc_pct` | GC content (%) |
| `conservation_pct` | % of target sequences containing this k-mer exactly |
| `target_seqs` | Number of target sequences used |
| `exact_offtarget` | Non-target SILVA sequences containing the k-mer exactly |
| `exact_bg_seqs` | Total background sequences scanned |
| `exact_specificity` | `1 − (exact_offtarget / exact_bg_seqs)` |
| `blast_total_hits` | Total BLAST hits across all SILVA sequences |
| `blast_target_hits` | BLAST hits within the target taxon |
| `blast_offtarget_hits` | BLAST hits outside the target taxon |
| `blast_specificity` | `target_hits / total_hits` from BLAST |
| `blast_top_offtarget` | Most common off-target taxon in BLAST results |

## Installation

```bash
pip install -e .
# BLAST+ is required for --blast scoring:
conda install -c bioconda blast
```

## Workflow

### 1. Download SILVA

```bash
bac-probes download --output-dir ./silva_db
```

Downloads `SILVA_138.1_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz` (~3 GB) into `./silva_db/`.
Gap characters (`.` / `-`) from the multiple sequence alignment are stripped automatically on read.

### 2. Build BLAST database *(optional but recommended)*

```bash
bac-probes build-db ./silva_db/SILVA_138.1_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz \
    --db-dir ./silva_db
```

Strips alignment gaps, then runs `makeblastdb`. Takes a few minutes. Required for `--blast` scoring.

### 3. Browse taxon names

```bash
bac-probes list-taxa ./silva_db/SILVA_138.1_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz \
    --level genus --top 50
```

Useful for finding the exact spelling expected by the `find` command.

### 4. Find specific k-mers

```bash
# Genus-level (default)
bac-probes find ./silva_db/SILVA_138.1_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz \
    Streptococcus \
    --level genus \
    --output streptococcus_probes.tsv

# Phylum-level, stricter conservation, no BLAST
bac-probes find ./silva_db/SILVA_138.1_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz \
    Firmicutes \
    --level phylum \
    --min-conservation 0.90 \
    --no-blast \
    --output firmicutes_probes.tsv

# Species-level, custom k-mer size
bac-probes find ./silva_db/SILVA_138.1_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz \
    "Escherichia coli" \
    --level species \
    --kmer-size 28 \
    --output ecoli_probes.tsv
```

## Key options

| Option | Default | Description |
|---|---|---|
| `--level` | `genus` | Taxonomic level: domain / phylum / class / order / family / genus / species |
| `-k / --kmer-size` | `32` | K-mer length in bp |
| `--min-conservation` | `0.80` | Minimum fraction of target sequences that must contain the k-mer |
| `--gc-min / --gc-max` | `0.35 / 0.65` | GC content bounds for valid probes |
| `--max-homopolymer` | `5` | Reject k-mers with homopolymer runs ≥ this length |
| `--blast / --no-blast` | `--blast` | Enable/disable BLAST near-match scoring |
| `--blast-identity` | `85.0` | Min % identity for a BLAST hit to count as off-target |
| `--top-n-blast` | `500` | Number of top candidates to submit to BLAST |
| `--threads` | `4` | BLAST CPU threads |

## Memory and runtime notes

- **Target sequences** are held in memory. Species/genus-level targets typically use < 100 MB.
- **Background sequences** are streamed from disk in a second pass — no full-database RAM required.
- SILVA NR99 has ~510,000 sequences. A typical `find` run takes 10–30 min on a laptop (most time is the two FASTA scan passes). BLAST adds a few additional minutes for the top candidates.
- For phylum-level targets with many sequences, conservation scanning may use 1–4 GB RAM.
