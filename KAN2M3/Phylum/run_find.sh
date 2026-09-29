#!/usr/bin/env bash
# Run bac-probes find (no BLAST) for each phylum in KAN2M3-phylum.csv, in parallel.
cd /d/Repos/Bac_probes
export PYTHONIOENCODING=utf-8 PYTHONUTF8=1
SILVA=silva_db/SILVA_138.1_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz
OUT=KAN2M3/Phylum
for t in Firmicutes Proteobacteria Bacteroidota Actinobacteriota Fusobacteriota Campylobacterota; do
  .venv/Scripts/bac-probes.exe find "$SILVA" "$t" --level phylum \
    -k 32 --min-conservation 0.25 --no-blast \
    --output "$OUT/find_runs/${t}_find.tsv" > "$OUT/logs/${t}.log" 2>&1 &
done
wait
echo ALL_DONE
