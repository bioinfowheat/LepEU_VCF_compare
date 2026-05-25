#!/usr/bin/env bash
# Tier 1 — per-file QC with bcftools stats, on raw AND normalised VCFs.
#
# Produces one .stats file per VCF (parsed later by lib/analyze_and_plot.py into
# the SNP / indel / multi-allelic count table and ts/tv, het/hom, depth summaries).
#
# Inputs : the 6 raw VCFs + the 6 normalised VCFs (step 01)
# Outputs: stats_raw/<base>.stats , stats_norm/<base>.stats
set -euo pipefail

VCFDIR=Pnapi_chr10_FR997704          # raw VCFs
NORMDIR=compare_out/norm             # normalised VCFs (from step 01)
OUT=compare_out
mkdir -p "$OUT/stats_raw" "$OUT/stats_norm" "$OUT/logs"

# raw
for v in "$VCFDIR"/*.vcf.gz; do
  base=$(basename "$v" .vcf.gz)
  bcftools stats -s - "$v" > "$OUT/stats_raw/${base}.stats" 2>> "$OUT/logs/stats_raw.log" &
done
wait

# normalised
for v in "$NORMDIR"/*.norm.vcf.gz; do
  base=$(basename "$v" .vcf.gz)
  bcftools stats -s - "$v" > "$OUT/stats_norm/${base}.stats" 2>> "$OUT/logs/stats_norm.log" &
done
wait

echo "Per-file stats written to $OUT/stats_raw and $OUT/stats_norm"
# Quick peek:  grep '^SN' compare_out/stats_raw/<file>.stats
