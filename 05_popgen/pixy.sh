#!/usr/bin/env bash
# Tier 4: pixy pi/fst/dxy on all 6 VCFs, both raw and normalized.
# Runs sequentially across VCFs but pixy uses 4 cores per call.
set -euo pipefail
export PATH=/Users/stockholmbutterflylab/miniforge3/envs/vcfcompare/bin:$PATH
ROOT=/Users/stockholmbutterflylab/sbl_claudecode/VCF_compare
VCFDIR=$ROOT/Pnapi_chr10_FR997704
POPMAP=$ROOT/popmap.tsv
OUT=$ROOT/compare_out

VCFS=(
  Pnapi_bwa_bcftools_cohort_FR997704.1.vcf.gz
  Pnapi_bwa_bcftools_hwe_FR997704.1.vcf.gz
  Pnapi_bwa_bcftools_nohwe_FR997704.1.vcf.gz
  Pnapi_ngm_bcftools_cohort_FR997704.1.vcf.gz
  Pnapi_ngm_bcftools_hwe_FR997704.1.vcf.gz
  Pnapi_ngm_bcftools_nohwe_FR997704.1.vcf.gz
)

for v in "${VCFS[@]}"; do
  base=$(basename "$v" .vcf.gz)
  for regime in raw norm; do
    if [ "$regime" = raw ]; then in=$VCFDIR/$v
    else                          in=$OUT/norm/${base}.norm.vcf.gz; fi
    outd=$OUT/pixy_${regime}/${base}
    mkdir -p "$outd"
    if [ -s "$outd/${base}_pi.txt" ]; then echo "skip pixy ${regime}/${base}"; continue; fi
    echo "[pixy ${regime}] $base"
    pixy --stats pi fst dxy \
         --vcf "$in" \
         --populations "$POPMAP" \
         --window_size 50000 \
         --n_cores 4 \
         --output_folder "$outd" \
         --output_prefix "$base" 2> "$OUT/logs/pixy.${regime}.${base}.log" || \
         echo "FAILED ${regime}/${base}"
  done
done
echo "pixy done"
