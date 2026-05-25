#!/usr/bin/env bash
# Per-site mean depth (across samples) for each method, via vcftools --site-mean-depth,
# on the RAW all-sites VCFs (so depth is sampled at ~every position, not just SNPs).
# Output: compare_out/depth/<method>.ldepth.mean  (CHROM POS MEAN_DEPTH VAR_DEPTH)
# These are binned to 50 kb windows by depth_windows.py.
set -euo pipefail
export PATH=/Users/stockholmbutterflylab/miniforge3/envs/vcfcompare/bin:$PATH
ROOT=/Users/stockholmbutterflylab/sbl_claudecode/VCF_compare
VCFDIR=$ROOT/Pnapi_chr10_FR997704
OUT=$ROOT/compare_out/depth
mkdir -p "$OUT"

one(){
  v=$1; base=$(basename "$v" .vcf.gz)
  [ -s "$OUT/${base}.ldepth.mean" ] && { echo "skip $base"; return; }
  vcftools --gzvcf "$v" --site-mean-depth --out "$OUT/${base}" 2> "$OUT/${base}.log"
  echo "done $base"
}
export -f one; export OUT

ls "$VCFDIR"/*.vcf.gz | xargs -P 6 -I{} bash -c 'one "{}"'
echo "site-mean-depth done"
