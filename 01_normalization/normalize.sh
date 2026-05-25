#!/usr/bin/env bash
set -euo pipefail
export PATH=/Users/stockholmbutterflylab/miniforge3/envs/vcfcompare/bin:$PATH
ROOT=/Users/stockholmbutterflylab/sbl_claudecode/VCF_compare
VCFDIR=$ROOT/Pnapi_chr10_FR997704
REF=$ROOT/genome/Pnapi_FR997704.fa
OUT=$ROOT/compare_out

norm_one(){
  v=$1; base=$(basename "$v" .vcf.gz); out=$OUT/norm/${base}.norm.vcf.gz
  if [ -s "$out".tbi ]; then echo "skip $base"; return; fi
  bcftools norm -f "$REF" -m -any -c w --atom-overlaps . \
       -Oz -o "$out" "$VCFDIR/$v" 2> "$OUT/logs/norm.${base}.log"
  bcftools index -t "$out"
  echo "done $base"
}
export -f norm_one
export REF OUT VCFDIR

# Run 6 in parallel (these are mostly I/O + single-thread bcftools, ok to overlap on a multi-core mac)
ls "$VCFDIR"/*.vcf.gz | xargs -I{} -P 6 -n 1 bash -c 'norm_one "$(basename {})"'
