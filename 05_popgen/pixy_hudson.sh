#!/usr/bin/env bash
# Re-run pixy Fst with the Hudson estimator (--fst_type hudson) for all 6 methods
# x {raw, norm}, fst only. Output column will be avg_hudson_fst.
# WC (default) results already live in compare_out/pixy_{raw,norm}/.
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

run_one(){
  v=$1; regime=$2
  base=$(basename "$v" .vcf.gz)
  if [ "$regime" = raw ]; then in=$VCFDIR/$v
  else                          in=$OUT/norm/${base}.norm.vcf.gz; fi
  outd=$OUT/pixy_hudson_${regime}/${base}
  mkdir -p "$outd"
  [ -s "$outd/${base}_fst.txt" ] && { echo "skip $regime/$base"; return; }
  pixy --stats fst --fst_type hudson \
       --vcf "$in" --populations "$POPMAP" \
       --window_size 50000 --n_cores 4 \
       --output_folder "$outd" --output_prefix "$base" \
       2> "$OUT/logs/pixy_hudson.${regime}.${base}.log" \
    && echo "done $regime/$base" || echo "FAILED $regime/$base"
}
export -f run_one
export VCFDIR OUT POPMAP

# build the job list (12 jobs), run 4 in parallel (4 cores each)
for regime in raw norm; do
  for v in "${VCFS[@]}"; do echo "$v $regime"; done
done | xargs -P 4 -n 2 bash -c 'run_one "$0" "$1"'

echo "hudson pixy done"
