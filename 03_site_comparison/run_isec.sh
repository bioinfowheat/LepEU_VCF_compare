#!/usr/bin/env bash
# Tier 2 — SNP-site set intersections with bcftools isec.
#
# Three intersections are produced:
#   isec_raw   : 6-way, raw VCFs
#   isec_norm  : 6-way, normalised VCFs
#   isec_all12 : 12-way, the 6 raw + 6 normalised together (for direct raw-vs-norm
#                comparison in one UpSet / Jaccard)
#
# Each writes a sites.txt whose 5th column is a presence/absence BITMASK, one bit
# per input file, IN THE ORDER THE FILES ARE LISTED ON THE COMMAND LINE. That
# bitmask is what the UpSet (upset_*.Rmd) and Jaccard (jaccard.py) steps consume,
# so the file order below DEFINES the column order downstream — keep them in sync.
#
# We compare SNP sites only (variant positions), so we first make variant-only
# views of every VCF (the study VCFs are all-sites; isec on 13 Mb x 12 would be
# enormous and meaningless for a *site*-set comparison).
set -euo pipefail

VCFDIR=Pnapi_chr10_FR997704
NORMDIR=compare_out/norm
OUT=compare_out
mkdir -p "$OUT/var_raw" "$OUT/var_norm" "$OUT/logs"

# method order (defines bitmask column order)
METHODS=(bwa_cohort bwa_hwe bwa_nohwe ngm_cohort ngm_hwe ngm_nohwe)

# 1. variant-only (SNP) views
extract(){ # $1 in.vcf.gz  $2 out.vcf.gz
  [ -s "$2".tbi ] && return
  bcftools view -v snps -Oz -o "$2" "$1"; bcftools index -t "$2"
}
for m in "${METHODS[@]}"; do
  extract "$VCFDIR/Pnapi_${m/_//}"*.vcf.gz "" 2>/dev/null || true
done
# (explicit, robust loop — handles the exact filenames in this study)
for v in "$VCFDIR"/*.vcf.gz;        do extract "$v" "$OUT/var_raw/$(basename "$v" .vcf.gz).var.vcf.gz"; done
for v in "$NORMDIR"/*.norm.vcf.gz;  do extract "$v" "$OUT/var_norm/$(basename "$v" .vcf.gz).var.vcf.gz"; done

# 2a. 6-way isec, raw
bcftools isec -p "$OUT/isec_raw"  -Oz -n+1 --collapse none \
  $(for m in "${METHODS[@]}"; do echo "$OUT"/var_raw/Pnapi_"${m%_*}"_bcftools_"${m#*_}"_*.var.vcf.gz; done) \
  2> "$OUT/logs/isec_raw.log"

# 2b. 6-way isec, normalised
bcftools isec -p "$OUT/isec_norm" -Oz -n+1 --collapse none \
  $(for m in "${METHODS[@]}"; do echo "$OUT"/var_norm/Pnapi_"${m%_*}"_bcftools_"${m#*_}"_*.norm.var.vcf.gz; done) \
  2> "$OUT/logs/isec_norm.log"

# 2c. 12-way isec: 6 raw THEN 6 norm (order matters for the bitmask -> labels)
bcftools isec -p "$OUT/isec_all12" -Oz -n+1 --collapse none \
  $(for m in "${METHODS[@]}"; do echo "$OUT"/var_raw/Pnapi_"${m%_*}"_bcftools_"${m#*_}"_*.var.vcf.gz; done) \
  $(for m in "${METHODS[@]}"; do echo "$OUT"/var_norm/Pnapi_"${m%_*}"_bcftools_"${m#*_}"_*.norm.var.vcf.gz; done) \
  2> "$OUT/logs/isec_all12.log"

for d in isec_raw isec_norm isec_all12; do
  echo "$d: $(wc -l < "$OUT/$d/sites.txt") union sites"
done
