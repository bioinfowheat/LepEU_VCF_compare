#!/usr/bin/env bash
# Tier 3 (annotation): VEP failed on this Mac (dyld symbol error in the bioconda
# ensembl-vep build), so we use `bcftools csq` instead — same GFF + FASTA inputs,
# same taxonomy of consequences.
#
# FIX (v2): the previous version piped csq through `bcftools +split-vep`, which
# DROPPED every variant with no BCSQ tag (i.e. all intergenic / no-overlap SNPs).
# That made it look like ~98% of variants were intronic. We now use plain
# `bcftools query` so unannotated variants are kept and emitted as '.', which the
# parser buckets as 'intergenic / no annotated feature'.
set -euo pipefail
export PATH=/Users/stockholmbutterflylab/miniforge3/envs/vcfcompare/bin:$PATH
ROOT=/Users/stockholmbutterflylab/sbl_claudecode/VCF_compare
REF=$ROOT/genome/Pnapi_FR997704.fa
GFF=$ROOT/genome/Pnapi_FR997704.sorted.gff.gz
OUT=$ROOT/compare_out
mkdir -p "$OUT/vep"

for v in "$OUT"/var_norm/*.var.vcf.gz; do
  base=$(basename "$v" .var.vcf.gz)
  out=$OUT/vep/${base}.csq.tsv
  echo "[csq] $base"
  bcftools csq -f "$REF" -g "$GFF" -p s --ncsq 16 -Ou "$v" 2> "$OUT/logs/csq.${base}.log" \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/BCSQ\n' - \
    > "$out"
  total=$(wc -l < "$out")
  inter=$(awk -F'\t' '$5=="."||$5==""{c++} END{print c+0}' "$out")
  echo "  done $base : total=$total  intergenic(.)=$inter"
done
echo "csq annotation (v2, intergenic-preserving) done"
