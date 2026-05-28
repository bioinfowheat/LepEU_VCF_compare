#!/usr/bin/env bash
# Folder-wise VCF site intersection -> UpSet plot.
#
# Usage: bash run_upset.sh <folder_of_vcfs>
#
# Acts on every *.vcf.gz in the given folder. Outputs land in ./compare_out/.
# Column 5 of compare_out/isec/sites.txt is a presence/absence bitmask, one bit
# per input file in the order recorded in compare_out/file_order.txt. That order
# is load-bearing — the UpSet step labels columns from it.

set -euo pipefail

VCFDIR=${1:?usage: bash run_upset.sh <folder_of_vcfs>}
OUT=${OUT:-compare_out}
mkdir -p "$OUT/var" "$OUT/isec" "$OUT/logs"

# 1. index any VCF missing a .tbi
for v in "$VCFDIR"/*.vcf.gz; do
  [ -s "$v.tbi" ] || bcftools index -t "$v"
done

# 2. SNP-only views (site comparison only cares about variant positions)
for v in "$VCFDIR"/*.vcf.gz; do
  o="$OUT/var/$(basename "$v" .vcf.gz).var.vcf.gz"
  bcftools view -v snps -Oz -o "$o" "$v"
  bcftools index -t "$o"
done

# 3. freeze the file order — defines bitmask column order downstream
ls "$OUT"/var/*.var.vcf.gz | tee "$OUT/file_order.txt"

# 4. N-way intersection (union, presence/absence bitmask in column 5 of sites.txt)
bcftools isec -p "$OUT/isec" -Oz -n+1 --collapse none \
  $(cat "$OUT/file_order.txt") 2> "$OUT/logs/isec.log"

echo "union sites: $(wc -l < "$OUT/isec/sites.txt")"

# 5. UpSet plot
Rscript -e '
  library(UpSetR)
  out <- Sys.getenv("OUT", "compare_out")
  sites <- read.table(file.path(out, "isec/sites.txt"), sep="\t")[, 5]
  labs  <- basename(readLines(file.path(out, "file_order.txt")))
  mat   <- do.call(rbind, lapply(strsplit(as.character(sites), ""), as.integer))
  colnames(mat) <- labs
  pdf(file.path(out, "upset.pdf"), width=10, height=6)
  upset(as.data.frame(mat), nsets=ncol(mat), order.by="freq")
  dev.off()
' OUT="$OUT"

echo "wrote $OUT/upset.pdf"
