# VCF intersection → UpSet plot (folder-wise)

Simplified, copy-paste workflow that takes **every `*.vcf.gz` in one folder** and produces an UpSet plot of shared SNP sites. Distilled from `03_site_comparison/`, with the Jaccard step omitted — UpSet only.

## What it does

1. Indexes every VCF in the input folder (if not already indexed).
2. Builds SNP-only views (site comparison only cares about variant positions; intersecting all-sites VCFs is huge and meaningless for a site-set comparison).
3. Runs `bcftools isec -n+1` across all files — emits a `sites.txt` whose 5th column is a presence/absence **bitmask**, one bit per input file in the order passed on the command line.
4. Feeds that bitmask into UpSetR to draw the plot.

The order of files passed to `isec` defines the bitmask column order downstream. Step 3 of the script freezes that order in `file_order.txt` so the R step can label columns correctly.

## Requirements

- `bcftools` (with `tabix`)
- `R` with the `UpSetR` package

## Usage

```bash
bash run_upset.sh path/to/folder_of_vcfs
```

Outputs land in `compare_out/`:

- `compare_out/var/` — SNP-only views of each input
- `compare_out/isec/sites.txt` — union of sites, bitmask in column 5
- `compare_out/file_order.txt` — file order that defines the bitmask columns
- `compare_out/upset.pdf` — the plot

## One-liners (if you'd rather paste them by hand)

```bash
VCFDIR=my_vcfs ; OUT=compare_out ; mkdir -p $OUT/var $OUT/isec $OUT/logs

# index any VCF missing a .tbi
for v in $VCFDIR/*.vcf.gz; do [ -s $v.tbi ] || bcftools index -t $v; done

# SNP-only views
for v in $VCFDIR/*.vcf.gz; do o=$OUT/var/$(basename $v .vcf.gz).var.vcf.gz; bcftools view -v snps -Oz -o $o $v && bcftools index -t $o; done

# freeze file order (defines bitmask column order)
ls $OUT/var/*.var.vcf.gz | tee $OUT/file_order.txt

# N-way intersection (union of all sites, presence/absence bitmask)
bcftools isec -p $OUT/isec -Oz -n+1 --collapse none $(cat $OUT/file_order.txt) 2> $OUT/logs/isec.log

# UpSet plot from column 5 of sites.txt
Rscript -e 'library(UpSetR); m<-read.table("compare_out/isec/sites.txt",sep="\t")[,5]; labs<-basename(readLines("compare_out/file_order.txt")); mat<-do.call(rbind,lapply(strsplit(m,""),as.integer)); colnames(mat)<-labs; pdf("compare_out/upset.pdf",w=10,h=6); upset(as.data.frame(mat),nsets=ncol(mat),order.by="freq"); dev.off()'
```

## Notes

- If your VCFs are already variant-only, the SNP-extraction step is a no-op and you can skip it.
- `-n+1` means "include every site present in at least one file" (the union). Switch to `-n=N` for "present in exactly N files" or `-n~1111...` for an explicit bitmask filter.
- `--collapse none` requires identical REF/ALT for a site to be considered shared — strict matching, which is what you want for a method-comparison UpSet.
