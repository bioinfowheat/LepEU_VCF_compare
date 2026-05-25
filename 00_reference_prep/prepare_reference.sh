#!/usr/bin/env bash
# Prepare a single-chromosome reference FASTA + GFF that MATCH the VCF contig name.
#
# Problem this solves
# -------------------
# The VCFs in this study use the ENA/INSDC accession  FR997704.1  for Pieris napi
# chromosome 10, but the downloaded genome (GCF_905475465.1_ilPieNapi1.2) uses the
# RefSeq accessions  NC_062xxx.1 / NW_025920xxx.1.  Chromosome 10 is NC_062243.1
# (confirmed by the GFF tag  chromosome=10 ).  bcftools/pixy/VEP all require the
# reference contig name to match the VCF, so we extract chr10 and RENAME its header
# to FR997704.1.  We do the same for the GFF.
#
# Adapt for your data: change CHR_REFSEQ / CHR_VCFNAME / CHR_LEN and the file paths.
set -euo pipefail

GENOME_FA=GCF_905475465.1_ilPieNapi1.2_genomic.fna   # full multi-FASTA
GENOME_GFF=GCF_905475465.1_ilPieNapi1.2_genomic.gff  # full GFF3
CHR_REFSEQ=NC_062243.1        # contig name in the reference
CHR_VCFNAME=FR997704.1        # contig name used by the VCFs (what we rename TO)
CHR_LEN=13261761              # length of that chromosome (from the .fai)

OUT_FA=Pnapi_${CHR_VCFNAME}.fa
OUT_GFF=Pnapi_${CHR_VCFNAME}.gff
OUT_GFF_SORTED=Pnapi_${CHR_VCFNAME}.sorted.gff.gz

# 1. Extract the one chromosome and rename the FASTA header to the VCF contig name
samtools faidx "$GENOME_FA" "$CHR_REFSEQ" \
  | sed "s/^>${CHR_REFSEQ}.*/>${CHR_VCFNAME}/" \
  > "$OUT_FA"
samtools faidx "$OUT_FA"          # -> Pnapi_FR997704.1.fa.fai

# 2. Pull only this chromosome's GFF records and rewrite column 1.
#    NB: FS=OFS='\t' is essential so spaces inside the attribute column survive.
{
  echo "##gff-version 3"
  echo "##sequence-region ${CHR_VCFNAME} 1 ${CHR_LEN}"
  awk -F'\t' -v OFS='\t' -v a="$CHR_REFSEQ" -v b="$CHR_VCFNAME" \
      '$1==a{$1=b; print}' "$GENOME_GFF"
} > "$OUT_GFF"

# 3. Sort + bgzip + tabix-index the GFF (required by bcftools csq / VEP --gff)
(grep '^#' "$OUT_GFF"; grep -v '^#' "$OUT_GFF" | sort -k1,1 -k4,4n) \
  | bgzip > "$OUT_GFF_SORTED"
tabix -p gff "$OUT_GFF_SORTED"

echo "Wrote: $OUT_FA (+ .fai), $OUT_GFF, $OUT_GFF_SORTED (+ .tbi)"
