# 00 — reference preparation

## Purpose

Produce a **single-chromosome reference FASTA + GFF whose contig name matches the
VCFs**. Three tools downstream (`bcftools norm`, `bcftools csq`, `pixy`) require the
reference contig name to be identical to the `#CHROM` value in the VCFs.

## The problem it solves

The VCFs label the chromosome with the **ENA/INSDC** accession `FR997704.1`, but the
downloaded genome (`GCF_905475465.1_ilPieNapi1.2`) uses **RefSeq** accessions
(`NC_062xxx.1`). *Pieris napi* chr 10 = `NC_062243.1`, confirmed by the GFF attribute
`chromosome=10`. We extract that record and **rename its header to `FR997704.1`** so
everything lines up. (Renaming the FASTA is cleaner than carrying a `bcftools annotate
--rename-chrs` alias through every step.)

## Inputs

- `GCF_905475465.1_ilPieNapi1.2_genomic.fna` — full genome multi-FASTA
- `GCF_905475465.1_ilPieNapi1.2_genomic.gff` — full GFF3

## Run

```bash
bash prepare_reference.sh
```

## Key commands explained

```bash
# extract one contig and rename its FASTA header to the VCF's contig name
samtools faidx genome.fna NC_062243.1 | sed 's/^>NC_062243.1.*/>FR997704.1/' > chr.fa
samtools faidx chr.fa

# keep only that contig's GFF rows, rewrite column 1.
# FS=OFS='\t' is ESSENTIAL — otherwise awk turns spaces inside the
# attribute column (col 9) into tabs and breaks the GFF.
awk -F'\t' -v OFS='\t' '$1=="NC_062243.1"{$1="FR997704.1"; print}' genome.gff

# csq/VEP need the GFF sorted, bgzipped, tabix-indexed
( grep '^#' chr.gff; grep -v '^#' chr.gff | sort -k1,1 -k4,4n ) | bgzip > chr.sorted.gff.gz
tabix -p gff chr.sorted.gff.gz
```

## Outputs

- `Pnapi_FR997704.1.fa` (+ `.fai`) — renamed single-chromosome reference
- `Pnapi_FR997704.1.gff` — chromosome-only GFF
- `Pnapi_FR997704.1.sorted.gff.gz` (+ `.tbi`) — for `bcftools csq` / VEP

## How to find your chromosome's RefSeq id

```bash
grep -E "chromosome=" genome.gff | grep "chromosome=10;" | head -1   # -> NC_062243.1
# or look up the assembly report on NCBI to map ENA <-> RefSeq accessions
```

## Adapting

Edit `CHR_REFSEQ`, `CHR_VCFNAME`, `CHR_LEN`, and the two genome paths at the top of
`prepare_reference.sh`. If your VCF and reference already share contig names, skip this
step entirely.
