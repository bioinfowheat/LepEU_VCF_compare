#!/usr/bin/env bash
# =============================================================================
# CONFIG — edit everything in this file, then run:  bash run_all.sh
# =============================================================================
# This is the ONLY file you should need to edit to run the comparison on a new
# set of VCFs. run_all.sh sources it.

# --- tools -------------------------------------------------------------------
# Directory that contains bcftools, pixy, tabix, bgzip, python (the conda env's bin).
ENV_BIN="$HOME/miniforge3/envs/vcfcompare/bin"
# Where UpSetR is installed (see ../environment/README.md), and where pandoc lives.
R_LIBS_DIR="$HOME/R_libs_vcfcompare"
PANDOC_DIR="$ENV_BIN"

# --- reference (must match the VCF contig names!) ----------------------------
# If your reference contig names differ from the VCFs, fix that first with
# ../00_reference_prep/prepare_reference.sh.
REF="/path/to/reference_chrom.fa"               # indexed (.fai); required
GFF="/path/to/reference_chrom.sorted.gff.gz"    # bgzipped+tabix; OPTIONAL (annotation). "" to skip.
POPMAP="/path/to/popmap.tsv"                     # sample<TAB>population; OPTIONAL (pixy). "" to skip.

# --- output + parameters -----------------------------------------------------
OUTDIR="$PWD/compare_out"
WINDOW=50000        # sliding-window size for pixy (bp)
NCORES=4

# --- THE VCFs TO COMPARE -----------------------------------------------------
# List every VCF, and give each a SHORT label. The two arrays must be the same
# length and in the SAME ORDER. The label order defines the bit/column order in
# every downstream plot (UpSet, Jaccard), so choose labels you can read in a legend.
#
# Works for ANY number of VCFs (2, 6, 12, ...). They can be from different
# mappers, callers, filtering settings — whatever you are comparing.
VCFS=(
  "/data/callsetA.vcf.gz"
  "/data/callsetB.vcf.gz"
  "/data/callsetC.vcf.gz"
)
LABELS=(
  "callsetA"
  "callsetB"
  "callsetC"
)
