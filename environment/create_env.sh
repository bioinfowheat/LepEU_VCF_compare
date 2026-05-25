#!/usr/bin/env bash
# Create the conda/mamba environment used by every step of this pipeline.
#
# All command-line tools (bcftools, vcftools, pixy, bedtools, samtools, htslib,
# pandoc) plus the Python plotting stack live in one env. R + UpSetR are handled
# separately (see note at the bottom) because UpSetR is installed from CRAN.
#
# Usage:  bash create_env.sh
set -euo pipefail

ENV=vcfcompare

mamba create -n "$ENV" -y \
  -c bioconda -c conda-forge --strict-channel-priority \
  bcftools vcftools pixy bedtools samtools tabix htslib \
  pandoc \
  python=3.11 pandas matplotlib seaborn numpy scipy upsetplot

# ---------------------------------------------------------------------------
# NOTES
#
# 1. ensembl-vep: we originally requested ensembl-vep here, but on macOS / Apple
#    Silicon the bioconda build fails to load (a dyld "missing symbol called"
#    error in one of its perl bundles). We therefore use `bcftools csq` for
#    annotation (step 04), which takes the same GFF + FASTA and yields the same
#    consequence taxonomy. On Linux, or via the official VEP Docker image, VEP
#    proper should work — add `ensembl-vep perl-bioperl` to the create line.
#
# 2. R + UpSetR (used in step 03 for the UpSet plots) are installed from system R
#    + CRAN, not conda, so versions are decoupled from the env:
#
#       R_LIBS=$HOME/R_libs_vcfcompare
#       mkdir -p "$R_LIBS"
#       Rscript -e '.libPaths("'"$R_LIBS"'"); install.packages(c("UpSetR","data.table","rmarkdown"), repos="https://cloud.r-project.org")'
#
#    rmarkdown needs pandoc; we install pandoc into the conda env above and point
#    R at it with  Sys.setenv(RSTUDIO_PANDOC="<env>/bin")  (see step 03).
#
# 3. PLINK2 (optional, for PCA/IBS/ROH — not run here) is not in the bioconda
#    channel set above; install separately if needed.
# ---------------------------------------------------------------------------
echo "Environment '$ENV' created. Activate with:  mamba activate $ENV"
