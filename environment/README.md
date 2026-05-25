# environment — tools & versions

All command-line tools and the Python plotting stack live in one conda/mamba env
(`vcfcompare`). R + UpSetR are installed separately from CRAN (step 03 uses them).

## Install

```bash
bash create_env.sh
mamba activate vcfcompare
```

## What gets installed

| tool | used in | role |
|------|---------|------|
| bcftools (≥1.22) | 01,02,03,04 | norm, stats, isec, csq |
| htslib / tabix / bgzip | all | bgzip + index VCFs and the GFF |
| samtools | 00 | extract & index the chromosome FASTA |
| vcftools | 02 (optional) | per-sample missingness/het/depth |
| pixy | 05 | π, dₓy, Fₛₜ on all-sites VCFs |
| bedtools | (optional) | repeat / callable-region intersects |
| pandoc | 03 | required by R `rmarkdown::render` |
| python 3.11 + pandas/numpy/scipy/matplotlib/seaborn/upsetplot | 03,06,07,lib | analysis + figures + report |

## R + UpSetR (for the UpSet plots in step 03)

System R + CRAN, kept separate from the conda env:

```bash
R_LIBS=$HOME/R_libs_vcfcompare
mkdir -p "$R_LIBS"
Rscript -e '.libPaths("'"$R_LIBS"'"); install.packages(c("UpSetR","data.table","rmarkdown"), repos="https://cloud.r-project.org")'
```

`rmarkdown` needs pandoc; point R at the env's pandoc:
`Sys.setenv(RSTUDIO_PANDOC="<conda-env>/bin")` (the render commands in step 03 do this).

## Known issue: Ensembl VEP on macOS

The bioconda `ensembl-vep` build aborts on Apple Silicon with
`dyld: missing symbol called`. We substitute `bcftools csq` (step 04), which takes the
same GFF + FASTA and produces the same consequence classes. On Linux / VEP-Docker,
add `ensembl-vep perl-bioperl` to `create_env.sh` and use VEP at step 04 instead.
