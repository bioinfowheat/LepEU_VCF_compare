# 02 — per-file QC

## Purpose

Characterise each VCF on its own (before any cross-comparison) so a single anomalous
file is caught early. Run on **both** raw and normalised VCFs so the effect of
normalisation on the counts is visible.

## Inputs

- 6 raw VCFs + 6 normalised VCFs (step 01)

## Run

```bash
bash per_file_stats.sh
```

Core command:
```bash
bcftools stats -s - input.vcf.gz > input.stats     # -s - = per-sample stats too
```

## Outputs

- `compare_out/stats_raw/<base>.stats`
- `compare_out/stats_norm/<base>.stats`

Each `.stats` file contains (grep the prefixes):

| prefix | content |
|--------|---------|
| `SN` | summary numbers: # records, SNPs, indels, multi-allelic sites |
| `TSTV` | transition/transversion ratio |
| `PSC` | per-sample: nHet, nHom, mean depth, missingness |
| `AF`, `QUAL`, `DP` | allele-frequency, quality, depth distributions |

The SNP / indel / multi-allelic counts are parsed by `lib/analyze_and_plot.py` into the
QC table in the report, and by `06_norm_vs_raw/norm_vs_raw.py` into the raw-vs-norm
count bars.

## Interpretation (this test case)

- Within a mapper, the three filters give very similar counts; `nohwe` is largest
  (nothing removed).
- **BWA calls ~10% more SNPs than NGM** across all filters.
- After normalisation the *multi-allelic sites* count drops to 0 (all decomposed) and
  SNP record counts rise — the signature of `-m -any`.

## Optional per-sample diagnostics (vcftools)

```bash
vcftools --gzvcf in.vcf.gz --missing-indv --out base
vcftools --gzvcf in.vcf.gz --het         --out base
vcftools --gzvcf in.vcf.gz --depth       --out base
```

## Adapting

Edit `VCFDIR`, `NORMDIR`, `OUT` at the top of `per_file_stats.sh`.
