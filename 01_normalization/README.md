# 01 — normalization

## Purpose

Rewrite every VCF into a **single canonical representation** so that biologically
identical variants written differently by different tools stop looking different. This
is the prerequisite for any honest set comparison (step 03).

## Why it matters

The same variant can be written several valid ways: an indel placed at different
positions inside a repeat, a multi-allelic site packed into one row vs split across
rows, an MNP vs adjacent SNPs, untrimmed REF/ALT. Without normalisation, `bcftools
isec` counts these as *different* sites and the methods look far more discordant than
they are. In this study normalisation **raises the fraction of SNP sites shared by all
six methods from 68.8% to 71.9%** — i.e. it improves cross-method concordance. (The
union *count* goes slightly up, 1.71 M → 1.79 M, because `-m -any` splits multi-allelic
SNP sites into more bi-allelic records; read the shared-fraction, not the union size.)

## What `bcftools norm -m -any` does

| operation | effect |
|-----------|--------|
| **left-align** (`-f REF`) | shift indels to the leftmost equivalent position in a repeat |
| **trim / parsimony** | strip shared prefixes/suffixes so REF/ALT are minimal |
| **split multi-allelics** (`-m -any`) | one row per ALT allele → all records become bi-allelic |
| **`--atom-overlaps .`** | resolve overlapping deletions deterministically |
| **`-c w`** | warn (don't crash) if a record's REF disagrees with the FASTA |

## Inputs

- 6 VCFs in `Pnapi_chr10_FR997704/`
- reference FASTA from step 00

## Run

```bash
bash normalize.sh          # normalises all 6 in parallel, then tabix-indexes
```

Core command:
```bash
bcftools norm -f REF -m -any -c w --atom-overlaps . -Oz -o out.norm.vcf.gz in.vcf.gz
bcftools index -t out.norm.vcf.gz
```

## Outputs

- `compare_out/norm/<base>.norm.vcf.gz` (+ `.tbi`) — one per input VCF

## ⚠️ Important downstream caveat

`-m -any` (splitting multi-allelics) is correct for **site comparison and
annotation**, but it **inflates π** if the split VCF is then fed to `pixy`, because
pixy counts the pairwise differences of each split bi-allelic record separately. See
`06_norm_vs_raw/` for the ~12% effect we measured. **Rule of thumb: normalise for
isec/annotation; use the original un-split all-sites VCFs for π/dₓy/Fₛₜ.**

## Adapting

Edit `VCFDIR`, `REF`, `OUT` at the top of `normalize.sh`. The script skips files that
already have a `.tbi`, so it is safe to re-run.
