# 06 — what normalisation actually does (raw vs normalised)

## Purpose

Make the effect of `bcftools norm -m -any` **visible and quantified**, head-to-head
within each method, at two levels: record counts and the downstream π estimate.

## Inputs

- `compare_out/stats_{raw,norm}/*.stats` (step 02)
- `compare_out/pixy_{raw,norm}/**/*_pi.txt` (step 05)

## Run

```bash
python norm_vs_raw.py
```

## Outputs (→ `figures/`)

| figure | what it shows |
|--------|---------------|
| `norm_vs_raw_counts.png` | per-method SNP / multi-allelic / indel record counts, raw vs norm. Multi-allelic bars collapse to **0** after norm (all decomposed); SNP/indel counts rise — normalisation in action. |
| `norm_vs_raw_pi.png` | all methods pooled: windowed π raw (x) vs norm (y). Points sit **above** y=x — normalisation systematically inflates π by ~12%. |
| `pi_norm_vs_raw_per_method.png` | the same, one panel per method, each with slope, Spearman ρ, and mean % inflation. |
| `pi_norm_vs_raw_per_method.tsv` | the numbers behind the per-method panels. |
| `norm_vs_raw_dxy.png` | between-population dₓy raw (x) vs norm (y). Inflated by ~11.5% too — the **same** multi-allelic-splitting mechanism as π, since dₓy is also an absolute difference count. |

## The key finding

Splitting multi-allelic sites with `-m -any`:
- is **correct** for site comparison (step 03) and annotation (step 04);
- but **inflates π** when the split VCF is fed to pixy. Across methods the effect is
  remarkably uniform: **slope ≈ 1.12, ρ ≈ 0.998, +11–12% mean π, 97–99% of windows
  higher**. Neither mapper nor filter changes its size — it is purely a representation
  artefact: pixy counts the differences of each split bi-allelic record separately, so
  a multi-allelic position contributes more than once.

**Practical rule:** normalise (split) for isec/annotation; use the **un-split all-sites
VCFs** for π/dₓy/Fₛₜ.

## Interpretation summary

- π **shape** is robust to normalisation (ρ ≈ 0.998) — window ranking is unchanged.
- π **magnitude** is not (+~12%) — so absolute π values must be reported with the
  normalisation/multiallelic-handling choice stated.

## Adapting

Edit `--root` (default `compare_out`). The reference set of methods is the
`METHODS` list at the top of `norm_vs_raw.py`.
