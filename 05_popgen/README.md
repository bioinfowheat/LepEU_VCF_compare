# 05 — population-genomic estimators (pixy)

## Purpose

Compute the windowed population-genetic summaries that end up in a paper — **π**
(within-population nucleotide diversity), **dₓy** (between-population divergence), and
**Fₛₜ** (differentiation) — for every method, so you can see whether the choice of
mapper/filter changes the biological conclusion.

## Why pixy

`pixy` correctly handles **missing data and invariant sites** by using an all-sites
VCF (variant + invariant) as the denominator. This avoids the well-known upward bias
of variant-only π estimators. The study VCFs are already all-sites, so they can be fed
to pixy directly.

## Inputs

- the VCFs (raw and/or normalised — see caveat below), tabix-indexed
- `popmap.tsv` — two columns: `sample_id <TAB> population`

## Run

```bash
bash pixy.sh        # runs all 6 methods × {raw,norm}, 50 kb windows
```

Core command:
```bash
pixy --stats pi fst dxy \
     --vcf in.vcf.gz \
     --populations popmap.tsv \
     --window_size 50000 \
     --n_cores 4 \
     --output_folder out/ --output_prefix base
```

## Outputs

Per method per regime, in `compare_out/pixy_{raw,norm}/<base>/`:

| file | key columns |
|------|-------------|
| `<base>_pi.txt`  | `pop, window_pos_1, avg_pi, no_sites, count_diffs, count_comparisons` |
| `<base>_fst.txt` | `pop1, pop2, window_pos_1, avg_wc_fst, no_snps` |
| `<base>_dxy.txt` | `pop1, pop2, window_pos_1, avg_dxy` |

Plotting (overlaid per-method tracks; Spearman correlation matrices; top-1% Fₛₜ outlier
Jaccard) is done by `lib/analyze_and_plot.py`; method-vs-method π scatter grids by
`lib/v3_additions.py` / `07_integrative`.

## ⚠️ Caveat: don't pre-split multi-allelics for pixy

Running pixy on the **normalised** (`-m -any`) VCFs **inflates π by ~12%** versus the
raw VCFs, because splitting a multi-allelic site into per-allele rows makes pixy count
that position's pairwise differences more than once (verified: count_diffs +6.6%,
count_comparisons +0.5% on a test window). **Use the raw / un-split all-sites VCFs for
π/dₓy/Fₛₜ.** See `06_norm_vs_raw/` for the full demonstration.

## Interpretation (this test case)

- **π:** windowed π correlates 0.96–0.999 across all methods — the *shape* of the
  diversity landscape is robust. (Absolute magnitude is sensitive to normalisation.)
- **dₓy:** between-population divergence; correlates 0.97–0.999 across methods — robust
  in shape **like π, and unlike Fₛₜ**, because both π and dₓy are absolute difference
  counts rather than ratios. Subject to the same ~11.5% normalisation inflation as π
  (see `06_norm_vs_raw/`).
- **Fₛₜ:** cross-method correlation only 0.34–0.88, and top-1% outlier (selection-
  candidate) windows overlap just 0–50% across methods. **Fₛₜ scans are method-sensitive
  — never report selection candidates from a single VCF.**

Plots: `dxy_per_method_{raw,norm}.png` (overlaid tracks) and
`dxy_correlation_{raw,norm}.png` (cross-method Spearman) are produced by
`lib/analyze_and_plot.py` alongside the π and Fₛₜ figures.

**Combined genome scan:** `combined_tracks_{raw,norm}.png` stacks all three
measures (π mean of A,P → Fₛₜ → dₓy) sharing the x-axis, all six methods overlaid,
so diversity, differentiation and divergence can be read against each other along the
chromosome. In this test case π and dₓy crash together at ≈3 Mb and ≈7 Mb, coinciding
with the Fₛₜ peaks — the classic low-diversity/reduced-recombination signature.

## Fₛₜ estimator: Weir & Cockerham vs Hudson

pixy's default Fₛₜ is **Weir & Cockerham (1984)** (`avg_wc_fst`). The **Hudson
(1992; Bhatia et al. 2013)** estimator (`--fst_type hudson`, column `avg_hudson_fst`)
uses a ratio-of-averages formulation that is less biased under unequal sample sizes
and is often preferred for genome scans.

```bash
bash pixy_hudson.sh          # re-runs Fst only, with --fst_type hudson, into pixy_hudson_{raw,norm}/
python fst_hudson_compare.py # WC-vs-Hudson scatter, Hudson tracks + cross-method correlation
```

Core change vs `pixy.sh`:
```bash
pixy --stats fst --fst_type hudson --vcf in.vcf.gz --populations popmap.tsv \
     --window_size 50000 --output_folder out/ --output_prefix base
```

Outputs (→ `figures/`): `fst_wc_vs_hudson_{raw,norm}.png` (per-method scatter),
`fst_hudson_per_method_{raw,norm}.png` (Hudson tracks),
`fst_hudson_correlation_{raw,norm}.png` (Hudson cross-method Spearman), and
`fst_estimator_summary.tsv`.

**Findings (this test case):** Hudson is systematically slightly higher than WC at
differentiated windows (mean Hudson−WC +0.01 to +0.045), and the two **agree well for
the cohort/hwe call sets (Pearson 0.85–0.99) but diverge for `nohwe` and `ngm_cohort`
(Pearson 0.12–0.55)** — the divergence concentrated at the high-Fₛₜ outlier windows
that matter for selection scans. Hudson's cross-method correlation (~0.58–0.94) has a
higher floor than WC's (0.34–0.88), i.e. it is somewhat more robust to the calling
method. **Always state which estimator you used.**

## Mean read depth per window (QC overlay)

Average per-sample depth in the same 50 kb windows — the key QC track to read
*alongside* π/Fₛₜ/dₓy, because depth anomalies flag mapping artefacts.

```bash
bash depth_windows.sh     # vcftools --site-mean-depth on the raw all-sites VCFs
python depth_windows.py   # bins to 50 kb windows -> depth_per_window.tsv + depth_per_method.png
```

Core command:
```bash
vcftools --gzvcf raw_allsites.vcf.gz --site-mean-depth --out out
# -> out.ldepth.mean : CHROM POS MEAN_DEPTH VAR_DEPTH  (MEAN_DEPTH = mean FORMAT/DP over samples)
```

Run it on the **raw all-sites** VCFs so depth is sampled at ~every base (not just
SNPs). Depth is a property of the mapping, so the three filters share one curve per
mapper (the `.ldepth.mean` files are byte-identical within a mapper here).

`depth_windows.py` also writes `depth_per_window.tsv`, which `lib/analyze_and_plot.py`
picks up to add a **depth panel on top of the combined genome scan**
(`combined_tracks_{raw,norm}.png`).

**Findings (this test case):** BWA averages **8.94×** vs NGM **8.16×** (~9% more) —
the direct cause of BWA calling ~10% more SNPs (more reads → more callable sites). And
a **depth spike to ~22× at ≈7 Mb** (≈2.5× background) marks a collapsed-repeat /
paralog region; that window is exactly where Fₛₜ peaks and π/dₓy collapse, so the depth
track reclassifies that &ldquo;differentiation peak&rdquo; as a likely mapping
artefact. **Check the depth track before trusting any Fₛₜ outlier**; a practical filter
is to mask windows with depth &gt; ~2× the chromosome median (and very low-depth ones).

## Adapting

Edit `VCFDIR`, `NORMDIR`, `POPMAP`, `--window_size`, `--n_cores` in `pixy.sh`. Replace
the placeholder population labels in `popmap.tsv` with your real ones (they propagate
into output filenames and figure legends).
