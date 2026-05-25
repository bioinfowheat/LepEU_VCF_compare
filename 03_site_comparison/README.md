# 03 — site-set comparison (isec → UpSet + Jaccard)

## Purpose

Quantify **which SNP sites the methods agree on**, three ways:
1. multi-way intersection counts (`bcftools isec`),
2. **UpSet plots** (R / UpSetR) — the readable alternative to a 6-/12-way Venn,
3. **pairwise Jaccard** similarity matrices (Python).

Run for the 6 raw VCFs, the 6 normalised VCFs, and all 12 together.

## Inputs

- 6 raw + 6 normalised VCFs (steps 00–01)

## Step 3.1 — intersections

```bash
bash run_isec.sh
```

What it does:
- makes **variant-only (SNP) views** of every VCF (the study VCFs are all-sites; isec
  over 13 Mb × 12 would be meaningless for a *site*-set comparison);
- runs three `bcftools isec -n+1 --collapse none`:
  `isec_raw` (6-way), `isec_norm` (6-way), `isec_all12` (12-way: 6 raw THEN 6 norm).

**Output of interest:** each `sites.txt` has columns `CHR POS REF ALT MASK`, where
`MASK` is a presence/absence **bitstring**, one bit per input file *in the order they
were listed on the command line*. That order DEFINES the column order used by the
UpSet and Jaccard steps — keep them in sync if you change file order.

```
FR997704.1  1318  G  A  111110     # called by files 1-5, not file 6
```

## Step 3.2 — UpSet plots (R)

```bash
# point R at UpSetR's lib + the env's pandoc, then render:
Rscript -e '.libPaths("~/R_libs_vcfcompare");
            Sys.setenv(RSTUDIO_PANDOC="<conda-env>/bin");
            rmarkdown::render("upset_tier2.Rmd")'   # 6-way raw + 6-way norm
Rscript -e '... rmarkdown::render("upset_all12.Rmd")'  # 12-way raw+norm in one plot
```

- `upset_tier2.Rmd` → `upsetR_raw.png`, `upsetR_norm.png` (set bars coloured by mapper:
  blue = BWA, purple = NGM).
- `upset_all12.Rmd` → `upsetR_all12.png` (set bars coloured by regime: grey = raw,
  dark blue = norm — for reading the normalisation effect directly).

**How to read UpSet:** each top bar = number of SNP sites belonging to *exactly* the
method-combination marked by the filled dots beneath it; left bars = each method's
total set size.

## Step 3.3 — pairwise Jaccard (Python)

```bash
python jaccard.py            # writes jaccard_raw/norm/all12 .png + .tsv
```

**Jaccard(A,B) = |A ∩ B| / |A ∪ B|** = (sites called by both) / (sites called by
either). Range 0–1; scale-free, so a method that simply calls more variants doesn't
automatically look more similar to everything. We restrict to SNP sites, so this is
*site* concordance (agreement on **where** a SNP is), not genotype concordance.

## Outputs (→ `figures/`)

`upsetR_raw.png`, `upsetR_norm.png`, `upsetR_all12.png`,
`jaccard_raw.{png,tsv}`, `jaccard_norm.{png,tsv}`, `jaccard_all12.{png,tsv}`,
plus the compact `isec_sites_raw.png` (from the Python upsetplot in `lib/`).

## Interpretation (this test case)

- Largest UpSet bar = core shared by all (~1.18 M); next two largest are
  **mapper-specific** (~200 k all-BWA-only, ~96 k all-NGM-only). The filter axis makes
  only small private sets → **the mapper is the dominant source of method-specific SNPs.**
- Jaccard: within-mapper 0.92–0.98, between-mapper 0.72–0.78. Normalisation lifts every
  value ~0.01–0.02 but does **not** close the BWA-vs-NGM gap → that gap is real.
- 12×12 Jaccard: same-method raw-vs-norm cells ~0.93–0.95 (normalisation keeps the site
  set nearly intact); a BWA set is no closer to an NGM set after normalisation.

## Adapting

`METHODS` and the file globbing at the top of `run_isec.sh` define the set and order.
`jaccard.py` and the `.Rmd` files key off the same order — update all three together.
