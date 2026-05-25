# lib — shared plotting engine & report builder

These two scripts are used by several steps and assemble the final report. They read
the outputs already written under `compare_out/` by steps 02–07.

## `analyze_and_plot.py`

The shared plotting engine. Parses the bcftools-stats, isec, csq, and pixy outputs and
writes most of the PNG/TSV figures:

- Tier 1 QC count tables
- Tier 2 site-overlap (Python `upsetplot`) + 6×6 Jaccard heatmaps
- Tier 3 consequence stacked bars (all classes **and** coding-only zoom)
- Tier 4 π & Fₛₜ tracks, Spearman correlation matrices, top-1% Fₛₜ outlier Jaccard
- Tier 5 norm-vs-raw counts + pooled π scatter

```bash
python lib/analyze_and_plot.py --root compare_out
```

> The R UpSetR plots (step 03) and the 12×12 Jaccard / per-method π grids
> (`v3_additions.py`, step 06/07) are produced separately and are higher-quality for
> those specific figures; `analyze_and_plot.py` also emits a simpler Python `upsetplot`
> version (`isec_sites_raw.png`) kept for a quick view.

## `v3_additions.py`

The 12×12 all-VCF Jaccard and the per-method raw-vs-norm π grid (also available as the
standalone `03_site_comparison/jaccard.py` and `06_norm_vs_raw/norm_vs_raw.py`).

```bash
python lib/v3_additions.py
```

## `make_report.py`

Assembles every figure + summary table under `compare_out/plots/` into one
self-contained HTML report (images embedded as base64), with the narrative, method
notes, and interpretation for each tier.

```bash
python lib/make_report.py        # -> reports/VCF_comparison_results_v3.html
```

## Typical end-to-end order

```bash
# after steps 01–05 have populated compare_out/:
python lib/analyze_and_plot.py
python lib/v3_additions.py
python 07_integrative/integrative_analysis.py
Rscript -e 'rmarkdown::render("03_site_comparison/upset_tier2.Rmd")'
Rscript -e 'rmarkdown::render("03_site_comparison/upset_all12.Rmd")'
python lib/make_report.py
```
