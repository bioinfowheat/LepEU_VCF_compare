# LepEU_VCF_compare

A reproducible framework for **comparing a set of VCF files generated from the same
sequencing data by different methods**, at two levels:

1. **the SNPs called** (counts, set overlap, genomic context), and
2. **the downstream population-genomic measures** they produce (π, dₓy, Fₛₜ).

It was built for a *Pieris napi* (Lepidoptera) test case — chromosome FR997704.1
(= RefSeq `NC_062243.1`, chr 10), 14 samples in two populations of 7 (**A** vs **P**) —
but every step is parameterised so it can be re-pointed at any VCF set.

The six VCFs differ along two axes (same reads, same caller throughout):

| axis | levels | meaning |
|------|--------|---------|
| **mapper** | `bwa`, `ngm` | BWA-MEM vs NextGenMap read alignment |
| **filtering** | `cohort`, `hwe`, `nohwe` | baseline cohort genotyping; with a Hardy–Weinberg filter; without it |

So the grid is `{bwa,ngm} × {cohort,hwe,nohwe}` = 6 VCFs. Because reads and caller are
held constant, every difference is attributable to **read placement (mapper)** or the
**filtering decision** — which is exactly what this assessment isolates.

---

## Workflow

```
                        genome FASTA + GFF (full)
                                  │  00_reference_prep
                                  ▼
                  single-chromosome FASTA + GFF, renamed
                          to match the VCF contig
                                  │
        6 VCFs ───────────────────┼──────────────────────────────
          │                       │                               │
          │ 01_normalization      │                               │
          ▼                       ▼                               ▼
   6 normalised VCFs        02_per_file_qc                  (all-sites VCFs)
          │                  bcftools stats                       │
          ├──────────────┬───────────────┬──────────────┐        │
          ▼              ▼               ▼              ▼        ▼
   03_site_comparison  04_annotation  05_popgen     06_norm_vs_raw
   isec + UpSetR +     bcftools csq    pixy          counts + π raw-vs-norm
   Jaccard (6 & 12-way) (+ extracted   (π/dxy/Fst)
          │              GFF)             │
          └──────────────┴───────────────┴──────────────┐
                                                         ▼
                                              07_integrative
                                       (are discordant SNPs enriched in
                                        consequence classes? Ts/Tv quality)
                                                         │
                                                         ▼
                                      reports/VCF_comparison_results_v3.html
                                       (lib/make_report.py assembles everything)
```

## Repository layout

| folder | step | what it does |
|--------|------|--------------|
| `environment/` | — | conda/mamba env + R/UpSetR install |
| `00_reference_prep/` | 0 | extract one chromosome, fix the RefSeq-vs-ENA contig-name mismatch, prep GFF |
| `01_normalization/` | 1 | `bcftools norm -m -any` (left-align, trim, split multi-allelics) |
| `02_per_file_qc/` | 2 | `bcftools stats` per VCF (raw + normalised) |
| `03_site_comparison/` | 3 | `bcftools isec` → UpSetR plots (R) + pairwise Jaccard (Python), 6-way & 12-way |
| `04_annotation/` | 4 | `bcftools csq` genomic-context annotation using the GFF |
| `05_popgen/` | 5 | `pixy` windowed π, dₓy, Fₛₜ |
| `06_norm_vs_raw/` | 6 | direct raw-vs-normalised comparison (counts + π) |
| `07_integrative/` | 7 | enrichment of method-discordant SNPs + Ts/Tv quality signal |
| `lib/` | — | shared plotting engine (`analyze_and_plot.py`) + report builder (`make_report.py`) |
| `run_on_your_data/` | — | **generalised, config-driven re-run** on *your own* VCFs (any number); detailed walkthrough of the intersection + UpSetR (RMD) steps |
| `figures/` | — | all generated PNGs and summary TSVs |
| `reports/` | — | the assembled HTML report(s) |

Each step folder has its own `README.md` with the exact commands, parameter
explanations, inputs/outputs, and notes on adapting it to your own data.

## Running it on your own VCFs

If you just want to compare your own call sets (any number, any methods), go to
**[`run_on_your_data/`](run_on_your_data/)** — edit one config file and run one script.
That folder also contains the most detailed walkthrough of the SNP-set intersection
and the UpSetR (RMarkdown) analysis. The per-step folders below document the methods
in depth on the demo dataset.

## Reproducing everything / starting over on a new computer

**[`REPRODUCE.md`](REPRODUCE.md)** is the complete operator's manual: full environment
setup, the exact run order of every script, how the HTML report is assembled, how to
adapt it to a new VCF set, a ready-to-paste briefing for a **fresh Claude Code instance**,
and a "gotchas & lessons learned" section (contig-name fix, all-sites VCFs for pixy, the
multi-allelic π inflation, csq-vs-VEP, depth/goleft, etc.). Read it before reproducing.

## Quick start (demo dataset)

```bash
# 0. tools
bash environment/create_env.sh && mamba activate vcfcompare

# 1. reference (edit the contig names/paths at the top first)
bash 00_reference_prep/prepare_reference.sh

# 2. the pipeline (each script has paths/variables at the top to edit)
bash 01_normalization/normalize.sh
bash 02_per_file_qc/per_file_stats.sh
bash 03_site_comparison/run_isec.sh
python 03_site_comparison/jaccard.py
bash 04_annotation/annotate_csq.sh
bash 05_popgen/pixy.sh
python 06_norm_vs_raw/norm_vs_raw.py
python 07_integrative/integrative_analysis.py

# 3. UpSet plots (R) and the final report
Rscript -e 'rmarkdown::render("03_site_comparison/upset_tier2.Rmd")'
Rscript -e 'rmarkdown::render("03_site_comparison/upset_all12.Rmd")'
python lib/make_report.py
```

## Headline findings (this test case)

- **Mapper ≫ filter.** BWA-vs-NGM is the dominant axis of disagreement (UpSet + Jaccard);
  the cohort/hwe/nohwe filter is a minor perturbation. Within-mapper site-set Jaccard
  0.92–0.98, between-mapper 0.72–0.78.
- **Absolute diversity (π, dₓy) robust in shape, not magnitude.** Windowed π and dₓy
  both correlate >0.96–0.97 across methods, but splitting multi-allelics with
  `norm -m -any` before pixy inflates both by ~12% (see `06_norm_vs_raw/`).
  **Use un-split all-sites VCFs for π/dₓy/Fₛₜ.**
- **Fₛₜ is fragile.** Cross-method correlation 0.34–0.88; top-1% outlier
  (selection-candidate) windows overlap only 0–50% across methods — never trust a
  selection scan from a single VCF. As a ratio it cancels the normalisation inflation
  but amplifies the calls where methods disagree.
- **Annotation robust.** Genomic-context breakdown near-identical across all six methods.
- **Discordant SNPs are lower quality.** Method-private SNPs are depleted in
  synonymous/missense, enriched in apparent high-impact classes, and have a falling
  Ts/Tv toward the random floor — i.e. disagreement concentrates in the least
  trustworthy calls.

## Data note

The input VCFs (~560 MB each) and the genome are **not** included (see `.gitignore`).
The repo carries the code, the generated figures/tables, and the HTML report so the
results are inspectable without re-running. Point the scripts at your own data by
editing the path variables at the top of each script.

## Annotation tool note (VEP → bcftools csq)

Ensembl VEP from bioconda fails to load on macOS/Apple-Silicon (a dyld symbol error in
its perl bundle). We use `bcftools csq` instead — same GFF + FASTA inputs, same
consequence taxonomy. On Linux or via the VEP Docker image, VEP proper can be swapped
back in at step 04.
