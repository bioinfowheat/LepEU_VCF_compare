# 07 — integrative analysis

## Purpose

Move past *that* / *where* the methods disagree to the biologically pointed question:
**what kind of SNPs are the method-discordant ones?** Are they a random sample of the
genome, or enriched for particular consequence classes — and do they look like real
variants or artefacts?

Also contains the **one-vs-five windowed-π comparison** (a reference method plotted
against the other five) used in report section 4c.

## Inputs

- `compare_out/isec_norm/sites.txt` — the 6-bit sharing bitmask per SNP (step 03)
- `compare_out/vep/*.csq.tsv` — per-site consequence (step 04)
- `compare_out/pixy_{raw,norm}/**/*_pi.txt` — windowed π (step 05)

## Run

```bash
python integrative_analysis.py
```

## What it computes

### A. One-vs-five π scatter  (`pi_method_compare_{raw,norm}.png`)
Fix one reference method on the x-axis, plot each of the other five on the y-axis,
within a regime. Each panel reports a through-origin slope (1.0 = identical magnitude)
and Spearman ρ (1.0 = identical ranking). Change the reference with the `reference=`
argument to `pi_one_vs_rest()`.

### B. Consequence-class enrichment of discordant SNPs  (`enrichment_consequence.png`)
Every normalised union SNP is labelled by its **sharing category**:
`shared by all 6`, `BWA-only (all 3)`, `NGM-only (all 3)`, `singleton (1 method)`,
`other partial`. We cross that with each site's consequence class (most-severe, from
the union of the csq files) and show, per category, the **log2 fold-change of each
class's proportion vs the shared-by-all-6 core** (red = over-represented, blue = under).

### C. Ts/Tv quality signal  (`tstv_by_sharing.png`)
Transition/transversion ratio per sharing category, computed straight from REF/ALT.
Genuine SNPs sit well above the random expectation of **0.5**; false positives trend
toward it because sequencing/mapping errors have no transition bias.

## Outputs (→ `figures/`)

`pi_method_compare_raw.png`, `pi_method_compare_norm.png`,
`enrichment_consequence.png` (+ `.tsv`, `_fractions.tsv`),
`tstv_by_sharing.png` (+ `.tsv`), `sharing_category_counts.tsv`.

## Interpretation (this test case)

- **Enrichment:** discordant SNPs are strongly **depleted in synonymous (−2.6 to −3.1
  log2FC) and missense (−1.3 to −2.3)** — ordinary coding SNPs in well-mapped exons are
  the most reproducible. They are **enriched in apparent high-impact classes**
  (splice-donor/acceptor, stop/start gained-lost, non-coding RNA). ⇒ *the variants you'd
  get most excited about are disproportionately the method-unstable ones — confirm
  high-impact calls across mappers.* (Those classes have small N, so noisy; the
  synonymous/missense depletion is on hundreds of thousands of sites and is solid.)
- **Ts/Tv:** monotonic decline — shared 0.82 > BWA-only 0.76 > NGM-only 0.64 >
  singletons 0.62 — toward the 0.5 floor. The more method-restricted a SNP, the more it
  looks like noise. (Absolute values are low because this is genome-wide, dominated by
  intergenic/intronic sites; the *ordering* is the signal.)

## Next angles (not yet implemented)

- per-category minor-allele-frequency spectra (artefacts cluster at very low MAF);
- genotype-level concordance (hap.py / vcfeval) stratified by MAF;
- intersect discordant SNPs with a repeat / low-complexity mask to test the hypothesis
  that mapper-private SNPs sit in repetitive sequence.

## Adapting

Edit paths at the top; change the reference method in `__main__`; the sharing-category
logic is in `categorize()` and severity ordering in `SEVERITY`.
