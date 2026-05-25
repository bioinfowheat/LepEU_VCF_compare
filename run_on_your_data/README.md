# run_on_your_data — re-run the whole comparison on new VCFs

This folder is the **"bring your own VCFs"** entry point. The per-step folders
(`00_…`–`07_…`) document each method in depth and are wired to the demo dataset;
the scripts here are **generalised** so you can compare *any* number of VCFs
(2, 6, 12, …) from any combination of mappers / callers / filters by editing a
single config file.

```
run_on_your_data/
├── config.sh            # ← edit this (your VCF paths + labels + reference + options)
├── run_all.sh           # runs the whole pipeline using config.sh
├── upset_generic.Rmd    # UpSet plot for any N call sets (RMarkdown / UpSetR)
└── jaccard_generic.py   # pairwise Jaccard heatmap for any N call sets
```

---

## 1. Prerequisites (once)

```bash
# tools: see ../environment/README.md
bash ../environment/create_env.sh && mamba activate vcfcompare
# R + UpSetR (for the RMD):
mkdir -p ~/R_libs_vcfcompare
Rscript -e '.libPaths("~/R_libs_vcfcompare"); install.packages(c("UpSetR","data.table","rmarkdown"), repos="https://cloud.r-project.org")'
```

Two things that trip people up:

1. **Contig names must match.** Every VCF's `#CHROM` value must equal a contig name
   in your reference FASTA. If they differ (e.g. RefSeq `NC_…` vs ENA `FR…`), fix it
   first with `../00_reference_prep/prepare_reference.sh`.
2. **For pixy you need ALL-SITES VCFs** (variant + invariant), and ideally the
   un-split (raw) ones — see `../05_popgen/README.md`. The intersection / UpSet /
   Jaccard / annotation steps work on ordinary variant VCFs.

## 2. Configure

Edit `config.sh`. The only thing that really matters is the two parallel arrays:

```bash
VCFS=(   /data/bwa.vcf.gz   /data/ngm.vcf.gz   /data/minimap2.vcf.gz )
LABELS=( bwa                ngm                minimap2              )
```

`VCFS` and `LABELS` must be the **same length and same order**. The label order is
the single source of truth for the bit/column order in every plot — pick short,
legend-friendly names.

## 3. Run

```bash
bash run_all.sh
```

It runs, per VCF: `bcftools norm` → `bcftools stats` → variant-only views; then the
N-way `bcftools isec`; then the UpSet (R) and Jaccard (Python); then (if a GFF /
popmap are configured) `bcftools csq` annotation and `pixy`. Annotation and pixy are
skipped automatically when those inputs aren't set.

Outputs land in `$OUTDIR` (default `./compare_out`):

| path | what |
|------|------|
| `isec/sites.txt` | one row per union SNP site; **column 5 = N-bit presence mask** |
| `plots/upset.png` | the multi-way intersection |
| `plots/jaccard.png` + `.tsv` | pairwise Jaccard similarity |
| `vep/<label>.csq.tsv` | per-site consequence (if GFF set) |
| `pixy/<label>/…` | π / Fₛₜ / dₓy windows (if popmap set) |
| `labels.txt` | the label order, consumed by the UpSet + Jaccard steps |

---

## 4. The intersection of all the SNP calls — in detail

This is the heart of the comparison. `bcftools isec -n+1 --collapse none` takes your
variant-only VCFs **in the order listed in `LABELS`** and writes `isec/sites.txt`:

```
chr   pos    ref  alt   mask
chr1  1318   G    A     10110     <- present in callset 1, 3, 4 ; absent in 2, 5
chr1  1391   G    A     00100     <- private to callset 3
```

- **`mask` is an N-bit string**, one bit per VCF, in `LABELS` order (bit 1 = first
  label). `1` = the SNP was called by that VCF, `0` = it was not.
- This single column encodes the *entire* Venn structure. Everything downstream
  (UpSet, Jaccard) is just a different summary of this column.

Why we run it on the **normalised** VCFs: two callers can write the same variant
differently (multi-allelic packing, indel placement, MNP vs SNPs). `bcftools norm
-m -any` (step 01) puts them in one canonical form so `isec` compares like with like.
Run it on raw VCFs too if you want to *measure* how much disagreement is merely
representational (the demo does both).

Quick things you can do straight from `sites.txt`:

```bash
# how many SNPs are shared by ALL call sets? (mask of all 1s)
n=$(head -1 compare_out/labels.txt >/dev/null; wc -l < compare_out/labels.txt)
allones=$(printf '1%.0s' $(seq 1 "$n"))
awk -v m="$allones" '$5==m' compare_out/isec/sites.txt | wc -l

# the most common sharing patterns
cut -f5 compare_out/isec/sites.txt | sort | uniq -c | sort -rn | head

# SNPs private to call set #2 (only bit 2 set), for n=3 -> 010
awk '$5=="010"' compare_out/isec/sites.txt | wc -l
```

`bcftools isec` also writes `0000.vcf.gz`, `0001.vcf.gz`, … — the actual records
private to / shared by each input combination, if you want the variants themselves
(see `isec/README.txt` for which file is which).

---

## 5. The RMD (UpSetR) analysis — in detail

`upset_generic.Rmd` turns the bitmask into the UpSet plot. It is **parameterised**, so
you never edit the Rmd itself — `run_all.sh` passes everything in. To run it by hand
(e.g. to tweak the figure), call `rmarkdown::render` with `params`:

```bash
# Use ABSOLUTE paths: rmarkdown::render() chdir's into the Rmd's folder, so
# relative paths would be resolved from there. (run_all.sh handles this for you.)
OUT="$PWD/compare_out"
Rscript -e '
.libPaths("~/R_libs_vcfcompare")
Sys.setenv(RSTUDIO_PANDOC="'"$HOME"'/miniforge3/envs/vcfcompare/bin")   # rmarkdown needs pandoc
rmarkdown::render("upset_generic.Rmd",
  params = list(
    sites  = "'"$OUT"'/isec/sites.txt",
    labels = "'"$OUT"'/labels.txt",
    outdir = "'"$OUT"'/plots"),
  knit_root_dir = "'"$PWD"'",
  output_file = "upset.html", output_dir = "'"$OUT"'", quiet = TRUE)'
```

What the Rmd does, step by step (all in the `setup` chunk):

1. `readLines(params$labels)` → the call-set names, in order.
2. `data.table::fread(params$sites)` → the isec table; it asserts that every mask is
   exactly `length(labels)` bits (a guard against a labels/VCF-count mismatch).
3. `strsplit(mask, "")` → a 0/1 **membership matrix**, one column per call set.
4. `UpSetR::upset(M, sets=rev(labels), nintersects=40, order.by="freq")` → the plot,
   also written to `plots/upset.png` for embedding elsewhere.

**Reading the UpSet plot:** each top bar = the number of SNP sites belonging to
*exactly* the combination of call sets marked by the filled dots beneath it; the
horizontal bars on the lower-left are each call set's total size. It is the readable
substitute for an N-way Venn diagram (which would need 2ᴺ−1 regions).

**Customising:**

- `nintersects` (default 40): how many intersection bars to show.
- `order.by`: `"freq"` (biggest intersections first) or `"degree"` (by how many sets
  participate).
- **Colour set bars by group.** If your labels fall into groups (say two mappers),
  pass a colour vector the same length as `labels`, in `rev(labels)` order:
  ```r
  sets.bar.color = rev(ifelse(grepl("^bwa", labels), "#2874a6", "#6c3483"))
  ```
  (The demo's `../03_site_comparison/upset_tier2.Rmd` and `upset_all12.Rmd` show this
  two-colour and a 12-set raw-vs-normalised version respectively — copy from there.)

`jaccard_generic.py` is the pairwise companion (every pair's
`|A∩B|/|A∪B|`); it reads the same `sites.txt` + `labels.txt` and writes
`plots/jaccard.png` + `jaccard.tsv`. It auto-drops the cell annotations once N > 14
so large grids stay legible.

---

## 6. Troubleshooting

| symptom | fix |
|---------|-----|
| `isec` warns about unindexed files | the scripts tabix-index automatically; if you supply your own, `bcftools index -t` each first |
| UpSet/Jaccard error about mask length | `labels.txt` count ≠ number of VCFs in the isec — they must match; re-run `run_all.sh` so both are regenerated together |
| `rmarkdown` can't find pandoc | set `RSTUDIO_PANDOC` to the conda env's `bin` (done for you in `run_all.sh`) |
| pixy fails | the VCF is probably variant-only; pixy needs an all-sites VCF (see `../05_popgen`) |
| empty / tiny intersections | check the VCFs really share a coordinate system (same reference, same contig names) |
