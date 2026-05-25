#!/usr/bin/env bash
# =============================================================================
# run_all.sh — run the whole comparison on ANY set of VCFs.
# Edit config.sh, then:  bash run_all.sh
# =============================================================================
# Steps: normalise -> per-file stats -> variant-only views -> N-way intersection
#        -> UpSet (R) + Jaccard (python) -> annotation (csq) -> pixy (pi/Fst/dxy)
# Annotation and pixy are skipped automatically if GFF / POPMAP are not set.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/config.sh"
export PATH="$ENV_BIN:$PATH"

n=${#VCFS[@]}
[ "${#LABELS[@]}" -eq "$n" ] || { echo "ERROR: VCFS (${#VCFS[@]}) and LABELS (${#LABELS[@]}) differ in length"; exit 1; }
echo ">>> Comparing $n VCFs"

mkdir -p "$OUTDIR"/{norm,stats_raw,stats_norm,var_raw,var_norm,isec,vep,pixy,plots,logs}

# labels.txt records the order; the UpSet + Jaccard scripts read it.
printf "%s\n" "${LABELS[@]}" > "$OUTDIR/labels.txt"

# ---------------------------------------------------------------------------
# 1. normalise, stats, variant-only SNP views (per VCF)
# ---------------------------------------------------------------------------
for i in "${!VCFS[@]}"; do
  v="${VCFS[$i]}"; lab="${LABELS[$i]}"
  echo ">>> [$((i+1))/$n] $lab"
  norm="$OUTDIR/norm/${lab}.norm.vcf.gz"
  bcftools norm -f "$REF" -m -any -c w --atom-overlaps . -Oz -o "$norm" "$v" 2> "$OUTDIR/logs/norm.${lab}.log"
  bcftools index -t "$norm"
  bcftools stats -s - "$v"    > "$OUTDIR/stats_raw/${lab}.stats"
  bcftools stats -s - "$norm" > "$OUTDIR/stats_norm/${lab}.stats"
  bcftools view -v snps -Oz -o "$OUTDIR/var_raw/${lab}.var.vcf.gz"  "$v";    bcftools index -t "$OUTDIR/var_raw/${lab}.var.vcf.gz"
  bcftools view -v snps -Oz -o "$OUTDIR/var_norm/${lab}.var.vcf.gz" "$norm"; bcftools index -t "$OUTDIR/var_norm/${lab}.var.vcf.gz"
done

# ---------------------------------------------------------------------------
# 2. N-way intersection on NORMALISED SNP sites  (the headline comparison)
#    sites.txt column 5 = an N-bit presence mask, one bit per VCF, in LABELS order.
# ---------------------------------------------------------------------------
varfiles=(); for lab in "${LABELS[@]}"; do varfiles+=("$OUTDIR/var_norm/${lab}.var.vcf.gz"); done
bcftools isec -p "$OUTDIR/isec" -Oz -n+1 --collapse none "${varfiles[@]}" 2> "$OUTDIR/logs/isec.log"
echo ">>> intersection: $(wc -l < "$OUTDIR/isec/sites.txt") union SNP sites"

# ---------------------------------------------------------------------------
# 3. UpSet plot (R / UpSetR) + pairwise Jaccard (python). Both read labels.txt.
# ---------------------------------------------------------------------------
# NB: rmarkdown::render chdir's to the Rmd folder, so we pass ABSOLUTE paths and
# also set knit_root_dir; that way relative paths would resolve sanely too.
Rscript -e ".libPaths('$R_LIBS_DIR'); Sys.setenv(RSTUDIO_PANDOC='$PANDOC_DIR'); \
  rmarkdown::render('$HERE/upset_generic.Rmd', \
    params=list(sites='$OUTDIR/isec/sites.txt', labels='$OUTDIR/labels.txt', outdir='$OUTDIR/plots'), \
    knit_root_dir='$PWD', output_file='upset.html', output_dir='$OUTDIR', quiet=TRUE)"
python "$HERE/jaccard_generic.py" --sites "$OUTDIR/isec/sites.txt" --labels "$OUTDIR/labels.txt" --outdir "$OUTDIR/plots"
echo ">>> wrote $OUTDIR/plots/upset.png and jaccard.png"

# ---------------------------------------------------------------------------
# 4. annotation (optional — needs a bgzipped+tabix GFF)
# ---------------------------------------------------------------------------
if [ -n "${GFF:-}" ] && [ -f "$GFF" ]; then
  for lab in "${LABELS[@]}"; do
    bcftools csq -f "$REF" -g "$GFF" -p s --ncsq 16 -Ou "$OUTDIR/var_norm/${lab}.var.vcf.gz" 2> "$OUTDIR/logs/csq.${lab}.log" \
      | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/BCSQ\n' - > "$OUTDIR/vep/${lab}.csq.tsv"
  done
  echo ">>> annotation written to $OUTDIR/vep/"
else
  echo ">>> GFF not set/found — skipping annotation"
fi

# ---------------------------------------------------------------------------
# 5. pixy pi/Fst/dxy (optional — needs popmap + ALL-SITES VCFs)
#    NB: pixy needs variant+invariant sites. Feed the ORIGINAL all-sites VCFs,
#    not the variant-only views, and prefer the un-split (raw) VCFs (see ../05_popgen).
# ---------------------------------------------------------------------------
if [ -n "${POPMAP:-}" ] && [ -f "$POPMAP" ]; then
  for i in "${!VCFS[@]}"; do
    lab="${LABELS[$i]}"; v="${VCFS[$i]}"
    mkdir -p "$OUTDIR/pixy/${lab}"
    pixy --stats pi fst dxy --vcf "$v" --populations "$POPMAP" \
         --window_size "$WINDOW" --n_cores "$NCORES" \
         --output_folder "$OUTDIR/pixy/${lab}" --output_prefix "$lab" \
         2> "$OUTDIR/logs/pixy.${lab}.log" || echo "    pixy failed for $lab (is it an all-sites VCF?)"
  done
  echo ">>> pixy written to $OUTDIR/pixy/"
else
  echo ">>> POPMAP not set/found — skipping pixy"
fi

echo ">>> DONE. Key outputs:"
echo "    $OUTDIR/isec/sites.txt        (N-bit presence mask per SNP)"
echo "    $OUTDIR/plots/upset.png       (multi-way intersection)"
echo "    $OUTDIR/plots/jaccard.png     (pairwise similarity)"
