# 04 — genomic-context annotation

## Purpose

Classify every SNP by the genomic feature it falls in (intergenic, intronic, UTR,
synonymous, missense, splice, stop/start…), per method, so you can ask whether the
methods agree on **where in the genome** their SNPs are, and feed the per-site
consequence into the integrative enrichment analysis (step 07).

## Tool: `bcftools csq` (VEP substitute)

We use `bcftools csq` rather than Ensembl VEP because VEP's bioconda build fails on
macOS (see `environment/README.md`). `csq` takes the same GFF + FASTA and produces an
equivalent consequence taxonomy. Swap VEP back in on Linux/Docker if preferred.

## Inputs

- normalised variant-only VCFs (`compare_out/var_norm/*.norm.var.vcf.gz`, from step 03)
- reference FASTA + sorted/indexed GFF (step 00)

## Run

```bash
bash annotate_csq.sh
```

Core command:
```bash
bcftools csq -f REF -g chr.sorted.gff.gz -p s --ncsq 16 -Ou in.vcf.gz \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/BCSQ\n' - > out.csq.tsv
```

| flag | meaning |
|------|---------|
| `-g` | bgzipped+tabix GFF3 with gene models |
| `-p s` | phase handling: treat as un-phased, take a conservative per-sample call |
| `--ncsq 16` | allow up to 16 consequences per record |

## ⚠️ Pitfall this script fixes

An earlier version piped `csq` through `bcftools +split-vep`, which **silently dropped
every variant with no consequence tag — i.e. all intergenic SNPs** — making it look
like ~98% of SNPs were intronic. **`bcftools query` keeps them**, emitting `.` for
unannotated (intergenic) sites, which the parser then labels
*intergenic / no annotated feature*. With the fix: ~39% intergenic, ~60% intronic,
~1% coding — as expected for a genome.

**`intronic` ≠ `intergenic`:** intronic = inside an annotated gene's intron (csq emits
an `intron` tag); intergenic = outside every transcript (no tag, shown as `.`). They
are distinct classes.

## Outputs

- `compare_out/vep/<base>.csq.tsv` — one row per SNP: `CHROM POS REF ALT BCSQ`

The BCSQ field looks like `missense|GENE|TRANSCRIPT|protein_coding|+|29Y>29I|28791T>A`;
lines beginning `@POS` are phased-neighbour pointers; `.` = no overlap (intergenic).

Plotting (stacked bars: all-classes + a coding-only zoom with intron+intergenic
removed) is done by `lib/analyze_and_plot.py` → `vep_consequence_by_method.png` and
`vep_consequence_coding_only.png`.

## Interpretation (this test case)

- ~39% intergenic, ~60% intronic, ~1% coding — **near-identical across all six
  methods**, so annotation conclusions are robust to mapper/filter choice.
- Of coding/splice SNPs: ~55% synonymous, ~32% splice_region, ~13% missense.

## Adapting

Edit `REF`, `GFF`, input dir at the top of `annotate_csq.sh`. For VEP instead:
`vep --gff chr.sorted.gff.gz --fasta REF --vcf -i in.vcf.gz -o out.vep.vcf`.
