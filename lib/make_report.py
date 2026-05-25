"""
Build VCF_comparison_results_v3.html — revised report.
Changes vs v2:
  - Tier 2: ADD a 12-way UpSet (6 raw + 6 normalised) for direct norm-vs-raw comparison
  - Tier 4: ADD one-vs-five windowed-pi scatter grid (method-vs-method)
  - New Tier 7: integrative analysis — are method-discordant SNPs enriched in
    particular consequence classes? + Ts/Tv quality signal by sharing category
"""
import os, json, base64, glob
from pathlib import Path
import pandas as pd

ROOT = Path("/Users/stockholmbutterflylab/sbl_claudecode/VCF_compare")
OUT = ROOT / "compare_out"
PLOTS = OUT / "plots"
TARGET = ROOT / "VCF_comparison_results_v3.html"

def img(name, alt="", maxw="100%"):
    p = PLOTS / name
    if not p.exists(): return f"<p><em>missing image: {name}</em></p>"
    b64 = base64.b64encode(p.read_bytes()).decode()
    return f'<img alt="{alt}" src="data:image/png;base64,{b64}" style="max-width:{maxw};border:1px solid #ccc;">'

def tsv_table(name):
    p = PLOTS / name
    if not p.exists(): return f"<p><em>missing table: {name}</em></p>"
    return pd.read_csv(p, sep="\t").to_html(index=False, classes="restable", border=0)

def isec_summary():
    rows=[]
    for regime in ("raw","norm"):
        s = OUT/f"isec_{regime}"/"sites.txt"
        if s.exists():
            n=sum(1 for _ in open(s)); a6=0
            for line in open(s):
                p=line.rstrip("\n").split("\t")
                if len(p)>=5 and p[4]=="111111": a6+=1
            rows.append({"regime":regime,"union SNP sites":f"{n:,}",
                         "shared by all 6":f"{a6:,}",
                         "% shared by all 6":f"{100*a6/n:.1f}%"})
    return pd.DataFrame(rows).to_html(index=False, classes="restable", border=0)

def stats_table():
    rows=[]
    for regime in ("raw","norm"):
        for f in sorted(glob.glob(str(OUT/f"stats_{regime}"/"*.stats"))):
            method=os.path.basename(f).replace(".stats","").replace("Pnapi_","").replace("_bcftools","").replace("_FR997704.1","")
            d={"method":method,"regime":regime}
            for line in open(f):
                if line.startswith("SN\t"):
                    p=line.rstrip("\n").split("\t"); k=p[2].rstrip(":")
                    if k in {"number of SNPs","number of indels","number of multiallelic sites"}:
                        d[k.replace("number of ","")]=int(float(p[3]))
            rows.append(d)
    df=pd.DataFrame(rows)
    cols=["method","regime","SNPs","indels","multiallelic sites"]
    return df[[c for c in cols if c in df.columns]].to_html(index=False, classes="restable", border=0)

def load_summary():
    p=PLOTS/"summary.json"
    return json.load(open(p)) if p.exists() else {}

H=[]
def add(s): H.append(s)

add("""<!DOCTYPE html><html><head><meta charset="UTF-8">
<title>VCF comparison results v2 — Pnapi chr FR997704</title>
<style>
 body{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;max-width:1200px;margin:2em auto;padding:0 2em;color:#222;line-height:1.55;}
 h1{border-bottom:3px solid #333;padding-bottom:0.3em;}
 h2{color:#1b4f72;border-bottom:1px solid #ddd;padding-bottom:0.2em;margin-top:2.2em;}
 h3{color:#117864;margin-top:1.4em;}
 .restable{border-collapse:collapse;margin:1em 0;font-size:0.9em;}
 .restable td,.restable th{border:1px solid #ccc;padding:0.4em 0.7em;text-align:right;}
 .restable th{background:#1b4f72;color:white;text-align:center;}
 .restable tr:nth-child(even){background:#f8f8f8;}
 .key{background:#fef9e7;border-left:4px solid #b7950b;padding:0.7em 1em;margin:1em 0;}
 .info{background:#f4f9fd;border-left:4px solid #1b4f72;padding:0.7em 1em;margin:1em 0;}
 .warn{background:#fdedec;border-left:4px solid #c0392b;padding:0.7em 1em;margin:1em 0;}
 img{display:block;margin:0.6em auto;}
 code{background:#eee;color:#c0392b;padding:0.1em 0.3em;border-radius:3px;}
 .grid2{display:grid;grid-template-columns:1fr 1fr;gap:1em;align-items:start;}
 dl dt{font-weight:600;color:#1b4f72;margin-top:0.5em;}
</style></head><body>""")

add("<h1>VCF comparison — results (v3)</h1>")
add("<p><em>Pieris napi</em>, chromosome FR997704.1 (= RefSeq NC_062243.1, chr 10). "
    "Six VCFs from the same reads, 14 samples in two populations of 7 (A vs P).</p>")

# ---------------- TABLE OF CONTENTS -----------------------------------------
add("""
<nav style="background:#f4f9fd;border:1px solid #cfe0ef;border-radius:8px;padding:1em 1.4em;margin:1.5em 0;">
<div style="font-weight:700;color:#1b4f72;font-size:1.05em;margin-bottom:0.4em;">Contents</div>
<ol style="margin:0;padding-left:1.4em;line-height:1.9;">
  <li><a href="#t0">What is being compared, and why</a></li>
  <li><a href="#t1">Per-file QC (bcftools stats)</a></li>
  <li><a href="#t2">Site-set comparison (bcftools isec)</a>
      <ul style="list-style:none;padding-left:0.6em;margin:0.1em 0;color:#555;">
        <li>↳ <a href="#t2a">2a. Simple site overlap (raw)</a></li>
        <li>↳ <a href="#t2b">2b. Full multi-way intersection (UpSetR, 6- &amp; 12-way)</a></li>
        <li>↳ <a href="#t2c">2c. Pairwise Jaccard similarity (incl. all-12 matrix)</a></li>
      </ul></li>
  <li><a href="#t3">Genomic context of SNPs (bcftools csq)</a></li>
  <li><a href="#t4">Sliding-window population-genomic estimators (pixy)</a>
      <ul style="list-style:none;padding-left:0.6em;margin:0.1em 0;color:#555;">
        <li>↳ <a href="#t4a">4a. π — nucleotide diversity</a></li>
        <li>↳ <a href="#t4b">4b. Fₛₜ — population differentiation</a></li>
        <li>↳ <a href="#t4c">4c. dₓy — between-population divergence</a></li>
        <li>↳ <a href="#t4d">4d. One method vs the other five (π)</a></li>
        <li>↳ <a href="#t4e">4e. Normalised vs raw π, per method</a></li>
        <li>↳ <a href="#t4f">4f. All measures together (combined genome scan, with depth)</a></li>
        <li>↳ <a href="#t4g">4g. Fₛₜ estimator comparison — Weir &amp; Cockerham vs Hudson</a></li>
        <li>↳ <a href="#t4h">4h. Mean read depth per window — all six files</a></li>
      </ul></li>
  <li><a href="#t5">What normalisation actually did (raw vs normalised)</a></li>
  <li><a href="#t6">Integrative analysis — what kind of SNPs do methods disagree on?</a></li>
  <li><a href="#t7">Bottom line</a></li>
  <li><a href="#t8">Still on the table (not yet run)</a></li>
</ol>
</nav>
""")

# ---------------- INTRO -----------------------------------------------------
add('<h2 id="t0">0. What is being compared, and why</h2>')
add("""<p>All six VCFs were produced from the <strong>same sequencing reads</strong> and the
<strong>same variant caller</strong> (bcftools). They differ along two axes:</p>""")
add("""<dl>
<dt>Mapper (2 levels)</dt>
<dd><b>bwa</b> = BWA-MEM; <b>ngm</b> = NextGenMap. Different aligners place ambiguous / repetitive
reads differently, which changes which positions become variant.</dd>
<dt>Post-call filtering (3 levels)</dt>
<dd><b>cohort</b> = baseline cohort genotyping; <b>hwe</b> = sites passing a Hardy–Weinberg-equilibrium
filter removed/flagged; <b>nohwe</b> = the HWE filter not applied. HWE filtering removes sites whose
genotype proportions depart from HWE — which can be technical artefacts <em>or</em> real biology
(selection, population structure, sex-linkage).</dd>
</dl>
<p>So the 2 × 3 grid is: <code>{bwa,ngm} × {cohort,hwe,nohwe}</code>. Because reads and caller are
held constant, every difference between VCFs is attributable to (a) read placement (mapper) or
(b) the filtering decision — exactly the two things this assessment isolates.</p>""")

add("""<div class="info"><b>What normalisation is expected to do.</b> The same biological variant
can be written in several valid ways in a VCF (indels placed at different positions in a repeat,
multi-allelic sites packed into one row vs split across rows, MNPs vs adjacent SNPs, untrimmed
REF/ALT). <code>bcftools norm -f REF -m -any</code> rewrites every record into one canonical form:
<b>left-align</b> indels, <b>trim</b> REF/ALT to minimal representation, and <b>split</b> multi-allelic
sites into separate bi-allelic rows. The <em>expected</em> effect is that two VCFs which are
biologically identical but textually different will look identical after normalisation — so any
remaining disagreement is real, not cosmetic. Section 5 shows what it actually did here (including a
side-effect on π).</div>""")

# pipeline + VEP note
add("""<div class="warn"><b>Annotation tool note.</b> Ensembl VEP from bioconda fails to load on this
Mac (a dyld symbol error in its perl bundle — a known macOS issue). Tier 3 uses <code>bcftools csq</code>
instead: same GFF + FASTA inputs, same consequence taxonomy. <b>v2 fix:</b> the v1 annotation piped csq
through <code>+split-vep</code>, which silently dropped every variant with no consequence tag — i.e. all
intergenic SNPs — making it look like ~98% of SNPs were intronic. v2 keeps those variants and labels
them <em>intergenic / no annotated feature</em>.</div>""")

# ---------------- TIER 1 ----------------------------------------------------
add('<h2 id="t1">1. Per-file QC (bcftools stats)</h2>')
add("<p>SNP / indel / multi-allelic counts per VCF, raw and normalised.</p>")
add(stats_table())
add("""<div class="key"><b>Read this:</b> within each mapper, the three filters give very similar counts,
with <code>nohwe</code> largest (no sites removed). Between mappers, <b>BWA calls ~10% more SNPs than NGM</b>.
After normalisation the <em>multi-allelic sites</em> column drops to 0 (all decomposed) and SNP counts
rise — see Section 5.</div>""")

# ---------------- TIER 2 ----------------------------------------------------
add('<h2 id="t2">2. Site-set comparison (bcftools isec)</h2>')
add("<p>How much of the apparent disagreement between methods is real, and how much is just representation?</p>")
add(isec_summary())
add("""<div class="key"><b>Headline:</b> the key number is the <em>fraction of SNP sites shared by all six
methods</em>, which <b>rises from 68.8% (raw) to 71.9% (normalised)</b> — normalisation improves
cross-method concordance. The union <em>count</em> itself goes slightly up (1.71 M → 1.79 M), not down,
because <code>-m -any</code> splits multi-allelic SNP sites into more bi-allelic records; the right way to
read the effect is therefore the shared-fraction column, not the raw union size.</div>""")

add('<h3 id="t2a">2a. Simple site overlap (raw)</h3>')
add("<p>The compact view kept from v1 — a quick read on how the six sets pile up.</p>")
add(img("isec_sites_raw.png","site overlap raw"))

add('<h3 id="t2b">2b. Full multi-way intersection — UpSetR</h3>')
add("""<p>Produced in R with the <b>UpSetR</b> package (via <code>scripts/upset_tier2.Rmd</code>). Each
bar is the number of SNP sites belonging to <em>exactly</em> the method-combination marked by the filled
dots below it; the left-hand bars are each method's total set size (blue = BWA, purple = NGM). UpSet is
the readable alternative to a 6-way Venn diagram (which would need 2⁶−1 = 63 regions).</p>""")
add("<h4>Raw</h4>")
add(img("upsetR_raw.png","UpSetR raw"))
add("<h4>Normalised</h4>")
add(img("upsetR_norm.png","UpSetR norm"))
add("""<div class="key"><b>What the UpSet plot reveals:</b> the largest bar is the core shared by all six
(~1.18 M). The next two largest intersections are <em>mapper-specific</em>: ~200 k sites shared by all
three BWA settings but absent from NGM, and ~96 k shared by all three NGM settings but absent from BWA.
The filter axis (cohort/hwe/nohwe) produces only small private sets. <b>Conclusion: the mapper is the
dominant source of method-specific SNPs.</b></div>""")

add("<h3>2b-bis. All 12 call sets at once — raw vs normalised in one plot</h3>")
add("""<p>The same UpSet idea, now with <b>all 12 sets</b> (each of the six methods in both its raw and
normalised form). Set-size bars are coloured by regime: <span style="color:#7a8a99;">grey = raw</span>,
<span style="color:#1b4f72;">dark blue = normalised</span>. This lets you read the effect of normalisation
directly: for any given method, the raw and normalised sets should land in the same intersections for
ordinary bi-allelic SNPs, with normalisation adding rows where it split multi-allelic sites.</p>""")
add(img("upsetR_all12.png","UpSetR 12-way"))
add("""<div class="key"><b>How to read it:</b> the dominant bar is the core shared by all 12 (the same
SNP called by every method in both regimes). Intersections that include all 6 of one regime but not the
other isolate variants that exist <em>only after</em> (or <em>only before</em>) normalisation — i.e. the
representational differences. Because raw and normalised of the same method otherwise co-occur, the plot
confirms that normalisation reshapes representation rather than discovering or losing real variants.</div>""")

add('<h3 id="t2c">2c. Pairwise Jaccard similarity — what it is and how to read it</h3>')
add("""<p>The UpSet plot shows every intersection at once; the <b>Jaccard matrix</b> condenses the
pairwise picture into a single similarity score for every pair of methods.</p>
<p>For two SNP-site sets <i>A</i> and <i>B</i>, the <b>Jaccard index</b> is</p>
<p style="text-align:center;font-size:1.1em;">J(A,B) = |A ∩ B| / |A ∪ B|
&nbsp;=&nbsp; (sites called by <em>both</em>) / (sites called by <em>either</em>)</p>
<p>It ranges from 0 (no shared sites) to 1 (identical site sets). Unlike a raw count of shared sites,
Jaccard is normalised by set size, so a method that simply calls more variants does not automatically
look &ldquo;more similar&rdquo; to everything. In the heatmaps below, each cell is J for that pair;
the diagonal is 1 by definition. We restrict to SNP sites (variant positions), so this measures
<em>site concordance</em> — agreement on <em>where</em> a SNP is — not genotype concordance.</p>
<p><b>How to read it:</b> bright blocks along the BWA–BWA and NGM–NGM sub-squares mean within-mapper
agreement is high; dimmer off-diagonal blocks (BWA vs NGM) mean cross-mapper agreement is lower. The
gap between the two is the mapper effect, quantified.</p>""")
add('<div class="grid2">')
add("<div><h4>Raw</h4>"+img("isec_jaccard_raw.png","jaccard raw")+"</div>")
add("<div><h4>Normalised</h4>"+img("isec_jaccard_norm.png","jaccard norm")+"</div>")
add("</div>")
add("""<div class="key"><b>Numbers:</b> within-mapper Jaccard is 0.92–0.98; between-mapper is 0.72–0.78.
Normalisation nudges every value up by ~0.01–0.02 but does <em>not</em> close the BWA-vs-NGM gap — that
gap is real biology/mapping, not representation.</div>""")

add("<h4>All 12 VCFs in one Jaccard matrix</h4>")
add("""<p>The same pairwise Jaccard, now across all 12 call sets at once (raw block top-left, normalised
block bottom-right, separated by the white lines). This puts every raw-vs-normalised and cross-mapper
comparison on a single colour scale.</p>""")
add(img("jaccard_all12.png","12-way jaccard", maxw="900px"))
add("""<div class="key"><b>Three things to read off it:</b> (1) the same method's raw and normalised sets
are highly similar (the <code>x_raw</code> vs <code>x_norm</code> cells are ~0.93–0.95) — normalisation
keeps the site set nearly intact; (2) the four 6×6 mapper blocks (BWA-raw, NGM-raw, BWA-norm, NGM-norm)
repeat the within-mapper-high / between-mapper-low pattern regardless of regime; (3) crucially, a BWA set
is no more similar to an NGM set <em>after</em> normalisation than before — confirming once more that the
BWA-vs-NGM gap is biological/mapping, not a representation artefact that normalisation could fix.</div>""")

# ---------------- TIER 3 ----------------------------------------------------
add('<h2 id="t3">3. Genomic context of SNPs (bcftools csq + extracted GFF)</h2>')
add("""<p>Every SNP classified by the GFF feature it falls in. <b>Important distinction:</b>
<em>intronic</em> = the SNP is inside an annotated gene but in an intron; <em>intergenic / no annotated
feature</em> = the SNP is outside every annotated transcript (bcftools csq emits no consequence tag for
these). The two are <b>separate classes</b> — intronic does <em>not</em> include intergenic. (The v1
plot accidentally dropped the intergenic class entirely; that is fixed here.)</p>""")
add("<h3>3a. All classes (including intergenic)</h3>")
add(img("vep_consequence_by_method.png","csq all classes"))
add("""<div class="key">Across all methods, roughly <b>~39% of SNPs are intergenic</b> and
<b>~60% intronic</b>, leaving ~1% in coding/splice classes. The proportions are nearly identical across
the six methods — i.e. the methods agree on <em>where in the genome</em> SNPs fall, even where they
disagree on the exact set.</div>""")
add("<h3>3b. Coding &amp; splice classes only (intron + intergenic removed)</h3>")
add("<p>Zooming in on the functional ~1% so the rare, biologically important classes are visible.</p>")
add(img("vep_consequence_coding_only.png","csq coding only"))
add("""<div class="key">Of coding/splice SNPs: ~55% synonymous, ~32% splice_region, ~13% missense, with
small numbers of splice donor/acceptor, stop/start gained-lost, and non-coding-RNA variants. Again,
strikingly consistent across methods — annotation conclusions are robust to mapper/filter choice.</div>""")

# ---------------- TIER 4 : SLIDING-WINDOW POPGEN ---------------------------
add('<h2 id="t4">4. Sliding-window population-genomic estimators (pixy, 50 kb windows)</h2>')
add("""<p>All three windowed estimators are presented together below, in the order <b>nucleotide
diversity (π) → Fₛₜ → dₓy</b>. For each, per-method tracks along the chromosome (raw and normalised) are
followed by the cross-method Spearman correlation matrix. Two π-specific comparison diagnostics
(method-vs-method, and raw-vs-normalised per method) close the section.</p>""")

# --- 4a. pi ---
add('<h3 id="t4a">4a. π — nucleotide diversity (within-population)</h3>')
add('<div class="grid2">')
add("<div><h4>raw</h4>"+img("pi_per_method_raw.png")+"</div>")
add("<div><h4>normalised</h4>"+img("pi_per_method_norm.png")+"</div>")
add("</div>")
add("<h4>Cross-method Spearman correlation of windowed π</h4>")
add('<div class="grid2">')
add("<div><h4>raw</h4>"+img("pi_correlation_raw.png")+"</div>")
add("<div><h4>normalised</h4>"+img("pi_correlation_norm.png")+"</div>")
add("</div>")
add("""<div class="key">π correlates 0.96–0.999 across all method pairs — the <em>shape</em> of the
diversity landscape is highly robust to method choice. (But absolute π is sensitive to normalisation —
see Section 5.)</div>""")

# --- 4b. Fst ---
add('<h3 id="t4b">4b. Fₛₜ — population differentiation (A vs P)</h3>')
add("<h4>Windowed Fₛₜ tracks — Weir &amp; Cockerham (default estimator)</h4>")
add('<div class="grid2">')
add("<div><h4>raw</h4>"+img("fst_per_method_raw.png")+"</div>")
add("<div><h4>normalised</h4>"+img("fst_per_method_norm.png")+"</div>")
add("</div>")
add("""<h4>Windowed Fₛₜ tracks — Hudson estimator (<code>--fst_type hudson</code>)</h4>""")
add("""<p>The same windowed Fₛₜ computed with the Hudson estimator, shown here alongside the default so the
two are directly comparable along the chromosome. The peaks fall at the same loci; Hudson runs a little
higher at them. The head-to-head per-window comparison and the cross-method robustness of Hudson are in
§4g.</p>""")
add('<div class="grid2">')
add("<div><h4>raw</h4>"+img("fst_hudson_per_method_raw.png","hudson tracks raw")+"</div>")
add("<div><h4>normalised</h4>"+img("fst_hudson_per_method_norm.png","hudson tracks norm")+"</div>")
add("</div>")
add("<h4>Cross-method Spearman correlation of windowed Fₛₜ (Weir &amp; Cockerham)</h4>")
add('<div class="grid2">')
add("<div><h4>raw</h4>"+img("fst_correlation_raw.png")+"</div>")
add("<div><h4>normalised</h4>"+img("fst_correlation_norm.png")+"</div>")
add("</div>")
add("""<div class="warn"><b>Fₛₜ is NOT robust.</b> Cross-method correlation ranges 0.34–0.88, and
within-mapper agreement is sometimes worse than between-mapper. Fₛₜ depends on rare alleles, which is
exactly where methods diverge.</div>""")
add("<h4>Top-1% Fₛₜ outlier-window agreement (Jaccard)</h4>")
add('<div class="grid2">')
add("<div><h4>raw</h4>"+tsv_table("fst_outlier_jaccard_raw.tsv")+"</div>")
add("<div><h4>normalised</h4>"+tsv_table("fst_outlier_jaccard_norm.tsv")+"</div>")
add("</div>")
add("""<div class="warn"><b>Most important caveat for downstream work:</b> the top-1% Fₛₜ windows
(candidate regions under selection) overlap only 0–50% between method pairs. The list of selection
candidates is heavily method-dependent — always report it across mappers, not from a single VCF.</div>""")
add("""<p><em>This Fₛₜ uses the default Weir &amp; Cockerham estimator; §4g compares it against the Hudson
estimator (<code>--fst_type hudson</code>) — the choice changes peak values and interacts with the calling
method.</em></p>""")

# --- 4c. dxy ---
add('<h3 id="t4c">4c. dₓy — between-population divergence (A vs P)</h3>')
add("""<p>dₓy is the average number of pairwise differences <em>between</em> the two populations per site —
the between-population analogue of π. Like π it is an absolute diversity measure (a raw difference count),
not a ratio, so we expect it to behave like π rather than like Fₛₜ.</p>""")
add('<div class="grid2">')
add("<div><h4>raw</h4>"+img("dxy_per_method_raw.png")+"</div>")
add("<div><h4>normalised</h4>"+img("dxy_per_method_norm.png")+"</div>")
add("</div>")
add("<h4>Cross-method Spearman correlation of windowed dₓy</h4>")
add('<div class="grid2">')
add("<div><h4>raw</h4>"+img("dxy_correlation_raw.png")+"</div>")
add("<div><h4>normalised</h4>"+img("dxy_correlation_norm.png")+"</div>")
add("</div>")
add("""<div class="key"><b>dₓy is robust, like π — and unlike Fₛₜ.</b> Cross-method correlation is
0.97–0.999 (within-mapper ≥0.995, between-mapper ~0.97–0.98), so the divergence landscape is essentially
method-independent in shape. This confirms the pattern: <em>absolute</em> diversity measures (π, dₓy) are
method-robust, whereas the <em>ratio</em> statistic Fₛₜ (§4b) is fragile because it amplifies the
denominator differences where methods disagree. Note dₓy is subject to the <em>same</em> ~11.5%
normalisation inflation as π — see §5c.</div>""")

# --- 4d. one-vs-five pi diagnostic ---
add('<h3 id="t4d">4d. One method vs the other five (windowed π, head-to-head)</h3>')
add("""<p>The same diagnostic as the raw-vs-norm scatter in §5b, but used to compare <em>methods</em>:
<code>bwa_cohort</code> is fixed on the x-axis and each of the other five methods is plotted on the y-axis,
within a single regime. Each panel reports the slope of a through-origin fit (1.0 = identical magnitude)
and the Spearman ρ (1.0 = identical ranking). This is the &ldquo;explore one reference against everything&rdquo;
view; the reference method is a single parameter in <code>scripts/integrative_analysis.py</code> if you
want to pivot on a different one.</p>""")
add("<h4>Within raw</h4>")
add(img("pi_method_compare_raw.png","pi one-vs-five raw"))
add("<h4>Within normalised</h4>")
add(img("pi_method_compare_norm.png","pi one-vs-five norm"))
add("""<div class="key"><b>What it shows:</b> against the two other BWA settings, slope ≈ 1.00 and ρ ≈ 0.999
— same mapper, so π is essentially identical regardless of the filter. Against the three NGM settings the
slope drops to ~0.91–0.96 and the scatter widens (ρ ≈ 0.96–0.97): NGM yields systematically slightly
<em>lower</em> windowed π than BWA. So the residual π difference between methods is, again, a mapper effect
— not a filtering effect.</div>""")

# --- 4e. per-method norm vs raw pi diagnostic ---
add('<h3 id="t4e">4e. Normalised vs raw π, one panel per method</h3>')
add("""<p>The §5b raw-vs-normalised π scatter, now broken out into one panel per method so the
normalisation effect can be inspected method-by-method. In every panel raw π is on the x-axis and
normalised π on the y-axis; each title gives the through-origin slope, Spearman ρ, and the mean per-window
π inflation.</p>""")
add(img("pi_norm_vs_raw_per_method.png","per-method raw vs norm pi"))
add("""<div class="key"><b>The effect is remarkably uniform across methods:</b> every panel shows slope
≈ 1.12, ρ ≈ 0.998, and a <b>+11–12% mean π inflation</b> with 97–99% of windows higher after
normalisation. Neither the mapper nor the filter changes the size of the effect — it is driven entirely by
<code>norm -m -any</code> splitting multi-allelic sites, which pixy then counts per allele. Reinforces the
practical rule in §5b: use the un-split all-sites VCFs for π/dₓy/Fₛₜ.</div>""")

# --- 4f. all three measures together (combined genome scan) ---
add('<h3 id="t4f">4f. All measures together — combined genome scan (with depth)</h3>')
add("""<p>The windowed estimators (§4a–4c, §4g) plus <b>mean read depth</b> (§4h), stacked in one figure and
<b>sharing the x-axis</b>, so they can be read against each other along the chromosome: <b>depth</b> on top
(the QC context), then <b>π</b> (mean of the two populations), <b>Fₛₜ (Weir &amp; Cockerham)</b>,
<b>Fₛₜ (Hudson)</b>, and <b>dₓy</b> (all A vs P). Both Fₛₜ estimators are shown so they can be compared
against each other and against π/dₓy; depth on top lets you check whether a signal coincides with a coverage
anomaly. All six methods are overlaid in every panel (depth is mapper-driven, so the three lines per mapper
coincide).</p>""")
add('<div class="grid2">')
add("<div><h4>raw</h4>"+img("combined_tracks_raw.png","combined tracks raw")+"</div>")
add("<div><h4>normalised</h4>"+img("combined_tracks_norm.png","combined tracks norm")+"</div>")
add("</div>")
add("""<div class="key"><b>What the aligned view shows:</b> (1) π and dₓy track each other almost exactly —
both are absolute diversity measures and both bundle tightly across the six methods; (2) the regions where
π and dₓy crash toward zero (≈3 Mb and ≈7 Mb) coincide with the sharp <b>Fₛₜ peaks</b> in <em>both</em>
estimator panels — the classic signature of low-diversity / reduced-recombination regions where
differentiation is inflated; (3) the two Fₛₜ panels agree on <em>where</em> the peaks are but Hudson reaches
<em>higher</em> at them (see §4g); (4) the method lines are visibly more spread in the Fₛₜ panels than in the
π/dₓy panels — the same robustness contrast quantified in §4a–4b, now visible at a glance;
<b>(5) the depth panel is the key QC overlay</b>: the ≈7 Mb Fₛₜ peak sits directly under a sharp
<b>depth spike</b> (≈22× vs the ≈9× background), the hallmark of a collapsed repeat / paralog where reads
from duplicated copies pile up and mis-call differentiation — i.e. that peak is most likely a mapping
artefact, not selection. This panel is the synthesis of section 4: it shows the biology (covarying diversity
and differentiation), the methodology (which statistics/estimators are method-stable), and the QC (which
signals coincide with coverage anomalies).</div>""")

# --- 4g. Fst estimator comparison: WC vs Hudson ---
add('<h3 id="t4g">4g. Fₛₜ estimator comparison — Weir &amp; Cockerham vs Hudson</h3>')
add("""<p>The Fₛₜ in §4b uses pixy's default <b>Weir &amp; Cockerham (1984)</b> estimator. pixy also offers
the <b>Hudson (1992; Bhatia et al. 2013)</b> estimator (<code>--fst_type hudson</code>), which uses a
ratio-of-averages formulation that is less biased under unequal sample sizes and is often preferred for
genome scans. Both were run on all six methods (raw and normalised); below, each window's WC value (x) is
plotted against its Hudson value (y), one panel per method.</p>""")
add("<h4>Within raw</h4>")
add(img("fst_wc_vs_hudson_raw.png","wc vs hudson raw"))
add("<h4>Within normalised</h4>")
add(img("fst_wc_vs_hudson_norm.png","wc vs hudson norm"))
add(tsv_table("fst_estimator_summary.tsv"))
add("""<div class="key"><b>The estimator choice is not cosmetic — and it interacts with the calling method.</b>
(1) Hudson is <em>systematically slightly higher</em> than WC at differentiated windows (mean Hudson−WC
+0.01 to +0.045; points sit above y=x), so reported peak-Fₛₜ values depend on which estimator you cite.
(2) Agreement is <b>method-dependent</b>: for the <code>cohort</code>/<code>hwe</code> call sets the two
track each other tightly (Pearson 0.85–0.99), but for <code>bwa_nohwe</code>, <code>ngm_cohort</code> and
<code>ngm_nohwe</code> they diverge markedly (Pearson 0.12–0.55; rank ρ down to ~0.53). The divergence is
concentrated at the high-Fₛₜ outlier windows — exactly the windows a selection scan cares about. (Pearson
is pulled down by those few high-leverage windows; Spearman is the more stable summary.)</div>""")
add("<h4>Does Hudson change the cross-method robustness picture?</h4>")
add('<div class="grid2">')
add("<div><h4>raw</h4>"+img("fst_hudson_correlation_raw.png","hudson corr raw")+"</div>")
add("<div><h4>normalised</h4>"+img("fst_hudson_correlation_norm.png","hudson corr norm")+"</div>")
add("</div>")
add("""<div class="key"><b>Hudson is somewhat more method-robust than WC.</b> Cross-method Spearman
correlation of windowed Hudson Fₛₜ ranges ~0.58–0.94 (raw), a higher floor than the WC equivalent in §4b
(0.34–0.88). So while Fₛₜ remains the most method-sensitive statistic overall, the Hudson estimator
tightens cross-method agreement — a reason to prefer it for between-method comparisons. <b>Bottom line:</b>
report which Fₛₜ estimator you used, and ideally show both; the WC-vs-Hudson gap is largest precisely at the
outlier windows that drive biological interpretation.</div>""")

# --- 4h. read depth per window ---
add('<h3 id="t4h">4h. Mean read depth per window — all six files</h3>')
add("""<p>Average per-sample read depth in the same 50 kb windows, from
<code>vcftools --gzvcf &lt;raw all-sites VCF&gt; --site-mean-depth</code> (which averages the per-sample
<code>FORMAT/DP</code> at every position across the 14 samples), binned to windows and computed
<b>independently for every one of the six VCFs</b> (no assumption that they share a depth profile). Because
the input is the all-sites VCF, depth is sampled at essentially every base, not just at SNPs.</p>""")
add("<h4>One panel per VCF</h4>")
add(img("depth_per_file.png","depth per file"))
add("<h4>Genome-wide mean depth per file, and all files overlaid</h4>")
add('<div class="grid2">')
add("<div>"+img("depth_mean_per_file.png","mean depth per file")+"</div>")
add("<div>"+img("depth_per_method.png","depth overlaid")+"</div>")
add("""<div class="key"><b>Empirically (md5-verified, not assumed):</b> the three filters within a mapper
produce <em>byte-identical</em> depth files — depth is set by the mapping, and cohort/hwe/nohwe do not change
it. So the six files collapse to two distinct profiles. <b>BWA 8.94× vs NGM 8.16×</b> (~9% more): different
data could of course differ by filter, which is exactly why each file is computed and shown separately
here.</div>""")
add("</div>")
add("""<div class="key"><b>Two things depth explains:</b> (1) the <b>BWA 8.94× vs NGM 8.16×</b> gap is the
direct cause of BWA calling ~10% more SNPs (§1): more reads mapped → more callable sites. (2) The
<b>≈7 Mb depth spike to ~22×</b> (≈2.5× the background, visible in every BWA and NGM panel) flags a
<b>collapsed-repeat / paralog region</b>: reads from duplicated copies pile up there, inflating apparent
heterozygosity differences. That window is exactly where Fₛₜ peaks (§4b) and π/dₓy collapse — so the depth
track reclassifies that &ldquo;differentiation peak&rdquo; as a likely <b>mapping artefact</b>, not a
selection signal. <b>This is why depth belongs in the scan</b>: it is the first thing to check before
believing an Fₛₜ outlier. A practical follow-up is to mask windows with depth &gt; ~2× the chromosome median
(and unusually low-depth windows) before any selection scan.</div>""")

add("<h4>Per-individual depth traces — bwa_hwe vs bwa_nohwe, overlaid</h4>")
add("""<p>The plots above average depth across the 14 samples. Here instead is the depth of <b>each
individual</b> (its own colour) along the chromosome, for two VCFs overlaid in one panel — <b>bwa_hwe</b>
(solid) and <b>bwa_nohwe</b> (dashed). Per-sample depth comes from <code>bcftools query -f
'%CHROM\\t%POS[\\t%DP]\\n'</code> binned to 50 kb windows.</p>""")
add(img("depth_per_individual.png","per-individual depth traces", maxw="100%"))
add("""<div class="key"><b>What this adds:</b> (1) the solid (hwe) and dashed (nohwe) lines lie on top of
each other for every individual — confirming, at the per-sample level, that the HWE filter does not change
coverage; (2) individuals vary in background coverage (the band spans ~5–15×), useful for spotting a
low-coverage sample that might drive missingness; (3) the ≈7 Mb spike is not uniform — it is driven by a
subset of individuals reaching ~40×, exactly the per-sample signature of reads from a duplicated/collapsed
region piling up in some samples, reinforcing that the Fₛₜ peak there (§4b) is a mapping artefact rather
than biology. (Generated by <code>05_popgen/depth_per_individual.py</code>; change the two VCFs at the top
to compare any pair.)</div>""")

# ---------------- TIER 5: NORM vs RAW --------------------------------------
add('<h2 id="t5">5. What normalisation actually did (raw vs normalised, head-to-head)</h2>')
add("<h3>5a. Record counts before vs after <code>bcftools norm -m -any</code></h3>")
add(img("norm_vs_raw_counts.png","norm vs raw counts"))
add("""<div class="key"><b>This is normalisation in action:</b> the <em>multi-allelic sites</em> bars
collapse to zero after normalisation (every multi-allelic row is decomposed into bi-allelic rows), while
<em>SNP</em> and <em>indel</em> record counts rise correspondingly (the split rows + left-aligned indels).
No variants are lost — they are re-expressed.</div>""")
add("<h3>5b. Side-effect: normalisation inflates π</h3>")
add(img("norm_vs_raw_pi.png","norm vs raw pi", maxw="640px"))
sm = load_summary().get("norm_vs_raw_pi",{})
if sm:
    add(f"""<div class="warn"><b>Watch out.</b> Splitting multi-allelic sites with <code>-m -any</code> is
    correct for site-set comparison and annotation, but it <b>systematically inflates π</b> when the
    split VCF is then fed to pixy: {sm.get('pi_pct_windows_norm_higher',0):.0f}% of windows have higher π
    after normalisation, by <b>{sm.get('pi_mean_relative_pct',0):.1f}% on average</b>. Mechanism: pixy
    counts the pairwise differences of each split bi-allelic record separately, so a multi-allelic
    position contributes its differences more than once (verified on a single window: count_diffs +6.6%,
    count_comparisons +0.5%). <b>Practical rule:</b> for pixy / diversity estimation, do <em>not</em>
    pre-split multi-allelics — keep them merged (or use a pixy-aware all-sites pipeline). Use the split,
    normalised VCFs for isec and annotation, but the un-split all-sites VCFs for π/dxy/Fst.</div>""")

add("<h3>5c. The same inflation hits dₓy</h3>")
add(img("norm_vs_raw_dxy.png","norm vs raw dxy", maxw="640px"))
if sm and sm.get("dxy_mean_relative_pct") is not None:
    add(f"""<div class="warn">dₓy is affected by the identical mechanism:
    <b>{sm.get('dxy_pct_windows_norm_higher',0):.0f}% of windows</b> have higher dₓy after normalisation,
    by <b>{sm.get('dxy_mean_relative_pct',0):.1f}% on average</b> — essentially the same magnitude as the π
    inflation. Both absolute diversity statistics (π and dₓy) are biased upward by multi-allelic splitting,
    so the &ldquo;use un-split all-sites VCFs&rdquo; rule applies to dₓy as well. (Fₛₜ, being a ratio of
    such quantities, largely cancels this out — its raw-vs-norm correlation stays high — but it is fragile
    for a different reason, §4b.)</div>""")

# ---------------- TIER 7: INTEGRATIVE --------------------------------------
add('<h2 id="t6">6. Integrative analysis — what kind of SNPs do the methods disagree on?</h2>')
add("""<p>The previous tiers established <em>that</em> methods disagree and <em>where</em> the disagreement
sits along the mapper/filter axes. This section asks the biologically pointed question: of the SNPs that
are method-discordant, are they a random sample of the genome, or are they enriched for particular kinds
of variant? We classify every normalised union SNP by its sharing pattern and cross it with (a) the
consequence class from Tier 3 and (b) the transition/transversion ratio.</p>""")
cnt = PLOTS/"sharing_category_counts.tsv"
if cnt.exists():
    _c = pd.read_csv(cnt, sep="\t"); _c.columns = ["sharing category","n SNP sites"]
    add(_c.to_html(index=False, classes="restable", border=0))
add("<h3>6a. Consequence-class enrichment of discordant SNPs</h3>")
add("""<p>For each sharing category we take the proportion of SNPs in each consequence class and compare it
to the shared-by-all-6 core, as a log2 fold-change (red = over-represented relative to the core, blue =
under). The core row is 0 by construction.</p>""")
add(img("enrichment_consequence.png","consequence enrichment"))
add("""<div class="key"><b>The clear, large-N result:</b> method-discordant SNPs are strongly
<b>depleted in synonymous (&minus;2.6 to &minus;3.1 log2FC) and missense (&minus;1.3 to &minus;2.3)</b>
variants. Ordinary coding SNPs in well-mapped exons are reproducibly called by every method — they are the
<em>most</em> trustworthy. <b>Conversely, discordant SNPs are enriched in apparent high-impact classes</b>:
splice-donor/acceptor (+0.8 to +2.0), stop/start gained-lost, and non-coding-RNA variants.
<b>Caution:</b> the very variants one gets excited about — a gained stop codon, a disrupted splice site —
are disproportionately the method-unstable ones, so confirm them across mappers before believing them.
(The stop/start and splice-donor classes have small absolute counts (tens–hundreds), so those particular
log2FC values are noisy; the synonymous/missense depletion is on hundreds of thousands of sites and is
rock-solid.)</div>""")
add("<h3>6b. Ts/Tv ratio by sharing category — a quality signal</h3>")
add("""<p>The transition/transversion ratio is a classic SNP-quality diagnostic: genuine SNPs sit above
the random expectation of 0.5 (transitions are mutationally favoured), whereas false positives sit closer
to random because sequencing/mapping errors have no transition bias. Computed directly from REF/ALT in
each sharing category.</p>""")
add(img("tstv_by_sharing.png","tstv by sharing", maxw="760px"))
add(tsv_table("tstv_by_sharing.tsv"))
add("""<div class="key"><b>Monotonic decline confirms the artefact gradient:</b> SNPs shared by all six
methods have the highest Ts/Tv (0.82), all-BWA-private next (0.76), then all-NGM-private (0.64), and
singletons lowest (0.62) — closest to the 0.5 random floor. The more method-restricted a SNP, the more it
looks like noise. (Absolute Ts/Tv is low because this is genome-wide, dominated by intergenic/intronic
sites; the <em>relative</em> ordering is the signal.) Independent corroboration of 6a: discordant calls
are enriched for likely artefacts.</div>""")
add("""<div class="info"><b>Other angles worth adding next</b> (not yet run): per-category <b>minor-allele-
frequency spectra</b> (artefacts cluster at very low MAF); <b>genotype-level concordance</b> with
hap.py/vcfeval stratified by MAF; and intersecting the discordant SNPs with a <b>repeat / low-complexity
mask</b> to test whether the mapper-private SNPs sit in repetitive sequence.</div>""")

# ---------------- CODA ------------------------------------------------------
add("""<h2 id="t7">7. Bottom line</h2>
<ul>
<li><b>Mapper &gt; filter.</b> BWA vs NGM is the dominant axis of disagreement (UpSet + Jaccard); the
cohort/hwe/nohwe filter is a minor perturbation.</li>
<li><b>Absolute diversity (π and dₓy) is robust in shape, not in magnitude.</b> Window-to-window π and dₓy
both correlate &gt;0.96–0.97 across methods, but both shift up ~12% if you split multi-allelics before pixy.
Use un-split all-sites VCFs for π/dₓy/Fst.</li>
<li><b>Fst is the fragile statistic.</b> Cross-method correlation 0.34–0.88, and outlier
(selection-candidate) windows agree only 0–50% across methods — never trust a selection scan from a single
VCF. (It is a ratio, so it cancels the normalisation inflation but amplifies the calls where methods
disagree.)</li>
<li><b>Annotation is robust.</b> The genomic-context breakdown is near-identical across all six methods.</li>
<li><b>Discordant SNPs are lower-quality.</b> Method-private SNPs are depleted in synonymous/missense, enriched
in apparent high-impact classes, and show a falling Ts/Tv toward the random floor — i.e. disagreement is
concentrated in the least trustworthy calls. Treat single-method high-impact variants with suspicion.</li>
<li><b>Depth explains the mapper gap and flags artefacts.</b> BWA averages 8.9× vs NGM 8.2× (→ BWA's extra
SNPs); and the ≈7 Mb Fₛₜ peak sits on a ~22× depth spike (collapsed repeat), so it is most likely a mapping
artefact. Overlay per-window depth (§4h) before trusting any Fₛₜ outlier.</li>
</ul>""")
add("""<h2 id="t8">Still on the table (not yet run)</h2>
<ul>
<li><b>VEP proper</b> (Docker/Linux) and <b>SnpEff</b> — independent cross-checks on the csq annotation.</li>
<li><b>repeats.bed intersection</b> — test whether the ~200k BWA-only / ~96k NGM-only SNPs sit in repeats.</li>
<li><b>callable.bed (mosdepth on BAMs)</b> — restrict all comparisons to regions covered by both mappers;
the single biggest remaining lever on the BWA-vs-NGM gap.</li>
<li><b>PLINK PCA / IBS / ROH</b> and <b>hap.py/vcfeval genotype concordance</b> — relational structure and
genotype-level (not just site-level) agreement.</li>
</ul>
<p style="color:#888;font-size:0.85em;">Generated by <code>scripts/make_results_html_v3.py</code> from <code>compare_out/</code>.
UpSet plots from <code>scripts/upset_tier2.Rmd</code> + <code>scripts/upset_all12.Rmd</code>;
integrative analysis from <code>scripts/integrative_analysis.py</code>; annotation from <code>scripts/tier3_csq.sh</code>.</p>""")

add("</body></html>")
TARGET.write_text("".join(H))
print("Wrote", TARGET, f"({TARGET.stat().st_size/1e6:.2f} MB)")
