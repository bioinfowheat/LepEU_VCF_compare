"""
Tier 5 + plotting: read outputs from Tiers 1-4, build comparison tables and figures.

Outputs (under compare_out/plots/):
  stats_summary.tsv          per-method ts/tv, het, indel ratio, # SNPs
  isec_sites_raw.png         UpSet plot of raw site sets
  isec_sites_norm.png        UpSet plot of normalised site sets
  isec_jaccard_raw.png       pairwise Jaccard heatmap (raw)
  isec_jaccard_norm.png      pairwise Jaccard heatmap (normalised)
  pi_per_method.png          windowed pi tracks, one line per method, per population
  fst_per_method.png         windowed Fst tracks
  pi_correlation_raw.png     Spearman correlation matrix of windowed pi across methods (raw)
  pi_correlation_norm.png    same (normalised)
  fst_correlation.png        Spearman correlation matrix of windowed Fst across methods
  vep_consequence_by_method.png  stacked bar of VEP consequence categories per method
"""
import argparse, os, glob, re, sys, json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from upsetplot import from_indicators, UpSet

METHODS = ["bwa_cohort","bwa_hwe","bwa_nohwe","ngm_cohort","ngm_hwe","ngm_nohwe"]
METHOD_COLORS = {
    "bwa_cohort":"#1f77b4","bwa_hwe":"#4a90c2","bwa_nohwe":"#7eb0d5",
    "ngm_cohort":"#6c3483","ngm_hwe":"#8e44ad","ngm_nohwe":"#b07cc6"}

def method_ls(m):
    """Linestyle for sliding-window tracks: the MIDDLE shade of each 3-step colour
    gradient (the *_hwe method) is drawn dashed so the three lines in a group
    (cohort=dark solid, hwe=medium dashed, nohwe=light solid) are easy to tell apart."""
    return "--" if str(m).endswith("_hwe") else "-"

def shortname(path):
    b = os.path.basename(path)
    m = re.search(r"Pnapi_(\w+)_bcftools_(\w+)_FR997704", b)
    return f"{m.group(1)}_{m.group(2)}" if m else b


# --- Tier 1: parse bcftools stats ----------------------------------------
def parse_bcftools_stats(path):
    """Return dict of summary numbers from a bcftools stats file."""
    out = {}
    with open(path) as f:
        for line in f:
            if line.startswith("SN\t"):
                p = line.strip().split("\t")
                key = p[2].rstrip(":").replace(" ", "_").lower()
                try: out[key] = float(p[3])
                except: out[key] = p[3]
    return out

def collect_stats(stats_dir, label):
    rows = []
    for f in sorted(glob.glob(os.path.join(stats_dir,"*.stats"))):
        d = parse_bcftools_stats(f)
        d["method"] = shortname(f)
        d["regime"] = label
        rows.append(d)
    return pd.DataFrame(rows)


# --- Tier 2: parse isec sites.txt ----------------------------------------
def parse_isec_sites(sites_path):
    """sites.txt columns: CHR POS REF ALT MASK (one bit per input)."""
    cols = ["chr","pos","ref","alt","mask"]
    df = pd.read_csv(sites_path, sep="\t", header=None, names=cols, dtype=str)
    mask = df["mask"].apply(lambda s: [int(c) for c in s.rstrip("\n")])
    M = pd.DataFrame(mask.tolist(), columns=METHODS).astype(bool)
    return df, M

def upset_plot(M, out_png, title):
    data = from_indicators(METHODS, data=M)
    fig = plt.figure(figsize=(12,6))
    try:
        UpSet(data, show_counts=True, sort_by="cardinality",
              min_subset_size=max(50, int(len(M)*0.001))).plot(fig=fig)
    except Exception as e:
        plt.close(fig)
        # Fallback: bar plot of top 20 intersection patterns
        fig,ax = plt.subplots(figsize=(12,5))
        keys = M.apply(lambda r: "".join("1" if v else "0" for v in r), axis=1)
        vc = keys.value_counts().head(20)
        ax.bar(range(len(vc)), vc.values, color="#1b4f72")
        ax.set_xticks(range(len(vc))); ax.set_xticklabels(vc.index, rotation=90, family="monospace")
        ax.set_ylabel("# sites")
        ax.set_title(title + " (fallback bar; columns = method-presence bitstring)")
    plt.suptitle(title)
    plt.savefig(out_png, dpi=140, bbox_inches="tight"); plt.close()

def jaccard_heatmap(M, out_png, title):
    A = M.values.astype(bool)
    n = A.shape[1]
    J = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            inter = (A[:,i]&A[:,j]).sum()
            uni   = (A[:,i]|A[:,j]).sum()
            J[i,j] = inter/uni if uni else np.nan
    fig,ax = plt.subplots(figsize=(7,6))
    sns.heatmap(J, annot=True, fmt=".3f", xticklabels=METHODS, yticklabels=METHODS,
                cmap="viridis", vmin=0.5, vmax=1.0, ax=ax)
    plt.title(title); plt.xticks(rotation=40, ha="right"); plt.tight_layout()
    plt.savefig(out_png, dpi=140); plt.close()
    return pd.DataFrame(J, index=METHODS, columns=METHODS)


# --- Tier 4: pixy outputs -------------------------------------------------
def load_pixy(pixy_dir, stat):
    """Concatenate <method>_<stat>.txt files from one regime (raw or norm)."""
    frames = []
    for f in sorted(glob.glob(os.path.join(pixy_dir,"**",f"*_{stat}.txt"), recursive=True)):
        d = pd.read_csv(f, sep="\t")
        m = re.search(r"Pnapi_(\w+)_bcftools_(\w+)_FR997704", os.path.basename(f))
        if m:
            d["method"] = f"{m.group(1)}_{m.group(2)}"
        else:
            d["method"] = os.path.basename(os.path.dirname(f))
        frames.append(d)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

def plot_pi_tracks(df_pi, out_png, title):
    if df_pi.empty: return
    fig, axes = plt.subplots(df_pi["pop"].nunique(), 1, figsize=(12,5), sharex=True, sharey=True)
    if df_pi["pop"].nunique() == 1: axes = [axes]
    for ax,(pop,sub) in zip(axes, df_pi.groupby("pop")):
        for m,s in sub.groupby("method"):
            s = s.sort_values("window_pos_1")
            ax.plot(s["window_pos_1"]/1e6, s["avg_pi"], label=m,
                    color=METHOD_COLORS.get(m,"k"), ls=method_ls(m), alpha=0.75, lw=1.0)
        ax.set_ylabel(f"π — pop {pop}")
    axes[-1].set_xlabel("Position on FR997704.1 (Mb)")
    axes[0].set_title(title)
    axes[0].legend(loc="upper right", fontsize=8, ncol=2)
    plt.tight_layout(); plt.savefig(out_png, dpi=140); plt.close()

def plot_fst_tracks(df_fst, out_png, title):
    if df_fst.empty: return
    pairs = df_fst[["pop1","pop2"]].drop_duplicates().values.tolist()
    fig, axes = plt.subplots(len(pairs), 1, figsize=(12, 3*len(pairs)), sharex=True)
    if len(pairs)==1: axes=[axes]
    for ax,(p1,p2) in zip(axes, pairs):
        sub = df_fst[(df_fst["pop1"]==p1)&(df_fst["pop2"]==p2)]
        for m,s in sub.groupby("method"):
            s = s.sort_values("window_pos_1")
            ax.plot(s["window_pos_1"]/1e6, s["avg_wc_fst"], label=m,
                    color=METHOD_COLORS.get(m,"k"), ls=method_ls(m), alpha=0.75, lw=1.0)
        ax.set_ylabel(f"Fst {p1} vs {p2}")
    axes[-1].set_xlabel("Position on FR997704.1 (Mb)")
    axes[0].set_title(title)
    axes[0].legend(loc="upper right", fontsize=8, ncol=2)
    plt.tight_layout(); plt.savefig(out_png, dpi=140); plt.close()

def plot_dxy_tracks(df_dxy, out_png, title):
    """dxy = between-population divergence; same (pop1,pop2,window) shape as Fst."""
    if df_dxy.empty: return
    pairs = df_dxy[["pop1","pop2"]].drop_duplicates().values.tolist()
    fig, axes = plt.subplots(len(pairs), 1, figsize=(12, 3*len(pairs)), sharex=True)
    if len(pairs)==1: axes=[axes]
    for ax,(p1,p2) in zip(axes, pairs):
        sub = df_dxy[(df_dxy["pop1"]==p1)&(df_dxy["pop2"]==p2)]
        for m,s in sub.groupby("method"):
            s = s.sort_values("window_pos_1")
            ax.plot(s["window_pos_1"]/1e6, s["avg_dxy"], label=m,
                    color=METHOD_COLORS.get(m,"k"), ls=method_ls(m), alpha=0.75, lw=1.0)
        ax.set_ylabel(f"dxy {p1} vs {p2}")
    axes[-1].set_xlabel("Position on FR997704.1 (Mb)")
    axes[0].set_title(title)
    axes[0].legend(loc="upper right", fontsize=8, ncol=2)
    plt.tight_layout(); plt.savefig(out_png, dpi=140); plt.close()

def plot_combined_tracks(pi, fst, dxy, out_png, title, hud=None, depth=None):
    """Stack the windowed measures sharing the x-axis so they can be read against
    each other along the chromosome. One line per method. π is averaged across the
    two populations to keep one line per method. If `hud` (Hudson Fst) is given, a
    separate Hudson-Fst panel is added below the default (Weir & Cockerham) Fst.
    If `depth` (per-window mean depth) is given, it is added as the TOP panel — so
    you can see whether π/Fst/dxy signals coincide with depth anomalies."""
    if pi.empty or fst.empty or dxy.empty: return
    have_hud = hud is not None and not hud.empty and "avg_hudson_fst" in hud.columns
    have_dep = depth is not None and not depth.empty and "mean_depth" in depth.columns
    npan = 3 + int(have_hud) + int(have_dep)
    fig, axes = plt.subplots(npan, 1, figsize=(13, 3.0*npan), sharex=True)
    k = 0
    # depth (TOP panel, optional) — QC context for the popgen signals below
    if have_dep:
        for m,s in depth.groupby("method"):
            s = s.sort_values("window_pos_1")
            axes[k].plot(s["window_pos_1"]/1e6, s["mean_depth"], label=m,
                         color=METHOD_COLORS.get(m,"k"), ls=method_ls(m), alpha=0.75, lw=1.0)
        axes[k].set_ylabel("mean depth (×)"); k+=1
    # π (mean over populations, per method)
    pim = (pi.groupby(["method","window_pos_1"], as_index=False)["avg_pi"].mean())
    for m,s in pim.groupby("method"):
        s = s.sort_values("window_pos_1")
        axes[k].plot(s["window_pos_1"]/1e6, s["avg_pi"], label=m,
                     color=METHOD_COLORS.get(m,"k"), ls=method_ls(m), alpha=0.75, lw=1.0)
    axes[k].set_ylabel("π (mean A,P)"); k+=1
    # Fst — Weir & Cockerham (default)
    for m,s in fst.groupby("method"):
        s = s.sort_values("window_pos_1")
        axes[k].plot(s["window_pos_1"]/1e6, s["avg_wc_fst"], label=m,
                     color=METHOD_COLORS.get(m,"k"), ls=method_ls(m), alpha=0.75, lw=1.0)
    axes[k].set_ylabel("Fst W&C (A vs P)"); axes[k].axhline(0, color="grey", lw=0.5, ls=":"); k+=1
    # Fst — Hudson (optional)
    if have_hud:
        for m,s in hud.groupby("method"):
            s = s.sort_values("window_pos_1")
            axes[k].plot(s["window_pos_1"]/1e6, s["avg_hudson_fst"], label=m,
                         color=METHOD_COLORS.get(m,"k"), ls=method_ls(m), alpha=0.75, lw=1.0)
        axes[k].set_ylabel("Fst Hudson (A vs P)"); axes[k].axhline(0, color="grey", lw=0.5, ls=":"); k+=1
    # dxy (A vs P)
    for m,s in dxy.groupby("method"):
        s = s.sort_values("window_pos_1")
        axes[k].plot(s["window_pos_1"]/1e6, s["avg_dxy"], label=m,
                     color=METHOD_COLORS.get(m,"k"), ls=method_ls(m), alpha=0.75, lw=1.0)
    axes[k].set_ylabel("dxy (A vs P)")
    axes[k].set_xlabel("Position on FR997704.1 (Mb)")
    axes[0].set_title(title)
    axes[0].legend(loc="upper right", fontsize=8, ncol=3)
    plt.tight_layout(); plt.savefig(out_png, dpi=140); plt.close()

def correlation_heatmap(df, value_col, key_cols, out_png, title):
    if df.empty: return
    pivot = df.pivot_table(index=key_cols, columns="method", values=value_col)
    corr = pivot.corr(method="spearman")
    fig,ax = plt.subplots(figsize=(7,6))
    sns.heatmap(corr, annot=True, fmt=".3f", cmap="RdYlGn",
                vmin=0.7, vmax=1.0, ax=ax)
    plt.title(title); plt.xticks(rotation=40, ha="right"); plt.tight_layout()
    plt.savefig(out_png, dpi=140); plt.close()
    return corr

def outlier_jaccard(df, value_col, key_cols, q=0.99):
    if df.empty: return pd.DataFrame()
    pivot = df.pivot_table(index=key_cols, columns="method", values=value_col)
    out = {m: set(pivot.index[pivot[m] >= pivot[m].quantile(q)]) for m in pivot.columns}
    M = list(out)
    J = pd.DataFrame(index=M, columns=M, dtype=float)
    for a in M:
        for b in M:
            u = len(out[a]|out[b]); J.loc[a,b] = len(out[a]&out[b])/u if u else np.nan
    return J


# --- Tier 3: bcftools csq consequence summary ----------------------------
# Format per line:  CHR  POS  REF  ALT  BCSQ
#   BCSQ examples:
#     missense|GENE|TX|protein_coding|+|29Y>29I|28791T>A
#     intron|GENE||protein_coding
#     @28791            (alias: this row belongs to a multi-base csq starting at 28791)
#     .                 (no overlap)
def parse_vep_summary(vep_dir):
    rows = []
    for f in sorted(glob.glob(os.path.join(vep_dir,"*.csq.tsv"))):
        method = shortname(f).replace(".norm","")
        counts = {}
        with open(f) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 5: continue
                bcsq = p[4]
                if not bcsq or bcsq == "." or bcsq.startswith("@"):
                    cat = "intergenic_or_no_overlap" if bcsq in (".","") else "linked_to_prior"
                    counts[cat] = counts.get(cat,0)+1; continue
                # BCSQ can contain multiple annotations comma-separated
                for ann in bcsq.split(","):
                    if not ann or ann.startswith("@"):
                        counts["linked_to_prior"] = counts.get("linked_to_prior",0)+1; continue
                    cat = ann.split("|",1)[0].lstrip("*") or "other"
                    counts[cat] = counts.get(cat,0)+1
        for cat,cnt in counts.items():
            rows.append({"method":method,"consequence":cat,"count":int(cnt)})
    return pd.DataFrame(rows)

# Ordered buckets + colours so the legend reads top→bottom in the same order as
# the stack, and so the two plots share a colour scheme.
CSQ_ORDER = ["intergenic / no annotated feature","intronic","UTR",
             "synonymous","splice_region","missense",
             "splice_donor / acceptor","stop/start gained-lost","non-coding RNA","other"]
CSQ_COLORS = {
    "intergenic / no annotated feature":"#bdc3c7",  # grey — outside any gene
    "intronic":"#85c1e9",                            # light blue — inside gene, non-coding
    "UTR":"#f7dc6f",
    "synonymous":"#82e0aa",
    "splice_region":"#f0b27a",
    "missense":"#e74c3c",
    "splice_donor / acceptor":"#cb4335",
    "stop/start gained-lost":"#7b241c",
    "non-coding RNA":"#bb8fce",
    "other":"#566573",
}

def bucket_csq(c):
    """Map a bcftools-csq first-token to a broad, clearly-named class.
    NOTE: 'intronic' = variant inside an annotated gene's intron. It does NOT
    include intergenic variants — those have no BCSQ tag and are their own class."""
    c = c.lower()
    if c == "intergenic_or_no_overlap":           return "intergenic / no annotated feature"
    if "intron" in c:                              return "intronic"
    if "utr" in c:                                 return "UTR"
    if "synonymous" in c:                          return "synonymous"
    if "splice_donor" in c or "splice_acceptor" in c: return "splice_donor / acceptor"
    if "splice_region" in c:                       return "splice_region"
    if "missense" in c:                            return "missense"
    if any(k in c for k in ["stop_","start_","frameshift","initiator","inframe"]):
                                                   return "stop/start gained-lost"
    if "non_coding" in c or "nc_transcript" in c:  return "non-coding RNA"
    if "linked_to_prior" in c:                     return "intronic"  # rare; phased neighbour
    return "other"

def _stacked_pct(pivot, out_png, title, ylab="% of all SNPs"):
    cols = [c for c in CSQ_ORDER if c in pivot.columns] + \
           [c for c in pivot.columns if c not in CSQ_ORDER]
    pivot = pivot[cols]
    pct = pivot.div(pivot.sum(axis=1), axis=0) * 100
    ax = pct.plot(kind="bar", stacked=True, figsize=(11,6),
                  color=[CSQ_COLORS.get(c,"#566573") for c in cols])
    ax.set_ylabel(ylab); ax.set_title(title)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha="right")
    # reverse legend so it matches the visual stacking order
    h,l = ax.get_legend_handles_labels()
    ax.legend(h[::-1], l[::-1], bbox_to_anchor=(1.02,1), loc="upper left", fontsize=9,
              title="consequence class")
    plt.tight_layout(); plt.savefig(out_png, dpi=140); plt.close()

def plot_vep(df_vep, out_full, out_coding):
    """Two plots: (1) ALL classes incl. intergenic; (2) coding/splice only
    (intron + intergenic removed) so the rare functional classes are visible."""
    if df_vep.empty: return None
    df_vep = df_vep.copy()
    df_vep["bucket"] = df_vep["consequence"].map(bucket_csq)
    g = df_vep.groupby(["method","bucket"], as_index=False)["count"].sum()
    pivot = g.pivot(index="method", columns="bucket", values="count").fillna(0)

    # Plot 1 — everything, % of all SNPs
    _stacked_pct(pivot, out_full,
                 "Genomic context of SNPs by method (all classes)",
                 "% of all SNPs")

    # Plot 2 — drop the two non-coding giants to reveal coding/splice classes
    drop = ["intergenic / no annotated feature","intronic"]
    coding = pivot.drop(columns=[c for c in drop if c in pivot.columns])
    _stacked_pct(coding, out_coding,
                 "Coding & splice consequences only (intron + intergenic removed)",
                 "% of coding/splice SNPs")
    return pivot


# --- norm vs raw: what does normalisation actually do? -------------------
def plot_norm_vs_raw(root, plots):
    """Compare raw vs normalised within each method: record/SNP/multiallelic
    counts (bar) and per-window pi (scatter, near-identity expected)."""
    # (a) counts from bcftools stats
    def load_counts(regime):
        rows=[]
        for f in sorted(glob.glob(str(root/f"stats_{regime}"/"*.stats"))):
            d = parse_bcftools_stats(f); d["method"]=shortname(f); d["regime"]=regime
            rows.append(d)
        return pd.DataFrame(rows)
    raw = load_counts("raw"); norm = load_counts("norm")
    if raw.empty or norm.empty: return
    m = raw.merge(norm, on="method", suffixes=("_raw","_norm"))
    metrics = [("number_of_snps","SNP records"),
               ("number_of_multiallelic_sites","multi-allelic sites"),
               ("number_of_indels","indel records")]
    fig, axes = plt.subplots(1, len(metrics), figsize=(15,5))
    x = np.arange(len(m)); w=0.38
    for ax,(col,label) in zip(axes, metrics):
        if f"{col}_raw" not in m or f"{col}_norm" not in m:
            ax.set_visible(False); continue
        ax.bar(x-w/2, m[f"{col}_raw"],  w, label="raw",        color="#aab7c4")
        ax.bar(x+w/2, m[f"{col}_norm"], w, label="normalised", color="#1b4f72")
        ax.set_xticks(x); ax.set_xticklabels(m["method"], rotation=40, ha="right", fontsize=8)
        ax.set_title(label); ax.legend()
    fig.suptitle("What normalisation does: counts before vs after  bcftools norm -m -any")
    plt.tight_layout(); plt.savefig(plots/"norm_vs_raw_counts.png", dpi=140); plt.close()

    # (b) per-window pi raw vs norm (one point per window per population per method)
    pr = load_pixy(root/"pixy_raw","pi"); pn = load_pixy(root/"pixy_norm","pi")
    if not pr.empty and not pn.empty:
        key=["method","pop","window_pos_1"]
        j = pr.merge(pn, on=key, suffixes=("_raw","_norm"))
        j = j[(j["avg_pi_raw"]>0)&(j["avg_pi_norm"]>0)].copy()
        j["signed_diff"]=j["avg_pi_norm"]-j["avg_pi_raw"]
        j["rel_pct"]=100*j["signed_diff"]/j["avg_pi_raw"]
        mean_rel = float(j["rel_pct"].mean())
        fig,ax = plt.subplots(figsize=(6.8,6.8))
        for meth,s in j.groupby("method"):
            ax.scatter(s["avg_pi_raw"], s["avg_pi_norm"], s=6, alpha=0.4,
                       color=METHOD_COLORS.get(meth,"k"), label=meth)
        lim = [0, max(j["avg_pi_raw"].max(), j["avg_pi_norm"].max())*1.05]
        ax.plot(lim, lim, "k--", lw=1.2, label="y = x (no change)")
        ax.set_xlim(lim); ax.set_ylim(lim)
        ax.set_xlabel("π per 50 kb window — raw VCF")
        ax.set_ylabel("π per 50 kb window — normalised VCF")
        ax.set_title("Normalisation systematically inflates π by ~%.1f%%\n"
                     "(points sit above y=x: norm -m -any splits multi-allelic\n"
                     "sites, which pixy then counts per allele)" % mean_rel)
        ax.legend(fontsize=7, markerscale=2, loc="lower right")
        plt.tight_layout(); plt.savefig(plots/"norm_vs_raw_pi.png", dpi=140); plt.close()
        j["abs_diff"]=j["signed_diff"].abs()
        out = {"pi_mean_abs_diff": float(j["abs_diff"].mean()),
               "pi_max_abs_diff": float(j["abs_diff"].max()),
               "pi_mean_relative_pct": mean_rel,
               "pi_pct_windows_norm_higher": float(100*(j["signed_diff"]>0).mean())}
    else:
        out = {}

    # (c) per-window dxy raw vs norm — does normalisation inflate dxy too?
    dr = load_pixy(root/"pixy_raw","dxy"); dn = load_pixy(root/"pixy_norm","dxy")
    if not dr.empty and not dn.empty:
        key=["method","pop1","pop2","window_pos_1"]
        jd = dr.merge(dn, on=key, suffixes=("_raw","_norm"))
        jd = jd[(jd["avg_dxy_raw"]>0)&(jd["avg_dxy_norm"]>0)].copy()
        jd["sd"]=jd["avg_dxy_norm"]-jd["avg_dxy_raw"]
        mean_rel_d = float((100*jd["sd"]/jd["avg_dxy_raw"]).mean())
        fig,ax = plt.subplots(figsize=(6.8,6.8))
        for meth,s in jd.groupby("method"):
            ax.scatter(s["avg_dxy_raw"], s["avg_dxy_norm"], s=6, alpha=0.4,
                       color=METHOD_COLORS.get(meth,"k"), label=meth)
        lim = [0, max(jd["avg_dxy_raw"].max(), jd["avg_dxy_norm"].max())*1.05]
        ax.plot(lim, lim, "k--", lw=1.2, label="y = x (no change)")
        ax.set_xlim(lim); ax.set_ylim(lim)
        ax.set_xlabel("dxy per 50 kb window — raw VCF")
        ax.set_ylabel("dxy per 50 kb window — normalised VCF")
        ax.set_title("Normalisation inflates dxy by ~%.1f%% too\n"
                     "(same multi-allelic-splitting mechanism as π)" % mean_rel_d)
        ax.legend(fontsize=7, markerscale=2, loc="lower right")
        plt.tight_layout(); plt.savefig(plots/"norm_vs_raw_dxy.png", dpi=140); plt.close()
        out["dxy_mean_relative_pct"] = mean_rel_d
        out["dxy_pct_windows_norm_higher"] = float(100*(jd["sd"]>0).mean())
    return out or None


# --- main ----------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="compare_out")
    args = ap.parse_args()
    root = Path(args.root); plots = root/"plots"
    plots.mkdir(parents=True, exist_ok=True)

    summary = {}

    # Tier 1
    for regime in ("raw","norm"):
        d = collect_stats(root/f"stats_{regime}", regime)
        if not d.empty:
            keep = ["method","regime","number_of_records","number_of_snps","number_of_indels",
                    "number_of_multiallelic_sites","number_of_multiallelic_snp_sites"]
            keep = [c for c in keep if c in d.columns]
            d[keep].to_csv(plots/f"stats_summary_{regime}.tsv", sep="\t", index=False)
            summary[f"stats_{regime}"] = d[keep].to_dict("records")

    # Tier 2
    for regime in ("raw","norm"):
        sites = root/f"isec_{regime}"/"sites.txt"
        if sites.exists():
            df, M = parse_isec_sites(sites)
            upset_plot(M, plots/f"isec_sites_{regime}.png",
                       f"Site overlap across methods — {regime}  (n={len(M):,})")
            J = jaccard_heatmap(M, plots/f"isec_jaccard_{regime}.png",
                                f"Pairwise Jaccard of site sets — {regime}")
            J.to_csv(plots/f"isec_jaccard_{regime}.tsv", sep="\t")
            summary[f"isec_{regime}_n_sites"] = int(len(M))
            summary[f"isec_{regime}_n_in_all_six"] = int(M.all(axis=1).sum())

    # Tier 4
    for regime in ("raw","norm"):
        pi = load_pixy(root/f"pixy_{regime}", "pi")
        if not pi.empty:
            plot_pi_tracks(pi, plots/f"pi_per_method_{regime}.png",
                           f"Nucleotide diversity π (50 kb windows) — {regime}")
            corr = correlation_heatmap(pi, "avg_pi", ["pop","window_pos_1"],
                                       plots/f"pi_correlation_{regime}.png",
                                       f"Spearman corr. of windowed π across methods — {regime}")
            if corr is not None: corr.to_csv(plots/f"pi_correlation_{regime}.tsv", sep="\t")
        fst = load_pixy(root/f"pixy_{regime}", "fst")
        if not fst.empty and "avg_wc_fst" in fst.columns:
            plot_fst_tracks(fst, plots/f"fst_per_method_{regime}.png",
                            f"Fst (50 kb windows) — {regime}")
            corr = correlation_heatmap(fst, "avg_wc_fst", ["pop1","pop2","window_pos_1"],
                                       plots/f"fst_correlation_{regime}.png",
                                       f"Spearman corr. of windowed Fst across methods — {regime}")
            if corr is not None: corr.to_csv(plots/f"fst_correlation_{regime}.tsv", sep="\t")
            J = outlier_jaccard(fst, "avg_wc_fst", ["pop1","pop2","window_pos_1"])
            J.to_csv(plots/f"fst_outlier_jaccard_{regime}.tsv", sep="\t")

        dxy = load_pixy(root/f"pixy_{regime}", "dxy")
        if not dxy.empty and "avg_dxy" in dxy.columns:
            plot_dxy_tracks(dxy, plots/f"dxy_per_method_{regime}.png",
                            f"dxy — between-population divergence (50 kb windows) — {regime}")
            corr = correlation_heatmap(dxy, "avg_dxy", ["pop1","pop2","window_pos_1"],
                                       plots/f"dxy_correlation_{regime}.png",
                                       f"Spearman corr. of windowed dxy across methods — {regime}")
            if corr is not None: corr.to_csv(plots/f"dxy_correlation_{regime}.tsv", sep="\t")

        # combined genome scan (depth / π / Fst-WC / Fst-Hudson / dxy stacked, shared x-axis)
        hud = load_pixy(root/f"pixy_hudson_{regime}", "fst")
        dep_tsv = plots/"depth_per_window.tsv"   # regime-independent (from raw all-sites)
        depth = pd.read_csv(dep_tsv, sep="\t") if dep_tsv.exists() else None
        if not pi.empty and not fst.empty and not dxy.empty:
            plot_combined_tracks(pi, fst, dxy,
                                 plots/f"combined_tracks_{regime}.png",
                                 f"depth, π, Fst (W&C + Hudson) and dxy along FR997704.1 (50 kb windows) — {regime}",
                                 hud=hud, depth=depth)

    # Tier 3
    df_vep = parse_vep_summary(root/"vep")
    if not df_vep.empty:
        df_vep.to_csv(plots/"vep_consequence_by_method.tsv", sep="\t", index=False)
        pivot = plot_vep(df_vep,
                         plots/"vep_consequence_by_method.png",
                         plots/"vep_consequence_coding_only.png")
        if pivot is not None:
            summary["vep_buckets"] = pivot.to_dict()

    # Norm vs raw: what does normalisation do?
    nv = plot_norm_vs_raw(root, plots)
    if nv: summary["norm_vs_raw_pi"] = nv

    with open(plots/"summary.json","w") as f:
        json.dump(summary, f, indent=2, default=str)
    print("Analysis complete. Outputs in", plots)

if __name__ == "__main__":
    main()
