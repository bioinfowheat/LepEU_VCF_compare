"""
Integrative analyses (v3 additions):

  A. One-vs-five pi scatter grid: pick a reference method, plot its windowed pi
     against each of the other five, within raw and within normalised.
     -> pi_method_compare_raw.png, pi_method_compare_norm.png

  B. Are method-DISCORDANT SNPs enriched in particular consequence classes?
     Classify every normalised union SNP by its isec sharing pattern
     (shared-by-all-6 / BWA-only / NGM-only / singleton / other), look up the
     consequence class of each site, and test enrichment vs the shared core.
     -> enrichment_consequence.png  + enrichment_consequence.tsv

  C. Ts/Tv ratio by sharing category (a SNP-quality signal: real SNPs ~2-3,
     artefacts trend to 0.5). Computed straight from REF/ALT in sites.txt.
     -> tstv_by_sharing.png + tstv_by_sharing.tsv
"""
import os, glob, re, json
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/Users/stockholmbutterflylab/sbl_claudecode/VCF_compare")
OUT  = ROOT/"compare_out"; PLOTS = OUT/"plots"
PLOTS.mkdir(exist_ok=True, parents=True)

METHODS = ["bwa_cohort","bwa_hwe","bwa_nohwe","ngm_cohort","ngm_hwe","ngm_nohwe"]
MCOL = {"bwa_cohort":"#1f77b4","bwa_hwe":"#4a90c2","bwa_nohwe":"#7eb0d5",
        "ngm_cohort":"#6c3483","ngm_hwe":"#8e44ad","ngm_nohwe":"#b07cc6"}

# ---- shared csq bucketing (kept consistent with analyze_and_plot.py) -------
SEVERITY = ["stop/start gained-lost","splice_donor / acceptor","missense",
            "splice_region","synonymous","UTR","non-coding RNA",
            "intronic","intergenic / no annotated feature"]
SEV_RANK = {c:i for i,c in enumerate(SEVERITY)}  # lower = more severe

def bucket_csq(c):
    c=c.lower()
    if c in (".",""):                                 return "intergenic / no annotated feature"
    if "intron" in c:                                 return "intronic"
    if "utr" in c:                                     return "UTR"
    if "synonymous" in c:                             return "synonymous"
    if "splice_donor" in c or "splice_acceptor" in c: return "splice_donor / acceptor"
    if "splice_region" in c:                          return "splice_region"
    if "missense" in c:                               return "missense"
    if any(k in c for k in ["stop_","start_","frameshift","initiator","inframe"]):
                                                      return "stop/start gained-lost"
    if "non_coding" in c or "nc_transcript" in c:     return "non-coding RNA"
    if c.startswith("@"):                             return "intronic"
    return "intergenic / no annotated feature"

def shortname(path):
    m=re.search(r"Pnapi_(\w+)_bcftools_(\w+)_FR997704", os.path.basename(path))
    return f"{m.group(1)}_{m.group(2)}" if m else os.path.basename(path)

# ---------------------------------------------------------------------------
# A. one-vs-five pi scatter
# ---------------------------------------------------------------------------
def load_pixy(regime, stat="pi"):
    frames=[]
    for f in glob.glob(str(OUT/f"pixy_{regime}"/"**"/f"*_{stat}.txt"), recursive=True):
        d=pd.read_csv(f, sep="\t"); d["method"]=shortname(f); frames.append(d)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

def pi_one_vs_rest(regime, reference="bwa_cohort"):
    pi = load_pixy(regime,"pi")
    if pi.empty: return
    others=[m for m in METHODS if m!=reference]
    wide = pi.pivot_table(index=["pop","window_pos_1"], columns="method", values="avg_pi")
    wide = wide.dropna(subset=[reference])
    fig, axes = plt.subplots(1,5, figsize=(20,4.3), sharex=True, sharey=True)
    lim=[0, np.nanmax(wide.values)*1.05]
    for ax,m in zip(axes, others):
        sub=wide.dropna(subset=[m])
        ax.scatter(sub[reference], sub[m], s=7, alpha=0.4, color=MCOL[m])
        ax.plot(lim,lim,"k--",lw=1)
        # robust slope through origin + spearman
        x=sub[reference].values; y=sub[m].values
        slope = np.sum(x*y)/np.sum(x*x) if np.sum(x*x)>0 else np.nan
        rho = pd.Series(x).corr(pd.Series(y), method="spearman")
        ax.set_title(f"{m}\nslope={slope:.3f}  ρ={rho:.3f}", fontsize=9)
        ax.set_xlabel(f"π — {reference}"); ax.set_xlim(lim); ax.set_ylim(lim)
    axes[0].set_ylabel("π — other method")
    fig.suptitle(f"Windowed π: {reference} vs the other five methods  ({regime})  "
                 f"— points above y=x mean the other method calls higher π", fontsize=12)
    plt.tight_layout(); plt.savefig(PLOTS/f"pi_method_compare_{regime}.png", dpi=140); plt.close()

# ---------------------------------------------------------------------------
# B + C. sharing category -> consequence enrichment + Ts/Tv
# ---------------------------------------------------------------------------
TS = {("A","G"),("G","A"),("C","T"),("T","C")}  # transitions

def build_pos_consequence():
    """pos -> most severe consequence bucket, from the union of all 6 csq files."""
    pos2sev={}
    for f in glob.glob(str(OUT/"vep"/"*.csq.tsv")):
        with open(f) as fh:
            for line in fh:
                p=line.rstrip("\n").split("\t")
                if len(p)<5: continue
                pos=int(p[1]); bcsq=p[4]
                # take the most severe annotation listed
                best=None
                for ann in (bcsq.split(",") if bcsq not in (".","") else ["."]):
                    b=bucket_csq(ann.split("|",1)[0] if ann not in (".","") else ".")
                    if best is None or SEV_RANK[b]<SEV_RANK[best]: best=b
                if pos not in pos2sev or SEV_RANK[best]<SEV_RANK[pos2sev[pos]]:
                    pos2sev[pos]=best
    return pos2sev

def categorize(mask):
    """mask = 6-bit string [bwa_cohort,bwa_hwe,bwa_nohwe,ngm_cohort,ngm_hwe,ngm_nohwe]"""
    bits=[c=="1" for c in mask]
    n=sum(bits)
    if n==6:                       return "shared by all 6"
    if mask=="111000":             return "BWA-only (all 3)"
    if mask=="000111":             return "NGM-only (all 3)"
    if n==1:                       return "singleton (1 method)"
    if all(bits[:3]) and not any(bits[3:]): return "BWA-only (all 3)"
    if all(bits[3:]) and not any(bits[:3]): return "NGM-only (all 3)"
    return "other partial"

def sharing_analysis():
    sites = OUT/"isec_norm"/"sites.txt"
    if not sites.exists(): return
    pos2sev = build_pos_consequence()
    rows=[]
    with open(sites) as fh:
        for line in fh:
            p=line.rstrip("\n").split("\t")
            if len(p)<5: continue
            pos=int(p[1]); ref=p[2]; alt=p[3]; mask=p[4]
            cat=categorize(mask)
            cons=pos2sev.get(pos,"intergenic / no annotated feature")
            tstv = "ts" if (ref,alt) in TS else ("tv" if len(ref)==1 and len(alt)==1 else "other")
            rows.append((cat,cons,tstv))
    df=pd.DataFrame(rows, columns=["category","consequence","tstv"])

    # ---- B. consequence enrichment vs the shared core ----
    order=["shared by all 6","BWA-only (all 3)","NGM-only (all 3)","singleton (1 method)","other partial"]
    ct = pd.crosstab(df["category"], df["consequence"])
    ct = ct.reindex([c for c in order if c in ct.index])
    frac = ct.div(ct.sum(axis=1), axis=0)            # within-category proportions
    base = frac.loc["shared by all 6"]               # baseline = the core
    enr = np.log2(frac.div(base).replace(0,np.nan))  # log2 fold-change vs core
    cols=[c for c in SEVERITY if c in enr.columns]
    enr=enr[cols]; frac=frac[cols]
    enr.to_csv(PLOTS/"enrichment_consequence.tsv", sep="\t")
    frac.to_csv(PLOTS/"enrichment_consequence_fractions.tsv", sep="\t")

    fig,ax=plt.subplots(figsize=(11,4.8))
    im=ax.imshow(enr.values, cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=35, ha="right", fontsize=9)
    ax.set_yticks(range(len(enr.index))); ax.set_yticklabels(enr.index, fontsize=10)
    for i in range(enr.shape[0]):
        for j in range(enr.shape[1]):
            v=enr.values[i,j]
            if not np.isnan(v):
                ax.text(j,i,f"{v:+.1f}",ha="center",va="center",fontsize=8,
                        color="white" if abs(v)>1.2 else "black")
    cb=plt.colorbar(im,ax=ax); cb.set_label("log2 fold-change vs shared core")
    ax.set_title("Are method-discordant SNPs enriched in any consequence class?\n"
                 "(red = over-represented vs the shared-by-all-6 core; blue = under)")
    plt.tight_layout(); plt.savefig(PLOTS/"enrichment_consequence.png", dpi=140); plt.close()

    # ---- C. Ts/Tv by sharing category ----
    tt = df[df["tstv"].isin(["ts","tv"])].groupby("category")["tstv"].value_counts().unstack().fillna(0)
    tt["ts/tv"]= tt.get("ts",0)/tt.get("tv",1)
    tt = tt.reindex([c for c in order if c in tt.index])
    tt.to_csv(PLOTS/"tstv_by_sharing.tsv", sep="\t")
    fig,ax=plt.subplots(figsize=(8.5,4.8))
    bars=ax.bar(range(len(tt)), tt["ts/tv"], color=["#27ae60","#2874a6","#6c3483","#c0392b","#7f8c8d"][:len(tt)])
    ax.set_xticks(range(len(tt))); ax.set_xticklabels(tt.index, rotation=20, ha="right")
    ax.axhline(0.5, ls="--", color="grey", lw=1)
    ax.text(len(tt)-0.5, 0.55, "Ts/Tv = 0.5 (random / no signal)", ha="right", fontsize=8, color="grey")
    ax.set_ylabel("Ts/Tv ratio")
    ax.set_title("SNP-quality signal by sharing category\n(real SNPs ~2–3; method-private artefacts trend toward 0.5)")
    for b,v in zip(bars, tt["ts/tv"]):
        ax.text(b.get_x()+b.get_width()/2, v+0.03, f"{v:.2f}", ha="center", fontsize=9)
    plt.tight_layout(); plt.savefig(PLOTS/"tstv_by_sharing.png", dpi=140); plt.close()

    # counts table for the report
    counts = df["category"].value_counts().reindex([c for c in order if c in df["category"].unique()])
    counts.to_csv(PLOTS/"sharing_category_counts.tsv", sep="\t", header=["n_sites"])
    return {"tstv": tt["ts/tv"].to_dict(), "category_counts": counts.to_dict()}

if __name__=="__main__":
    for regime in ("raw","norm"):
        pi_one_vs_rest(regime, reference="bwa_cohort")
    summ = sharing_analysis()
    print("integrative analysis done:", json.dumps(summ, indent=1, default=str) if summ else "no sites")
