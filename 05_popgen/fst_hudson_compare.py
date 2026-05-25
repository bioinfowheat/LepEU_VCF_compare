#!/usr/bin/env python3
"""
Compare the two pixy Fst estimators: Weir & Cockerham (default, 'wc') vs Hudson.
WC results: compare_out/pixy_{regime}/<m>/<m>_fst.txt        (col avg_wc_fst)
Hudson    : compare_out/pixy_hudson_{regime}/<m>/<m>_fst.txt (col avg_hudson_fst)

Outputs (compare_out/plots/):
  fst_wc_vs_hudson_{raw,norm}.png   per-method scatter, WC (x) vs Hudson (y)
  fst_hudson_per_method_{raw,norm}.png   Hudson windowed tracks (all methods)
  fst_hudson_correlation_{raw,norm}.png  Hudson cross-method Spearman matrix
  fst_estimator_summary.tsv         per method/regime: slope, Pearson, Spearman, mean diff
"""
import glob, os, re, json
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, seaborn as sns

ROOT = Path("/Users/stockholmbutterflylab/sbl_claudecode/VCF_compare")
OUT = ROOT/"compare_out"; PLOTS = OUT/"plots"; PLOTS.mkdir(parents=True, exist_ok=True)
METHODS = ["bwa_cohort","bwa_hwe","bwa_nohwe","ngm_cohort","ngm_hwe","ngm_nohwe"]
MCOL = {"bwa_cohort":"#1f77b4","bwa_hwe":"#4a90c2","bwa_nohwe":"#7eb0d5",
        "ngm_cohort":"#6c3483","ngm_hwe":"#8e44ad","ngm_nohwe":"#b07cc6"}

def shortname(p):
    m=re.search(r"Pnapi_(\w+)_bcftools_(\w+)_FR997704", os.path.basename(p))
    return f"{m.group(1)}_{m.group(2)}" if m else os.path.basename(p)

def load_fst(pixy_dir):
    fr=[]
    for f in glob.glob(str(pixy_dir/"**"/"*_fst.txt"), recursive=True):
        d=pd.read_csv(f, sep="\t"); d["method"]=shortname(f); fr.append(d)
    return pd.concat(fr, ignore_index=True) if fr else pd.DataFrame()

KEY=["method","pop1","pop2","window_pos_1"]
summary=[]

for regime in ("raw","norm"):
    wc = load_fst(OUT/f"pixy_{regime}")
    hud = load_fst(OUT/f"pixy_hudson_{regime}")
    if wc.empty or hud.empty:
        print(f"missing fst for {regime} (wc={wc.empty}, hud={hud.empty})"); continue
    hud_col = "avg_hudson_fst" if "avg_hudson_fst" in hud.columns else \
              [c for c in hud.columns if "fst" in c.lower()][0]
    j = wc[KEY+["avg_wc_fst"]].merge(hud[KEY+[hud_col]], on=KEY)
    j = j.dropna(subset=["avg_wc_fst", hud_col])

    # --- per-method scatter: WC vs Hudson ---
    fig, axes = plt.subplots(2,3, figsize=(15,9.5), sharex=True, sharey=True)
    lo = min(j["avg_wc_fst"].min(), j[hud_col].min())
    hi = max(j["avg_wc_fst"].max(), j[hud_col].max())
    lim=[lo-0.02, hi+0.02]
    for ax,m in zip(axes.ravel(), METHODS):
        s=j[j.method==m]; x=s["avg_wc_fst"].values; y=s[hud_col].values
        ax.scatter(x,y,s=8,alpha=.45,color=MCOL[m])
        ax.plot(lim,lim,"k--",lw=1.2)
        slope=np.sum(x*y)/np.sum(x*x) if np.sum(x*x)>0 else np.nan
        pear=pd.Series(x).corr(pd.Series(y))
        spear=pd.Series(x).corr(pd.Series(y),method="spearman")
        md=float((y-x).mean())
        summary.append({"regime":regime,"method":m,"slope_hud_vs_wc":slope,
                        "pearson":pear,"spearman":spear,"mean_hud_minus_wc":md})
        ax.set_title(f"{m}\nr={pear:.3f}  ρ={spear:.3f}  Δμ={md:+.4f}",fontsize=9)
        ax.set_xlim(lim); ax.set_ylim(lim)
    for ax in axes[1,:]: ax.set_xlabel("Weir & Cockerham Fst (per window)")
    for ax in axes[:,0]: ax.set_ylabel("Hudson Fst (per window)")
    fig.suptitle(f"Fst estimator comparison — Weir & Cockerham (x) vs Hudson (y) — {regime}",fontsize=13)
    plt.tight_layout(); plt.savefig(PLOTS/f"fst_wc_vs_hudson_{regime}.png",dpi=140); plt.close()

    # --- Hudson tracks (per method) ---
    pairs=hud[["pop1","pop2"]].drop_duplicates().values.tolist()
    fig,axes=plt.subplots(len(pairs),1,figsize=(12,3*len(pairs)),sharex=True)
    if len(pairs)==1: axes=[axes]
    for ax,(p1,p2) in zip(axes,pairs):
        sub=hud[(hud.pop1==p1)&(hud.pop2==p2)]
        for m,s in sub.groupby("method"):
            s=s.sort_values("window_pos_1")
            ax.plot(s["window_pos_1"]/1e6, s[hud_col], label=m, color=MCOL.get(m,"k"), alpha=.75, lw=1.0)
        ax.set_ylabel(f"Hudson Fst {p1} vs {p2}"); ax.axhline(0,color="grey",lw=0.5,ls=":")
    axes[-1].set_xlabel("Position on FR997704.1 (Mb)")
    axes[0].set_title(f"Hudson Fst (50 kb windows) — {regime}")
    axes[0].legend(fontsize=8,ncol=3,loc="upper right")
    plt.tight_layout(); plt.savefig(PLOTS/f"fst_hudson_per_method_{regime}.png",dpi=140); plt.close()

    # --- Hudson cross-method Spearman correlation ---
    piv=hud.pivot_table(index=["pop1","pop2","window_pos_1"],columns="method",values=hud_col)
    corr=piv.corr(method="spearman")
    fig,ax=plt.subplots(figsize=(7,6))
    sns.heatmap(corr,annot=True,fmt=".3f",cmap="RdYlGn",vmin=0.3,vmax=1.0,ax=ax)
    plt.title(f"Spearman corr. of windowed Hudson Fst across methods — {regime}")
    plt.xticks(rotation=40,ha="right"); plt.tight_layout()
    plt.savefig(PLOTS/f"fst_hudson_correlation_{regime}.png",dpi=140); plt.close()
    corr.to_csv(PLOTS/f"fst_hudson_correlation_{regime}.tsv",sep="\t")

pd.DataFrame(summary).to_csv(PLOTS/"fst_estimator_summary.tsv",sep="\t",index=False)
print("fst hudson comparison done"); print(json.dumps(summary,indent=1,default=str)[:1500])
