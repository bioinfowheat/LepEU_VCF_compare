"""
Two v3 add-ons:
  1. 12x12 pairwise Jaccard of SNP-site sets across all 12 VCFs (6 raw + 6 norm),
     same framework as section 2c -> jaccard_all12.png / .tsv
  2. Per-method raw-vs-normalised windowed-pi scatter grid (6 panels) for section 4c
     -> pi_norm_vs_raw_per_method.png
"""
import os, glob, re, json
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

ROOT = Path("/Users/stockholmbutterflylab/sbl_claudecode/VCF_compare")
OUT  = ROOT/"compare_out"; PLOTS = OUT/"plots"

METHODS = ["bwa_cohort","bwa_hwe","bwa_nohwe","ngm_cohort","ngm_hwe","ngm_nohwe"]
MCOL = {"bwa_cohort":"#1f77b4","bwa_hwe":"#4a90c2","bwa_nohwe":"#7eb0d5",
        "ngm_cohort":"#6c3483","ngm_hwe":"#8e44ad","ngm_nohwe":"#b07cc6"}
LABELS12 = [f"{m}_raw"  for m in METHODS] + [f"{m}_norm" for m in METHODS]

def shortname(path):
    m=re.search(r"Pnapi_(\w+)_bcftools_(\w+)_FR997704", os.path.basename(path))
    return f"{m.group(1)}_{m.group(2)}" if m else os.path.basename(path)

# ---------------------------------------------------------------------------
# 1. 12x12 Jaccard
# ---------------------------------------------------------------------------
def jaccard_all12():
    f = OUT/"isec_all12"/"sites.txt"
    if not f.exists():
        print("missing isec_all12/sites.txt"); return
    masks = []
    with open(f) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 5 and len(p[4]) == 12:
                masks.append(p[4])
    A = np.array([[c == "1" for c in m] for m in masks], dtype=bool)  # n x 12
    n = A.shape[1]
    J = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            inter = np.count_nonzero(A[:, i] & A[:, j])
            uni   = np.count_nonzero(A[:, i] | A[:, j])
            J[i, j] = inter/uni if uni else np.nan
    Jdf = pd.DataFrame(J, index=LABELS12, columns=LABELS12)
    Jdf.to_csv(PLOTS/"jaccard_all12.tsv", sep="\t")

    fig, ax = plt.subplots(figsize=(11, 9.2))
    sns.heatmap(Jdf, annot=True, fmt=".2f", cmap="viridis", vmin=0.5, vmax=1.0,
                annot_kws={"size":7}, square=True, ax=ax,
                cbar_kws={"label":"Jaccard index  |A∩B| / |A∪B|"})
    # separate the raw (first 6) and norm (last 6) blocks
    for k in (6,):
        ax.axhline(k, color="white", lw=3); ax.axvline(k, color="white", lw=3)
    ax.set_title("Pairwise Jaccard of SNP-site sets — all 12 VCFs\n"
                 "(top-left 6×6 = raw, bottom-right 6×6 = normalised)")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=8)
    plt.tight_layout(); plt.savefig(PLOTS/"jaccard_all12.png", dpi=140); plt.close()
    # report a couple of useful summaries
    same_method = np.mean([Jdf.loc[f"{m}_raw", f"{m}_norm"] for m in METHODS])
    print(f"jaccard_all12 done. mean(raw vs norm, same method) = {same_method:.3f}")
    return {"mean_same_method_raw_vs_norm": float(same_method)}

# ---------------------------------------------------------------------------
# 2. per-method raw vs norm pi grid
# ---------------------------------------------------------------------------
def load_pixy(regime, stat="pi"):
    frames=[]
    for f in glob.glob(str(OUT/f"pixy_{regime}"/"**"/f"*_{stat}.txt"), recursive=True):
        d=pd.read_csv(f, sep="\t"); d["method"]=shortname(f); frames.append(d)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

def pi_norm_vs_raw_per_method():
    pr = load_pixy("raw","pi"); pn = load_pixy("norm","pi")
    if pr.empty or pn.empty: print("missing pixy"); return
    key=["method","pop","window_pos_1"]
    j = pr.merge(pn, on=key, suffixes=("_raw","_norm"))
    j = j[(j["avg_pi_raw"]>0)&(j["avg_pi_norm"]>0)].copy()
    lim=[0, max(j["avg_pi_raw"].max(), j["avg_pi_norm"].max())*1.05]
    fig, axes = plt.subplots(2, 3, figsize=(15, 9.5), sharex=True, sharey=True)
    stats={}
    for ax, m in zip(axes.ravel(), METHODS):
        s = j[j["method"]==m]
        ax.scatter(s["avg_pi_raw"], s["avg_pi_norm"], s=8, alpha=0.45, color=MCOL[m])
        ax.plot(lim, lim, "k--", lw=1.2)
        x=s["avg_pi_raw"].values; y=s["avg_pi_norm"].values
        slope = np.sum(x*y)/np.sum(x*x) if np.sum(x*x)>0 else np.nan
        rho = pd.Series(x).corr(pd.Series(y), method="spearman")
        rel = float((100*(y-x)/x).mean())
        higher = float(100*np.mean(y>x))
        stats[m]={"slope":float(slope),"rho":float(rho),"mean_pct_inflation":rel,"pct_windows_higher":higher}
        ax.set_title(f"{m}\nslope={slope:.3f}  ρ={rho:.3f}  +{rel:.1f}% π", fontsize=10)
        ax.set_xlim(lim); ax.set_ylim(lim)
    for ax in axes[1,:]: ax.set_xlabel("π per 50 kb window — raw")
    for ax in axes[:,0]: ax.set_ylabel("π per 50 kb window — normalised")
    fig.suptitle("Normalisation effect on π, per method (raw x-axis vs normalised y-axis)\n"
                 "points above y=x ⇒ normalisation inflates π", fontsize=13)
    plt.tight_layout(); plt.savefig(PLOTS/"pi_norm_vs_raw_per_method.png", dpi=140); plt.close()
    pd.DataFrame(stats).T.to_csv(PLOTS/"pi_norm_vs_raw_per_method.tsv", sep="\t")
    print("pi_norm_vs_raw_per_method done")
    return stats

if __name__=="__main__":
    s1 = jaccard_all12()
    s2 = pi_norm_vs_raw_per_method()
    print(json.dumps({"jaccard_all12":s1,"pi_per_method":s2}, indent=1, default=str))
