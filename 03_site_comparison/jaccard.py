#!/usr/bin/env python3
"""
Tier 2 — pairwise Jaccard similarity of SNP-site sets.

Jaccard(A,B) = |A ∩ B| / |A ∪ B|  = (sites called by both) / (sites called by either)
A scale-free measure (0 = disjoint, 1 = identical) of how much two methods agree on
*where* SNPs are. Computed from the bitmask in bcftools-isec sites.txt.

Produces:
  jaccard_raw.png/.tsv     6x6, raw VCFs            (from compare_out/isec_raw)
  jaccard_norm.png/.tsv    6x6, normalised VCFs     (from compare_out/isec_norm)
  jaccard_all12.png/.tsv   12x12, raw + normalised  (from compare_out/isec_all12)

Usage:  python jaccard.py [--root compare_out] [--outdir compare_out/plots]
"""
import argparse, numpy as np, pandas as pd
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, seaborn as sns

METHODS = ["bwa_cohort","bwa_hwe","bwa_nohwe","ngm_cohort","ngm_hwe","ngm_nohwe"]
LABELS12 = [f"{m}_raw" for m in METHODS] + [f"{m}_norm" for m in METHODS]

def read_mask_matrix(sites_txt, nbits):
    masks=[]
    with open(sites_txt) as fh:
        for line in fh:
            p=line.rstrip("\n").split("\t")
            if len(p)>=5 and len(p[4])==nbits:
                masks.append(p[4])
    return np.array([[c=="1" for c in m] for m in masks], dtype=bool)

def jaccard_matrix(A):
    n=A.shape[1]; J=np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            inter=np.count_nonzero(A[:,i]&A[:,j]); uni=np.count_nonzero(A[:,i]|A[:,j])
            J[i,j]=inter/uni if uni else np.nan
    return J

def plot(J, labels, out_png, title, block=None):
    df=pd.DataFrame(J, index=labels, columns=labels)
    sq = len(labels)<=6
    fig,ax=plt.subplots(figsize=(7,6) if sq else (11,9.2))
    sns.heatmap(df, annot=True, fmt=".2f" if not sq else ".3f", cmap="viridis",
                vmin=0.5, vmax=1.0, square=True, ax=ax,
                annot_kws={"size":7 if not sq else 9},
                cbar_kws={"label":"Jaccard index"})
    if block:
        ax.axhline(block, color="white", lw=3); ax.axvline(block, color="white", lw=3)
    ax.set_title(title)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=8)
    plt.tight_layout(); plt.savefig(out_png, dpi=140); plt.close()
    return df

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--root", default="compare_out")
    ap.add_argument("--outdir", default="compare_out/plots")
    a=ap.parse_args(); root=Path(a.root); out=Path(a.outdir); out.mkdir(parents=True, exist_ok=True)

    for regime in ("raw","norm"):
        f=root/f"isec_{regime}"/"sites.txt"
        if f.exists():
            A=read_mask_matrix(f,6); J=jaccard_matrix(A)
            df=plot(J, METHODS, out/f"jaccard_{regime}.png",
                    f"Pairwise Jaccard of SNP-site sets — {regime}")
            df.to_csv(out/f"jaccard_{regime}.tsv", sep="\t")

    f12=root/"isec_all12"/"sites.txt"
    if f12.exists():
        A=read_mask_matrix(f12,12); J=jaccard_matrix(A)
        df=plot(J, LABELS12, out/"jaccard_all12.png",
                "Pairwise Jaccard — all 12 VCFs (top-left 6×6 raw, bottom-right 6×6 norm)",
                block=6)
        df.to_csv(out/"jaccard_all12.tsv", sep="\t")
        same=np.mean([df.loc[f"{m}_raw",f"{m}_norm"] for m in METHODS])
        print(f"mean Jaccard(raw vs norm, same method) = {same:.3f}")
    print("jaccard done ->", out)

if __name__=="__main__":
    main()
