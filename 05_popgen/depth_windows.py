#!/usr/bin/env python3
"""
Bin vcftools --site-mean-depth output into 50 kb windows (matching pixy) and plot
a per-window mean-depth track, one line per method.

Input : compare_out/depth/<method>.ldepth.mean  (CHROM POS MEAN_DEPTH VAR_DEPTH)
Output: compare_out/plots/depth_per_window.tsv   (method, window_pos_1, mean_depth, n_sites)
        compare_out/plots/depth_per_method.png    (sliding-window track)

NB: MEAN_DEPTH is the per-site mean across the 14 samples (incl. zero-depth samples),
so the window value is the average per-sample coverage in that window. Depth is a
property of the mapping, so within a mapper the cohort/hwe/nohwe lines are identical.
"""
import glob, os, re
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/Users/stockholmbutterflylab/sbl_claudecode/VCF_compare")
OUT = ROOT/"compare_out"; PLOTS = OUT/"plots"; PLOTS.mkdir(parents=True, exist_ok=True)
WIN = 50000
METHODS = ["bwa_cohort","bwa_hwe","bwa_nohwe","ngm_cohort","ngm_hwe","ngm_nohwe"]
METHOD_COLORS = {"bwa_cohort":"#1f77b4","bwa_hwe":"#4a90c2","bwa_nohwe":"#7eb0d5",
                 "ngm_cohort":"#6c3483","ngm_hwe":"#8e44ad","ngm_nohwe":"#b07cc6"}
def method_ls(m): return "--" if str(m).endswith("_hwe") else "-"
def shortname(p):
    m=re.search(r"Pnapi_(\w+)_bcftools_(\w+)_FR997704", os.path.basename(p))
    return f"{m.group(1)}_{m.group(2)}" if m else os.path.basename(p)

def bin_file(path):
    """Stream the (large) per-site file in chunks; accumulate sum + count per window."""
    sums, counts = {}, {}
    for chunk in pd.read_csv(path, sep="\t", usecols=["POS","MEAN_DEPTH"],
                             dtype={"POS":"int64","MEAN_DEPTH":"float64"}, chunksize=2_000_000):
        w = ((chunk["POS"]-1)//WIN)*WIN + 1
        g = chunk.groupby(w)["MEAN_DEPTH"]
        for win, ssum in g.sum().items():  sums[win]   = sums.get(win,0)+ssum
        for win, cnt  in g.count().items(): counts[win] = counts.get(win,0)+cnt
    rows = [{"window_pos_1":w, "mean_depth":sums[w]/counts[w], "n_sites":counts[w]}
            for w in sorted(sums)]
    return pd.DataFrame(rows)

def main():
    frames=[]
    for f in sorted(glob.glob(str(OUT/"depth"/"*.ldepth.mean"))):
        m = shortname(f)
        d = bin_file(f); d["method"]=m; frames.append(d)
        print(f"binned {m}: {len(d)} windows, mean depth {d['mean_depth'].mean():.2f}x")
    if not frames: print("no depth files"); return
    df = pd.concat(frames, ignore_index=True)
    df.to_csv(PLOTS/"depth_per_window.tsv", sep="\t", index=False)

    present = [m for m in METHODS if m in set(df["method"])] + \
              [m for m in df["method"].unique() if m not in METHODS]

    # (1) overlay — all files on one axis
    fig, ax = plt.subplots(figsize=(13, 4.2))
    for m in present:
        s = df[df["method"]==m].sort_values("window_pos_1")
        ax.plot(s["window_pos_1"]/1e6, s["mean_depth"], label=m,
                color=METHOD_COLORS.get(m,"k"), ls=method_ls(m), alpha=0.8, lw=1.0)
    ax.set_xlabel("Position on FR997704.1 (Mb)")
    ax.set_ylabel("mean depth per sample (×)")
    ax.set_title("Mean read depth in 50 kb windows (vcftools --site-mean-depth) — all files overlaid")
    ax.legend(loc="upper right", fontsize=8, ncol=3)
    plt.tight_layout(); plt.savefig(PLOTS/"depth_per_method.png", dpi=140); plt.close()

    # (2) faceted — one panel per FILE, so each is individually visible (no assumption)
    ncol = 3; nrow = int(np.ceil(len(present)/ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(16, 3.2*nrow), sharex=True, sharey=True)
    axes = np.atleast_2d(axes)
    means = {}
    for i, m in enumerate(present):
        ax = axes[i//ncol][i%ncol]
        s = df[df["method"]==m].sort_values("window_pos_1")
        mu = float(s["mean_depth"].mean()); means[m] = mu
        ax.plot(s["window_pos_1"]/1e6, s["mean_depth"],
                color=METHOD_COLORS.get(m,"k"), ls=method_ls(m), lw=1.0)
        ax.axhline(mu, color="grey", lw=0.7, ls=":")
        ax.set_title(f"{m}   (mean {mu:.2f}×)", fontsize=10)
    for j in range(len(present), nrow*ncol): axes[j//ncol][j%ncol].set_visible(False)
    for ax in axes[-1]: ax.set_xlabel("Position on FR997704.1 (Mb)")
    for r in range(nrow): axes[r][0].set_ylabel("mean depth (×)")
    fig.suptitle("Read depth per window — one panel per VCF (vcftools --site-mean-depth)", fontsize=13)
    plt.tight_layout(); plt.savefig(PLOTS/"depth_per_file.png", dpi=140); plt.close()

    # (3) genome-wide mean depth per file — bar (makes any between-file difference obvious)
    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    vals = [means[m] for m in present]
    ax.bar(range(len(present)), vals, color=[METHOD_COLORS.get(m,"k") for m in present])
    ax.set_xticks(range(len(present))); ax.set_xticklabels(present, rotation=30, ha="right")
    ax.set_ylabel("genome-wide mean depth per sample (×)")
    ax.set_title("Mean read depth per VCF file")
    for i,v in enumerate(vals): ax.text(i, v+0.05, f"{v:.2f}", ha="center", fontsize=9)
    plt.tight_layout(); plt.savefig(PLOTS/"depth_mean_per_file.png", dpi=140); plt.close()
    pd.DataFrame({"method":present, "genome_mean_depth":vals}).to_csv(
        PLOTS/"depth_mean_per_file.tsv", sep="\t", index=False)

    print("wrote depth_per_method.png, depth_per_file.png, depth_mean_per_file.png + tsvs")
    print("per-file genome mean depth:", {m: round(means[m],3) for m in present})

if __name__ == "__main__":
    main()
