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

    fig, ax = plt.subplots(figsize=(13, 4.2))
    for m, s in df.groupby("method"):
        s = s.sort_values("window_pos_1")
        ax.plot(s["window_pos_1"]/1e6, s["mean_depth"], label=m,
                color=METHOD_COLORS.get(m,"k"), ls=method_ls(m), alpha=0.8, lw=1.0)
    ax.set_xlabel("Position on FR997704.1 (Mb)")
    ax.set_ylabel("mean depth per sample (×)")
    ax.set_title("Mean read depth in 50 kb windows (vcftools --site-mean-depth)")
    ax.legend(loc="upper right", fontsize=8, ncol=3)
    plt.tight_layout(); plt.savefig(PLOTS/"depth_per_method.png", dpi=140); plt.close()
    print("wrote depth_per_method.png + depth_per_window.tsv")

if __name__ == "__main__":
    main()
