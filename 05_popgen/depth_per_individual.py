#!/usr/bin/env python3
"""
Per-INDIVIDUAL read-depth traces (50 kb windows), one PANEL per VCF (stacked, shared
axes). Default: bwa_hwe (top) and bwa_nohwe (bottom). Each of the 14 samples gets its
own colour, consistent across panels.

Depth comes from per-sample FORMAT/DP:  bcftools query -f '%CHROM\\t%POS[\\t%DP]\\n'
binned to 50 kb windows (mean DP per sample per window; missing '.' treated as 0).

Output: compare_out/plots/depth_per_individual.png  (+ depth_per_individual.tsv)
"""
import subprocess, io
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path("/Users/stockholmbutterflylab/sbl_claudecode/VCF_compare")
BIN  = "/Users/stockholmbutterflylab/miniforge3/envs/vcfcompare/bin"
VCFDIR = ROOT/"Pnapi_chr10_FR997704"; PLOTS = ROOT/"compare_out"/"plots"; PLOTS.mkdir(parents=True, exist_ok=True)
WIN = 50000

# (label, vcf path)  — one panel per VCF
VCFS = [
    ("bwa_hwe",   VCFDIR/"Pnapi_bwa_bcftools_hwe_FR997704.1.vcf.gz"),
    ("bwa_nohwe", VCFDIR/"Pnapi_bwa_bcftools_nohwe_FR997704.1.vcf.gz"),
]

def samples(vcf):
    return subprocess.run([f"{BIN}/bcftools","query","-l",str(vcf)],
                          capture_output=True, text=True, check=True).stdout.split()

def windowed_per_sample(vcf, samps):
    """Stream per-sample DP and bin to windows with awk; return long df."""
    awk = (r'{w=int(($2-1)/%d)*%d+1; for(i=3;i<=NF;i++){d=($i=="."?0:$i+0); '
           r'k=w SUBSEP (i-2); s[k]+=d; n[k]++}} '
           r'END{for(k in s){split(k,a,SUBSEP); print a[1]"\t"a[2]"\t"s[k]/n[k]}}') % (WIN, WIN)
    q = subprocess.Popen([f"{BIN}/bcftools","query","-f",
                          "%CHROM\\t%POS[\\t%DP]\\n", str(vcf)],
                         stdout=subprocess.PIPE)
    out = subprocess.run(["awk", awk], stdin=q.stdout, capture_output=True, text=True)
    q.wait()
    df = pd.read_csv(io.StringIO(out.stdout), sep="\t",
                     names=["window_pos_1","sidx","mean_depth"])
    df["sample"] = df["sidx"].map(lambda i: samps[int(i)-1])
    return df.drop(columns="sidx")

def main():
    samps = samples(VCFS[0][1])
    cmap = plt.get_cmap("tab20")
    scol = {s: cmap(i % 20) for i,s in enumerate(samps)}

    n = len(VCFS)
    fig, axes = plt.subplots(n, 1, figsize=(14, 3.6*n), sharex=True, sharey=True)
    axes = np.atleast_1d(axes)
    allrows = []
    for ax,(label, vcf) in zip(axes, VCFS):
        df = windowed_per_sample(vcf, samps); df["vcf"] = label; allrows.append(df)
        for s in samps:
            d = df[df["sample"]==s].sort_values("window_pos_1")
            ax.plot(d["window_pos_1"]/1e6, d["mean_depth"], color=scol[s], lw=0.9, alpha=0.85)
        ax.set_ylabel("mean depth in 50 kb (×)")
        ax.set_title(f"{label}   (overall mean {df['mean_depth'].mean():.2f}×, 14 individuals)", fontsize=11)
        print(f"{label}: {df['sample'].nunique()} samples, overall mean DP {df['mean_depth'].mean():.2f}")
    pd.concat(allrows, ignore_index=True).to_csv(PLOTS/"depth_per_individual.tsv", sep="\t", index=False)

    axes[-1].set_xlabel("Position on FR997704.1 (Mb)")
    fig.suptitle("Per-individual read-depth traces — one panel per VCF", fontsize=13)
    # shared legend (individuals) to the right
    axes[0].legend(handles=[Line2D([0],[0], color=scol[s], lw=2, label=s) for s in samps],
                   title="individual", ncol=1, fontsize=7, loc="upper left",
                   bbox_to_anchor=(1.01, 1.0))
    plt.tight_layout(); plt.savefig(PLOTS/"depth_per_individual.png", dpi=150, bbox_inches="tight"); plt.close()
    print("wrote depth_per_individual.png + .tsv")

if __name__ == "__main__":
    main()
