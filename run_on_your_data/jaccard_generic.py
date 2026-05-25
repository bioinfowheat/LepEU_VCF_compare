#!/usr/bin/env python3
"""
Pairwise Jaccard of SNP-site sets — generic (any number of VCFs).

Jaccard(A,B) = |A ∩ B| / |A ∪ B| = (sites called by both) / (sites called by either).
Reads the bcftools-isec sites.txt bitmask + a labels file (one label per line, in the
isec input order) and writes an NxN Jaccard heatmap + table.

Usage:
  python jaccard_generic.py --sites compare_out/isec/sites.txt \
                            --labels compare_out/labels.txt \
                            --outdir compare_out/plots
"""
import argparse
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, seaborn as sns

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sites", required=True)
    ap.add_argument("--labels", required=True)
    ap.add_argument("--outdir", default="compare_out/plots")
    a = ap.parse_args()
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)

    labels = [l.strip() for l in open(a.labels) if l.strip()]
    n = len(labels)

    # read the N-bit mask column (5th) into a boolean matrix
    masks = []
    with open(a.sites) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 5 and len(p[4]) == n:
                masks.append(p[4])
    if not masks:
        raise SystemExit(f"No rows with a {n}-bit mask found in {a.sites} "
                         f"(does labels.txt match the number of VCFs in the isec?)")
    A = np.array([[c == "1" for c in m] for m in masks], dtype=bool)

    J = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            inter = np.count_nonzero(A[:, i] & A[:, j])
            uni   = np.count_nonzero(A[:, i] | A[:, j])
            J[i, j] = inter / uni if uni else np.nan
    df = pd.DataFrame(J, index=labels, columns=labels)
    df.to_csv(out / "jaccard.tsv", sep="\t")

    annot = n <= 14                      # stop annotating when the grid gets dense
    fig, ax = plt.subplots(figsize=(max(7, 0.8*n), max(6, 0.7*n)))
    sns.heatmap(df, annot=annot, fmt=".2f", cmap="viridis", vmin=0.5, vmax=1.0,
                square=True, ax=ax, annot_kws={"size": 7},
                cbar_kws={"label": "Jaccard index  |A∩B| / |A∪B|"})
    ax.set_title(f"Pairwise Jaccard of SNP-site sets (n={n})")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=8)
    plt.tight_layout(); plt.savefig(out / "jaccard.png", dpi=140); plt.close()
    print(f"Wrote {out/'jaccard.png'} and {out/'jaccard.tsv'}  (n={n}, {A.shape[0]:,} sites)")

if __name__ == "__main__":
    main()
