#!/usr/bin/env python3
"""
Section 5/6 — what normalisation actually does (raw vs normalised, head-to-head).

Three figures:
  norm_vs_raw_counts.png          per-method record counts (SNP / multi-allelic /
                                  indel) before vs after  bcftools norm -m -any
  norm_vs_raw_pi.png              all-method pooled scatter of windowed pi (raw x,
                                  norm y) — reveals the systematic pi inflation
  pi_norm_vs_raw_per_method.png   the same, broken out into one panel per method

KEY FINDING this quantifies: splitting multi-allelic sites with -m -any is correct
for site comparison / annotation, but it INFLATES pi by ~12% when the split VCF is
fed to pixy (pixy counts each split bi-allelic record's differences separately).
Practical rule: use un-split all-sites VCFs for pi/dxy/Fst.

Usage:  python norm_vs_raw.py [--root compare_out]
"""
import argparse, glob, os, re
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

METHODS=["bwa_cohort","bwa_hwe","bwa_nohwe","ngm_cohort","ngm_hwe","ngm_nohwe"]
MCOL={"bwa_cohort":"#1f77b4","bwa_hwe":"#4a90c2","bwa_nohwe":"#7eb0d5",
      "ngm_cohort":"#6c3483","ngm_hwe":"#8e44ad","ngm_nohwe":"#b07cc6"}

def shortname(p):
    m=re.search(r"Pnapi_(\w+)_bcftools_(\w+)_FR997704", os.path.basename(p))
    return f"{m.group(1)}_{m.group(2)}" if m else os.path.basename(p)

def parse_stats(path):
    out={}
    for line in open(path):
        if line.startswith("SN\t"):
            p=line.strip().split("\t"); k=p[2].rstrip(":").replace(" ","_").lower()
            try: out[k]=float(p[3])
            except: out[k]=p[3]
    return out

def counts(root):
    def load(reg):
        rows=[]
        for f in sorted(glob.glob(str(root/f"stats_{reg}"/"*.stats"))):
            d=parse_stats(f); d["method"]=shortname(f); rows.append(d)
        return pd.DataFrame(rows)
    raw=load("raw"); norm=load("norm")
    if raw.empty or norm.empty: return
    m=raw.merge(norm,on="method",suffixes=("_raw","_norm"))
    metrics=[("number_of_snps","SNP records"),
             ("number_of_multiallelic_sites","multi-allelic sites"),
             ("number_of_indels","indel records")]
    fig,axes=plt.subplots(1,3,figsize=(15,5)); x=np.arange(len(m)); w=0.38
    for ax,(c,lab) in zip(axes,metrics):
        ax.bar(x-w/2,m[f"{c}_raw"],w,label="raw",color="#aab7c4")
        ax.bar(x+w/2,m[f"{c}_norm"],w,label="normalised",color="#1b4f72")
        ax.set_xticks(x); ax.set_xticklabels(m["method"],rotation=40,ha="right",fontsize=8)
        ax.set_title(lab); ax.legend()
    fig.suptitle("What normalisation does: counts before vs after  bcftools norm -m -any")
    plt.tight_layout(); plt.savefig(root/"plots"/"norm_vs_raw_counts.png",dpi=140); plt.close()

def load_pi(root,reg):
    fr=[]
    for f in glob.glob(str(root/f"pixy_{reg}"/"**"/"*_pi.txt"),recursive=True):
        d=pd.read_csv(f,sep="\t"); d["method"]=shortname(f); fr.append(d)
    return pd.concat(fr,ignore_index=True) if fr else pd.DataFrame()

def pi_pooled(root):
    pr=load_pi(root,"raw"); pn=load_pi(root,"norm")
    if pr.empty or pn.empty: return
    j=pr.merge(pn,on=["method","pop","window_pos_1"],suffixes=("_raw","_norm"))
    j=j[(j.avg_pi_raw>0)&(j.avg_pi_norm>0)]
    rel=float((100*(j.avg_pi_norm-j.avg_pi_raw)/j.avg_pi_raw).mean())
    lim=[0,max(j.avg_pi_raw.max(),j.avg_pi_norm.max())*1.05]
    fig,ax=plt.subplots(figsize=(6.8,6.8))
    for mth,s in j.groupby("method"):
        ax.scatter(s.avg_pi_raw,s.avg_pi_norm,s=6,alpha=.4,color=MCOL.get(mth,"k"),label=mth)
    ax.plot(lim,lim,"k--",lw=1.2,label="y = x")
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("π per 50 kb window — raw"); ax.set_ylabel("π per 50 kb window — normalised")
    ax.set_title(f"Normalisation systematically inflates π by ~{rel:.1f}%\n"
                 "(-m -any splits multi-allelics; pixy counts them per allele)")
    ax.legend(fontsize=7,markerscale=2,loc="lower right")
    plt.tight_layout(); plt.savefig(root/"plots"/"norm_vs_raw_pi.png",dpi=140); plt.close()

def pi_per_method(root):
    pr=load_pi(root,"raw"); pn=load_pi(root,"norm")
    if pr.empty or pn.empty: return
    j=pr.merge(pn,on=["method","pop","window_pos_1"],suffixes=("_raw","_norm"))
    j=j[(j.avg_pi_raw>0)&(j.avg_pi_norm>0)]
    lim=[0,max(j.avg_pi_raw.max(),j.avg_pi_norm.max())*1.05]
    fig,axes=plt.subplots(2,3,figsize=(15,9.5),sharex=True,sharey=True); stats={}
    for ax,mth in zip(axes.ravel(),METHODS):
        s=j[j.method==mth]; x=s.avg_pi_raw.values; y=s.avg_pi_norm.values
        ax.scatter(x,y,s=8,alpha=.45,color=MCOL[mth]); ax.plot(lim,lim,"k--",lw=1.2)
        slope=np.sum(x*y)/np.sum(x*x) if np.sum(x*x)>0 else np.nan
        rho=pd.Series(x).corr(pd.Series(y),method="spearman")
        rel=float((100*(y-x)/x).mean())
        stats[mth]={"slope":slope,"rho":rho,"mean_pct_inflation":rel}
        ax.set_title(f"{mth}\nslope={slope:.3f}  ρ={rho:.3f}  +{rel:.1f}% π",fontsize=10)
        ax.set_xlim(lim); ax.set_ylim(lim)
    for ax in axes[1,:]: ax.set_xlabel("π — raw")
    for ax in axes[:,0]: ax.set_ylabel("π — normalised")
    fig.suptitle("Normalisation effect on π, per method (points above y=x ⇒ inflation)",fontsize=13)
    plt.tight_layout(); plt.savefig(root/"plots"/"pi_norm_vs_raw_per_method.png",dpi=140); plt.close()
    pd.DataFrame(stats).T.to_csv(root/"plots"/"pi_norm_vs_raw_per_method.tsv",sep="\t")

def load_stat(root,reg,stat):
    fr=[]
    for f in glob.glob(str(root/f"pixy_{reg}"/"**"/f"*_{stat}.txt"),recursive=True):
        d=pd.read_csv(f,sep="\t"); d["method"]=shortname(f); fr.append(d)
    return pd.concat(fr,ignore_index=True) if fr else pd.DataFrame()

def dxy_pooled(root):
    """dxy raw vs norm — same multi-allelic-splitting inflation as pi."""
    dr=load_stat(root,"raw","dxy"); dn=load_stat(root,"norm","dxy")
    if dr.empty or dn.empty: return
    j=dr.merge(dn,on=["method","pop1","pop2","window_pos_1"],suffixes=("_raw","_norm"))
    j=j[(j.avg_dxy_raw>0)&(j.avg_dxy_norm>0)]
    rel=float((100*(j.avg_dxy_norm-j.avg_dxy_raw)/j.avg_dxy_raw).mean())
    lim=[0,max(j.avg_dxy_raw.max(),j.avg_dxy_norm.max())*1.05]
    fig,ax=plt.subplots(figsize=(6.8,6.8))
    for mth,s in j.groupby("method"):
        ax.scatter(s.avg_dxy_raw,s.avg_dxy_norm,s=6,alpha=.4,color=MCOL.get(mth,"k"),label=mth)
    ax.plot(lim,lim,"k--",lw=1.2,label="y = x"); ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("dxy per 50 kb window — raw"); ax.set_ylabel("dxy per 50 kb window — normalised")
    ax.set_title(f"Normalisation inflates dxy by ~{rel:.1f}% too\n(same mechanism as π)")
    ax.legend(fontsize=7,markerscale=2,loc="lower right")
    plt.tight_layout(); plt.savefig(root/"plots"/"norm_vs_raw_dxy.png",dpi=140); plt.close()

if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("--root",default="compare_out")
    a=ap.parse_args(); root=Path(a.root); (root/"plots").mkdir(parents=True,exist_ok=True)
    counts(root); pi_pooled(root); pi_per_method(root); dxy_pooled(root)
    print("norm_vs_raw done ->", root/"plots")
