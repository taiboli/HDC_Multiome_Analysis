"""
Collect S-LDSC results (17_sldsc.sh) for one or more batches.

Per trait x test annotation (row of the annotation in the .results of its model):
  Prop_SNPs, Prop_h2, Enrichment (+ se, p), tau (Coefficient, se, z), p_tau = one-sided P(tau > 0) test
  (the conditional test: contribution beyond baselineLD v2.2 + the batch's control annotation),
  tau_star = tau * M_ref * sd(annot) / h2g (Gazal et al. 2017; binary annot sd = sqrt(p (1 - p))).
Trait-level: h2g, h2 z from the baseline model (baselineLD only). Primary set: h2 z > 5 (standard S-LDSC filter).
LD reference: EUR (FinnGen, UKB EUR / meta), EAS (BBJ); PAGE and UKB AFR use EUR LD (ld_match = False).
FDR (BH) per batch x annotation family (peakspec, genespec, caQTLpeak, eGene, controls) over primary traits
with matched LD reference; p-values for other traits are reported without FDR.

Output: GWAS/ldsc/sldsc_<batch>.tsv  and GWAS/results/sldsc_<batch>.tsv (copy)
Usage: python 18_sldsc_summarize.py spec [qtl]
"""
import glob
import os
import re
import shutil
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm

BASE = os.environ["HDC_DIR"] + "/"
L = f"{BASE}/GWAS/ldsc/"
REF = os.environ["LDSC_REF_DIR"] + "/"
man = pd.read_csv(f"{BASE}/GWAS/gwas_manifest.tsv", sep="\t")[["trait_id", "source", "ancestry", "label", "type"]]


def ref_M(pop):
    d = REF + ("1000G_Phase3_baselineLD_v2.2_ldscores" if pop == "EUR" else "1000G_Phase3_EAS_baselineLD_v2.2_ldscores")
    suf = "M_5_50" if pop == "EUR" else "M"  # EAS run with --not-M-5-50
    return sum(float(open(f"{d}/baselineLD.{c}.l2.{suf}").read().split()[0]) for c in range(1, 23))


M = {"EUR": ref_M("EUR"), "EAS": ref_M("EAS")}


def h2_from_log(f):
    t = open(f).read()
    m = re.search(r"Total Observed scale h2: ([-\d.e]+) \(([-\d.e]+)\)", t)
    return (float(m.group(1)), float(m.group(2))) if m else (np.nan, np.nan)


def bh(p):
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out
    pp = p[ok]
    o = np.argsort(pp)
    q = np.minimum.accumulate((pp[o] * len(pp) / np.arange(1, len(pp) + 1))[::-1])[::-1]
    r = np.empty(len(pp))
    r[o] = np.minimum(q, 1)
    out[ok] = r
    return out


for batch in sys.argv[1:]:
    rows = []
    for tdir in sorted(glob.glob(f"{L}/results/{batch}/*/")):
        trait = tdir.rstrip("/").split("/")[-1]
        pop = "EAS" if trait.startswith("BBJ__") else "EUR"
        bl = f"{L}/results/spec/{trait}/baseline.log"
        h2, h2se = h2_from_log(bl) if os.path.exists(bl) else (np.nan, np.nan)
        for f in sorted(glob.glob(f"{tdir}/*.results")):
            name = os.path.basename(f)[:-8]
            if name == "baseline":
                continue
            r = pd.read_csv(f, sep="\t")
            r = r[r.Category.str.replace(r"L2_\d+$", "", regex=True) == name]
            if len(r) != 1:
                continue
            r = r.iloc[0]
            h2m, _ = h2_from_log(f[:-8] + ".log")
            p = r["Prop._SNPs"]
            tau, se = r["Coefficient"], r["Coefficient_std_error"]
            rows.append(dict(trait_id=trait, annot=name, ref_pop=pop, h2g=h2, h2g_se=h2se,
                             h2_z=h2 / h2se if h2se > 0 else np.nan,
                             prop_snps=p, prop_h2=r["Prop._h2"], enrichment=r["Enrichment"],
                             enrichment_se=r["Enrichment_std_error"], enrichment_p=r["Enrichment_p"],
                             tau=tau, tau_se=se, tau_z=r["Coefficient_z-score"],
                             p_tau=norm.sf(r["Coefficient_z-score"]),
                             tau_star=tau * M[pop] * np.sqrt(p * (1 - p)) / h2m if h2m > 0 else np.nan))
    res = pd.DataFrame(rows).merge(man, on="trait_id", how="left")
    res["ld_match"] = ~res.source.isin(["PAGE", "UKB_AFR"])
    res["family"] = res.annot.str.split("_").str[0]
    res["primary"] = res.ld_match & (res.h2_z > 5)
    res["FDR"] = np.nan
    for fam, idx in res[res.primary].groupby("family").groups.items():
        res.loc[idx, "FDR"] = bh(res.loc[idx, "p_tau"])
    res = res.sort_values("p_tau")
    res.to_csv(f"{L}/sldsc_{batch}.tsv", sep="\t", index=False, float_format="%.4g")
    shutil.copy(f"{L}/sldsc_{batch}.tsv", f"{BASE}/GWAS/results/sldsc_{batch}.tsv")
    print(batch, len(res), "rows;", res.trait_id.nunique(), "traits;", int(res.primary.sum()), "primary rows")
    print(res[res.primary].groupby("family").apply(lambda d: pd.Series(
        {"n": len(d), "FDR<0.05": int((d.FDR < 0.05).sum()), "FDR<0.1": int((d.FDR < 0.1).sum())})).to_string())
    cols = ["annot", "source", "label", "h2_z", "enrichment", "tau_star", "tau_z", "p_tau", "FDR"]
    top = res[res.primary].head(40).copy()
    top["label"] = top.label.str[:35]
    print(top[cols].to_string())
