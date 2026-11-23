"""
S-LDSC annotation files + LD scores for one batch x reference population x chromosome.

  make : binary annotation per set of 15_ldsc_sets.py over the 1000G Phase 3 <pop> plink SNPs
         -> GWAS/ldsc/annot/<batch>/<pop>/all.<chr>.annot.gz (CHR BP SNP CM + one column per set)
  (then ldsc.py --l2 --ld-wind-cm 1, printing the SNPs of the matching baselineLD v2.2 LD-score file so the
   regression can concatenate them with the baseline; see submit_15b_ldsc_l2.sh)
  split: per-set .annot.gz / .l2.ldscore.gz / .l2.M / .l2.M_5_50 under GWAS/ldsc/ldscores/<pop>/<set>/<set>.<chr>.*
         so each set can be tested on its own on top of baselineLD (+ the batch's control set)

Reference population: EUR for FinnGen / UKB EUR / UKB meta; EAS for BBJ. No AFR reference exists in the
Zenodo release (record 10515792), so AFR-rich GWAS (PAGE, UKB AFR) are run with EUR LD (LD mismatch, flagged).
The annotations themselves are genomic intervals (peaks, gene windows) and do not use the YRI donors' LD.
Usage: python 15b_ldsc_annot.py <make|split> <batch> <EUR|EAS> <chr>
"""
import glob
import os
import sys

import numpy as np
import pandas as pd

BASE = os.environ["HDC_DIR"] + "/"
REF = os.environ["LDSC_REF_DIR"] + "/"
L = f"{BASE}/GWAS/ldsc/"
mode, batch, pop, chrom = sys.argv[1:5]
BFILE = {"EUR": "1000G_Phase3_plinkfiles/1000G_EUR_Phase3_plink/1000G.EUR.QC.{}",
         "EAS": "1000G_Phase3_EAS_plinkfiles/*/1000G.EAS.QC.{}"}[pop]
bim = glob.glob(REF + BFILE.format(chrom) + ".bim")[0]
ad = f"{L}/annot/{batch}/{pop}/"
os.makedirs(ad, exist_ok=True)

if mode == "make":
    b = pd.read_csv(bim, sep="\t", header=None, names=["CHR", "SNP", "CM", "BP", "A1", "A2"])
    s = pd.read_csv(f"{L}/sets_{batch}.tsv.gz", sep="\t", dtype={"chr": str})
    s = s[s.chr == str(chrom)]
    out = b[["CHR", "BP", "SNP", "CM"]].copy()
    bp = b.BP.to_numpy()
    for name, d in s.groupby("annot", sort=False):
        # merge intervals, then SNP in [start, end]
        d = d.sort_values("start")
        st, en = d.start.to_numpy(), d.end.to_numpy()
        ends = np.maximum.accumulate(en)
        new = np.r_[True, st[1:] > ends[:-1]]
        ms, me = st[new], np.r_[ends[np.r_[new[1:], True]]]
        i = np.searchsorted(ms, bp, side="right") - 1
        out[name] = ((i >= 0) & (bp <= me[np.clip(i, 0, None)])).astype(np.int8)
    for name in pd.read_csv(f"{L}/sets_{batch}.tsv.gz", sep="\t", usecols=["annot"]).annot.unique():
        if name not in out:
            out[name] = np.int8(0)
    out.to_csv(f"{ad}/all.{chrom}.annot.gz", sep="\t", index=False)
    print(out.iloc[:, 4:].mean().round(4).to_string())

elif mode == "split":
    a = pd.read_csv(f"{ad}/all.{chrom}.annot.gz", sep="\t")
    l2 = pd.read_csv(f"{ad}/all.{chrom}.l2.ldscore.gz", sep="\t")
    M = open(f"{ad}/all.{chrom}.l2.M").read().split()
    M5 = open(f"{ad}/all.{chrom}.l2.M_5_50").read().split()
    names = list(a.columns[4:])
    assert len(M) == len(names) == l2.shape[1] - 3
    for k, n in enumerate(names):
        od = f"{L}/ldscores/{pop}/{n}/"
        os.makedirs(od, exist_ok=True)
        a[["CHR", "BP", "SNP", "CM", n]].to_csv(f"{od}/{n}.{chrom}.annot.gz", sep="\t", index=False)
        l2.iloc[:, [0, 1, 2, 3 + k]].rename(columns={l2.columns[3 + k]: f"{n}L2"}).to_csv(
            f"{od}/{n}.{chrom}.l2.ldscore.gz", sep="\t", index=False, float_format="%.4g")
        open(f"{od}/{n}.{chrom}.l2.M", "w").write(M[k] + "\n")
        open(f"{od}/{n}.{chrom}.l2.M_5_50", "w").write(M5[k] + "\n")
    print("split", len(names))
