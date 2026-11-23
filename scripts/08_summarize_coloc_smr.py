"""
Collect coloc / SMR outputs of 07_coloc_smr.py across cell types, modalities and chromosomes.

Writes GWAS/results/:
  coloc_all.tsv.gz          every reported feature x trait coloc (GWAS p<1e-6 in window), with trait annotation
  coloc_H4gt0.5.tsv         PP.H4 > 0.5
  coloc_summary.tsv         per cell type x modality x GWAS source: # tests, # PP.H4>0.5 / >0.8, # traits, # loci
  smr_gwsLeads.tsv.gz       SMR at lead QTL variants that are genome-wide significant for the trait (p<5e-8), with the
                            QTL z taken from the permutation p (pval_beta) to avoid lead-variant winner's curse;
                            BH FDR within cell type x modality x trait (FDR_SMR, lenient: few tests per group) and
                            across all SMR tests (FDR_SMR_global, recommended); joined with coloc PP.H4 when available
"""
import glob
import os

import numpy as np
import pandas as pd
from scipy.stats import chi2, norm

BASE = os.environ["HDC_DIR"] + "/"
R = f"{BASE}/GWAS/results/"
os.makedirs(R, exist_ok=True)
man = pd.read_csv(f"{BASE}/GWAS/gwas_manifest.tsv", sep="\t")[["trait_id", "source", "ancestry", "label", "type"]]


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    if n == 0:
        return p
    o = np.argsort(p)
    q = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1)
    return out


def loci(df, gap=500_000):
    """number of distinct loci: per trait, best SNPs clustered when < gap apart"""
    n = 0
    for _, d in df.groupby(["trait_id", "chr"]):
        p = np.sort(d.best_pos38.to_numpy())
        n += 1 + int((np.diff(p) > gap).sum())
    return n


coloc, smr = [], []
for mod in ("RNA", "ATAC"):
    for ct_dir in sorted(glob.glob(f"{BASE}/GWAS/coloc/{mod}/*/")):
        ct = ct_dir.rstrip("/").split("/")[-1]
        for f in sorted(glob.glob(f"{ct_dir}/chr*.coloc.parquet")):
            d = pd.read_parquet(f)
            d.insert(0, "modality", mod)
            d.insert(0, "celltype", ct)
            coloc.append(d)
        lead_files = sorted(glob.glob(f"{ct_dir}/chr*.lead.parquet"))
        if not lead_files:
            continue
        lead = pd.concat([pd.read_parquet(f).assign(chr=int(os.path.basename(f).split(".")[0][3:]))
                          for f in lead_files], ignore_index=True)
        traits = [c for c in lead.columns if "__" in c]
        for t in traits:
            zg = lead[t].to_numpy(float)
            sel = np.isfinite(zg) & (zg * zg > 29.716)
            if not sel.any():
                continue
            d = lead.loc[sel, ["chr", "phenotype_id", "variant_id", "pos38", "af", "pval_beta", "qval", "slope", "z_qtl"]].copy()
            zq, z = d.z_qtl.to_numpy(), zg[sel]
            # instrument strength from the permutation (beta-approximated) p-value: the lead variant's nominal
            # p is the minimum over thousands of cis tests (winner's curse), far too optimistic with n = 17
            zp = norm.isf(np.clip(d.pval_beta.to_numpy(), 1e-300, 1) / 2) * np.sign(d.slope.to_numpy())
            d["trait_id"], d["z_gwas"] = t, z
            d["p_SMR_nominal"] = chi2.sf(z * z * zq * zq / (z * z + zq * zq), 1)  # classic SMR, reference only
            d["p_SMR"] = chi2.sf(z * z * zp * zp / (z * z + zp * zp), 1)
            d["b_SMR_ratio_z"] = z / zq  # sign = direction of trait effect per increase in feature
            d["FDR_SMR"] = bh(d.p_SMR)
            d.insert(0, "modality", mod)
            d.insert(0, "celltype", ct)
            smr.append(d)
        print(mod, ct, flush=True)

coloc = pd.concat(coloc, ignore_index=True).merge(man, on="trait_id", how="left")
coloc.to_csv(f"{R}/coloc_all.tsv.gz", sep="\t", index=False, float_format="%.4g")
hits = coloc[coloc["PP.H4"] > 0.5].sort_values("PP.H4", ascending=False)
hits.to_csv(f"{R}/coloc_H4gt0.5.tsv", sep="\t", index=False, float_format="%.4g")

rows = []
for (ct, mod, src), d in coloc.groupby(["celltype", "modality", "source"]):
    h5, h8 = d[d["PP.H4"] > 0.5], d[d["PP.H4"] > 0.8]
    rows.append({"celltype": ct, "modality": mod, "source": src, "n_tests": len(d),
                 "n_features_tested": d.phenotype_id.nunique(), "n_H4>0.5": len(h5), "n_H4>0.8": len(h8),
                 "features_H4>0.5": h5.phenotype_id.nunique(), "traits_H4>0.5": h5.trait_id.nunique(),
                 "loci_H4>0.5": loci(h5) if len(h5) else 0, "mean_PP.H4": d["PP.H4"].mean(),
                 "mean_PP.H3": d["PP.H3"].mean()})
pd.DataFrame(rows).to_csv(f"{R}/coloc_summary.tsv", sep="\t", index=False, float_format="%.4g")

smr = pd.concat(smr, ignore_index=True).merge(man, on="trait_id", how="left")
smr["FDR_SMR_global"] = bh(smr.p_SMR)
smr = smr.merge(coloc[["celltype", "modality", "trait_id", "phenotype_id", "PP.H3", "PP.H4", "best_snp"]],
                on=["celltype", "modality", "trait_id", "phenotype_id"], how="left")
smr.to_csv(f"{R}/smr_gwsLeads.tsv.gz", sep="\t", index=False, float_format="%.4g")
print(pd.DataFrame(rows).groupby("modality")[["n_tests", "n_H4>0.5", "n_H4>0.8"]].sum())
print("SMR tests:", len(smr), "; FDR_SMR<0.05:", (smr.FDR_SMR < 0.05).sum(), "; FDR_SMR_global<0.05:",
      (smr.FDR_SMR_global < 0.05).sum(), "; FDR_SMR<0.05 & PP.H4>0.5:", ((smr.FDR_SMR < 0.05) & (smr["PP.H4"] > 0.5)).sum())
