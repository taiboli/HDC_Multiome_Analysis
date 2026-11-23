"""
Check the vectorized coloc in 07_coloc_smr.py against coloc::coloc.abf (R) for a sample of feature x trait pairs.
Usage: python 07v_validate_coloc.py <CellType> <RNA|ATAC> <chr> <outdir>   then   Rscript 07v_validate_coloc.R <outdir>
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm

BASE = os.environ["HDC_DIR"] + "/"
TAG = "noBootstrap_TMM_invnorm"
ct, mod, chrom, od = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
os.makedirs(od, exist_ok=True)
qdir = f"{BASE}/{'eQTL' if mod == 'RNA' else 'caQTL'}_Output_optCov/{ct}/"
res = pd.read_parquet(f"{BASE}/GWAS/coloc/{mod}/{ct}/chr{chrom}.coloc.parquet")
res = res.sort_values("PP.H4")
pick = res.iloc[np.unique(np.linspace(0, len(res) - 1, 10).astype(int))]
man = pd.read_csv(f"{BASE}/GWAS/gwas_manifest.tsv", sep="\t").set_index("trait_id")
u = pd.read_parquet(f"{BASE}/GWAS/variant_universe.parquet", columns=["vidx", "id"], filters=[("chr", "==", chrom)])
nom = pd.read_parquet(f"{qdir}/chr{chrom}_cis_nominal_{TAG}.cis_qtl_pairs.{chrom}.parquet",
                      columns=["phenotype_id", "variant_id", "pval_nominal", "slope"],
                      filters=[("phenotype_id", "in", list(pick.phenotype_id.unique()))])
nom = nom.merge(u, left_on="variant_id", right_on="id")
for i, (_, r) in enumerate(pick.iterrows()):
    g = pd.read_parquet(f"{BASE}/GWAS/harmonized/{r.trait_id}.parquet", filters=[("chr", "==", chrom)])
    qc = pd.read_csv(f"{BASE}/GWAS/harmonized/qc/{r.trait_id}.tsv", sep="\t").iloc[0]
    d = nom[nom.phenotype_id == r.phenotype_id].merge(g, on="vidx")
    d = d[np.isfinite(d.slope) & (d.slope != 0)]
    z = norm.isf(np.clip(d.pval_nominal, 1e-300, 1) / 2) * np.sign(d.slope)
    d["varbeta_qtl"] = (d.slope / z) ** 2
    d["varbeta_gwas"] = d.se.astype(float) ** 2
    d[["variant_id", "slope", "varbeta_qtl", "beta", "varbeta_gwas"]].to_csv(f"{od}/pair{i}.tsv", sep="\t", index=False)
    sdy = qc.sdY if man.loc[r.trait_id, "type"] == "quant" and np.isfinite(qc.sdY) else 1.0
    pd.DataFrame([dict(pair=i, trait=r.trait_id, feature=r.phenotype_id, type=man.loc[r.trait_id, "type"], sdY=sdy,
                       nsnps=r.nsnps, H3=r["PP.H3"], H4=r["PP.H4"])]).to_csv(
        f"{od}/meta.tsv", sep="\t", index=False, mode="a" if i else "w", header=not i)
print(pick[["trait_id", "phenotype_id", "nsnps", "PP.H3", "PP.H4"]])
