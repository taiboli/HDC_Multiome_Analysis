"""
Extract per-cell-type QTL z-scores for mashr (one chromosome, RNA or ATAC), from the optCov tensorQTL output.

strong set: for each feature (gene / peak), the variant with the smallest nominal p among the per-cell-type
            permutation leads ("top cross-cell-type variant"); z of that pair in every cell type
random set: deterministic pseudo-random pairs, identical in every cell type: features whose id hash
            falls in a FEAT_FRAC slice x variants whose rsID number % VAR_MOD == 0 (unbiased w.r.t. QTL status)
z = sign(slope) * Phi^-1(1 - p_nominal / 2)   (t -> z through the p-value, as in 07_coloc_smr.py)
Missing (feature not tested in a cell type) -> absent row; 14_mash.R fills Bhat = 0, Shat = 1e3.

Output: mash/extract/<MOD>_chr<N>_{strong,random}.parquet  (phenotype_id, variant_id, celltype, z)
Usage: python 13_mash_extract.py <RNA|ATAC> <chr>
"""
import os
import sys
import zlib

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
from scipy.stats import norm

BASE = os.environ["HDC_DIR"] + "/"
TAG = "noBootstrap_TMM_invnorm"
CTS = ["Ectoderm_1", "Ectoderm_3", "Endoderm", "iPSC", "Mesoderm", "Neural_crest", "Neuronal", "Retinal",
       "Unknown", "Vascular_endothelial"]
mod, chrom = sys.argv[1], sys.argv[2]
QD = "eQTL" if mod == "RNA" else "caQTL"
FEAT_FRAC, VAR_MOD = (0.25, 50) if mod == "RNA" else (0.02, 50)
OUT = f"{BASE}/mash/extract/"
os.makedirs(OUT, exist_ok=True)

# ---- strong pairs: top variant across cell-type permutation leads
leads = []
for ct in CTS:
    f = f"{BASE}/{QD}_Output_optCov/{ct}/chr{chrom}_cis_{TAG}.cis_qtl.txt.gz"
    if os.path.exists(f):
        leads.append(pd.read_csv(f, sep="\t", usecols=["phenotype_id", "variant_id", "pval_nominal"]))
leads = pd.concat(leads).sort_values("pval_nominal").drop_duplicates("phenotype_id")
strong = set(zip(leads.phenotype_id, leads.variant_id))
strong_feats = pa.array(sorted(leads.phenotype_id.unique()))

# ---- random features (stable across cell types)
all_feats = sorted(leads.phenotype_id.unique())
rand_feats = pa.array([f for f in all_feats if (zlib.crc32(f.encode()) % 10_000) < FEAT_FRAC * 10_000])

out_s, out_r = [], []
for ct in CTS:
    f = f"{BASE}/{QD}_Output_optCov/{ct}/chr{chrom}_cis_nominal_{TAG}.cis_qtl_pairs.{chrom}.parquet"
    if not os.path.exists(f):
        continue
    pf = pq.ParquetFile(f)
    for i in range(pf.num_row_groups):
        t = pf.read_row_group(i, columns=["phenotype_id", "variant_id", "pval_nominal", "slope"])
        ts = t.filter(pc.is_in(t["phenotype_id"], value_set=strong_feats))
        if ts.num_rows:
            d = ts.to_pandas()
            k = pd.Series(list(zip(d.phenotype_id, d.variant_id))).isin(strong).to_numpy()
            if k.any():
                out_s.append(d[k].assign(celltype=ct))
        tr = t.filter(pc.is_in(t["phenotype_id"], value_set=rand_feats))
        if tr.num_rows:
            d = tr.to_pandas()
            rs = pd.to_numeric(d.variant_id.str.replace("rs", "", regex=False), errors="coerce")
            out_r.append(d[(rs % VAR_MOD == 0).to_numpy()].assign(celltype=ct))
    print(ct, flush=True)


def finish(lst, name):
    d = pd.concat(lst, ignore_index=True)
    d = d[d.slope != 0]
    d["z"] = norm.isf(np.clip(d.pval_nominal, 1e-300, 1) / 2) * np.sign(d.slope)
    d[["phenotype_id", "variant_id", "celltype", "z"]].to_parquet(f"{OUT}/{mod}_chr{chrom}_{name}.parquet", index=False)
    print(name, len(d), d.groupby(["phenotype_id", "variant_id"]).ngroups, "pairs", flush=True)


finish(out_s, "strong")
finish(out_r, "random")
print("done")
