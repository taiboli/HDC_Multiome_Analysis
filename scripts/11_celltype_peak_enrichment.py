"""
Donor-number-independent test of cell-type relevance: GWAS signal in cell-type-specific HDC chromatin.
Uses all cells (pseudobulk over all 17 donors), so it does not depend on QTL power; includes Ectoderm_2.

Cell-type-specific peaks: from PseudoBulk/<CT>_ATAC_SumCount.tsv.gz (all 524k peaks), CPM per cell type
(summed over donors); specificity = CPM_ct / mean CPM of the other cell types; per cell type the top 10% of
peaks by specificity among peaks with CPM_ct above the median CPM_ct.
For each GWAS trait: mean chi2 of harmonized variants inside each cell type's specific peaks divided by the
mean chi2 of variants inside any HDC peak (background), SE by delete-one-block jackknife over 200 contiguous
genomic blocks of variants; also a z-score of the log-ratio across cell types (relative ranking).
Exploratory (no LD-score modelling); S-LDSC is the follow-up.

Output: GWAS/results/celltype_peak_enrichment.tsv, GWAS/celltype_specific_peaks.tsv.gz
"""
import glob
import os

import numpy as np
import pandas as pd
from scipy.stats import norm

BASE = os.environ["HDC_DIR"] + "/"
R = f"{BASE}/GWAS/results/"
NBLOCK = 200

cts = sorted(f.split("/")[-1].replace("_ATAC_SumCount.tsv.gz", "")
             for f in glob.glob(f"{BASE}/PseudoBulk/*_ATAC_SumCount.tsv.gz")
             if not f.split("/")[-1].startswith("AllCells"))  # pooled pseudobulk (01b) is not a cell type
cpm = {}
for ct in cts:
    x = pd.read_csv(f"{BASE}/PseudoBulk/{ct}_ATAC_SumCount.tsv.gz", sep="\t", index_col=0).sum(axis=1)
    cpm[ct] = x / x.sum() * 1e6
cpm = pd.DataFrame(cpm).fillna(0)
spec = {}
for ct in cts:
    other = cpm.drop(columns=ct).mean(axis=1)
    s = (cpm[ct] + 0.1) / (other + 0.1)
    s = s[cpm[ct] > cpm[ct].median()]
    spec[ct] = set(s.sort_values(ascending=False).index[: int(0.1 * len(cpm))])
peaks = pd.DataFrame({"peak": cpm.index})
p = peaks.peak.str.split("_", expand=True)
peaks["chr"] = pd.to_numeric(p[0].str.replace("chr", ""), errors="coerce")
peaks["start"], peaks["end"] = p[1].astype(int), p[2].astype(int)
for ct in cts:
    peaks[ct] = peaks.peak.isin(spec[ct])
peaks = peaks.dropna(subset=["chr"])
peaks["chr"] = peaks.chr.astype(int)
peaks.to_csv(f"{BASE}/GWAS/celltype_specific_peaks.tsv.gz", sep="\t", index=False)

# map universe variants to peaks (peaks are non-overlapping after sorting; take containing peak)
u = pd.read_parquet(f"{BASE}/GWAS/variant_universe.parquet", columns=["vidx", "chr", "pos38"])
in_peak = np.full(len(u), -1)
for c, pk in peaks.groupby("chr"):
    pk = pk.sort_values("start")
    uu = u[u.chr == c]
    j = np.searchsorted(pk.start.to_numpy(), uu.pos38.to_numpy(), side="right") - 1
    ok = (j >= 0) & (uu.pos38.to_numpy() <= pk.end.to_numpy()[np.maximum(j, 0)])
    in_peak[uu.index[ok]] = pk.index.to_numpy()[j[ok]]
u["peak_row"] = in_peak
u = u[u.peak_row >= 0].reset_index(drop=True)
mask = {ct: peaks[ct].to_numpy()[peaks.index.get_indexer(u.peak_row)] for ct in cts}
print(f"{len(u)} variants in HDC peaks", flush=True)

man = pd.read_csv(f"{BASE}/GWAS/gwas_manifest.tsv", sep="\t")
rows = []
for t in man.trait_id:
    if not os.path.exists(f"{BASE}/GWAS/harmonized/{t}.parquet"):
        continue
    g = pd.read_parquet(f"{BASE}/GWAS/harmonized/{t}.parquet", columns=["vidx", "beta", "se"])
    gi = g.set_index("vidx")
    z = (gi.beta / gi.se).reindex(u.vidx).to_numpy(float)
    ok = np.isfinite(z)
    chi = z[ok] ** 2
    blk = np.arange(ok.sum()) * NBLOCK // ok.sum()
    bg_sum, bg_n = np.bincount(blk, chi, NBLOCK), np.bincount(blk, None, NBLOCK)
    for ct in cts:
        m = mask[ct][ok]
        s_sum, s_n = np.bincount(blk, chi * m, NBLOCK), np.bincount(blk, m.astype(float), NBLOCK)
        est = np.log((s_sum.sum() / s_n.sum()) / (bg_sum.sum() / bg_n.sum()))
        jk = np.log(((s_sum.sum() - s_sum) / (s_n.sum() - s_n)) / ((bg_sum.sum() - bg_sum) / (bg_n.sum() - bg_n)))
        se = np.sqrt((NBLOCK - 1) / NBLOCK * np.sum((jk - jk.mean()) ** 2))
        rows.append(dict(trait_id=t, celltype=ct, n_var_specific=int(s_n.sum()), enrichment=np.exp(est),
                         log_enr=est, se_jk=se, p=2 * norm.sf(abs(est / se))))
    print(t, flush=True)
res = pd.DataFrame(rows).merge(man[["trait_id", "source", "ancestry", "label"]], on="trait_id")
res["z_across_celltypes"] = res.groupby("trait_id").log_enr.transform(lambda x: (x - x.mean()) / x.std())
pv = res.p.to_numpy()
o = np.argsort(pv)
q = np.minimum.accumulate((pv[o] * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1]
res["FDR"] = np.nan
res.loc[res.index[o], "FDR"] = np.minimum(q, 1)
res.to_csv(f"{R}/celltype_peak_enrichment.tsv", sep="\t", index=False, float_format="%.4g")
print(res.sort_values("p").head(30).to_string())
