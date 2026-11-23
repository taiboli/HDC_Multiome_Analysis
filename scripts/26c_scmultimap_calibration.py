"""
scMultiMap null calibration (23_scmultimap.py --perm: RNA cells permuted within donor, Neuronal and iPSC).
For a grid of minimum detection rates (gene, peak: fraction of the cell type's cells with count > 0), BH FDR 10% links
in the real vs permuted data. Under the permutation the bulk is calibrated (lambda ~0.95) but sparse features give
a positive heavy tail (co-occurring rare counts; normal approximation fails), so links are restricted to features
above the threshold where permuted / real links is small (chosen: gene >= 10%, peak >= 2%, both >= 100 cells;
used in 26_link_analysis.py).
Output: PeakAnnot/results/scmultimap_calibration.tsv
Usage: python 26c_scmultimap_calibration.py iPSC Neuronal
"""
import os
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp

PA=os.environ["HDC_DIR"] + "/PeakAnnot/"
def bh(p):
    p=np.asarray(p); o=np.argsort(p); m=len(p); q=np.empty(m); q[o]=np.minimum(np.minimum.accumulate((p[o]*m/np.arange(1,m+1))[::-1])[::-1],1); return q
fr=pd.read_csv(f"{PA}/sc/features_RNA.tsv.gz",sep="\t"); fa=pd.read_csv(f"{PA}/sc/features_ATAC.tsv.gz",sep="\t")
rows=[]
for ct in sys.argv[1:]:
    n=len(pd.read_csv(f"{PA}/sc/{ct}/cells.tsv.gz",sep="\t"))
    g=pd.Series(np.diff(sp.load_npz(f"{PA}/sc/{ct}/RNA.npz").indptr)/n,index=fr.id)
    p=pd.Series(np.diff(sp.load_npz(f"{PA}/sc/{ct}/ATAC.npz").indptr)/n,index=fa.id)
    D={}
    for t in [ct, ct+"_perm"]:
        d=pd.read_parquet(f"{PA}/scMultiMap/{t}.pairs.parquet",columns=["gene","peak","pval","test_stat"]).dropna(subset=["pval"])
        d["gd"]=g.reindex(d.gene).to_numpy(); d["pd"]=p.reindex(d.peak).to_numpy(); D[t]=d
    for gt in [0.01,0.05,0.1,0.2]:
        for pt in [0.01,0.02,0.05,0.1]:
            r={"ct":ct,"gene_det":gt,"peak_det":pt}
            for t,d in D.items():
                x=d[(d.gd>=gt)&(d["pd"]>=pt)]
                q=bh(x.pval.to_numpy())
                k="perm" if t.endswith("perm") else "real"
                r[f"{k}_pairs"]=len(x); r[f"{k}_fdr10"]=int((q<0.1).sum()); r[f"{k}_p1e5"]=int((x.pval<1e-5).sum())
            r["perm/real_fdr10"]=r["perm_fdr10"]/max(r["real_fdr10"],1); r["perm/real_p1e5"]=r["perm_p1e5"]/max(r["real_p1e5"],1)
            rows.append(r)
out = pd.DataFrame(rows)
out.to_csv(f"{PA}/results/scmultimap_calibration.tsv", sep="\t", index=False, float_format="%.4g")
print(out.to_string())
