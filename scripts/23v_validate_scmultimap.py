"""
Validate scmultimap_py.py against the original R functions (scMultiMap v1.0.1 R/scMultiMap_IRLS.R, R/scMultiMap_WLS*.R,
copied to PeakAnnot/scMultiMap/Rsrc/) on a real subset: Retinal cells, chr21 genes / peaks detected in >= 5% of cells,
bsample = donor. Writes the subset for 23v_validate_scmultimap.R, then compares after the R run (argument 'compare').
"""
import os
import sys

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scmultimap_py import irls, wls_block  # noqa: E402

BASE = os.environ["HDC_DIR"] + "/"
SC, V = f"{BASE}/PeakAnnot/sc/", f"{BASE}/PeakAnnot/scMultiMap/validate/"
ct = "Retinal"
cells = pd.read_csv(f"{SC}/{ct}/cells.tsv.gz", sep="\t")
R = sp.load_npz(f"{SC}/{ct}/RNA.npz").tocsc()
A = sp.load_npz(f"{SC}/{ct}/ATAC.npz").tocsc()
fr = pd.read_csv(f"{SC}/features_RNA.tsv.gz", sep="\t")
fa = pd.read_csv(f"{SC}/features_ATAC.tsv.gz", sep="\t")
n = len(cells)
gi = np.where((fr.chr == "chr21").to_numpy() & (np.diff(R.indptr) / n >= 0.05))[0][:80]
pi = np.where((fa.chr == "chr21").to_numpy() & (np.diff(A.indptr) / n >= 0.05))[0][:1500]
donors = pd.Index(sorted(cells.donor.unique()))
b = donors.get_indexer(cells.donor)
sd_g, sd_p = cells.depth_RNA.to_numpy(float), cells.depth_ATAC.to_numpy(float)
if sys.argv[1:] == ["write"]:
    sio.mmwrite(f"{V}/gene.mtx", R[:, gi])
    sio.mmwrite(f"{V}/peak.mtx", A[:, pi])
    pd.DataFrame({"donor": cells.donor, "sd_g": sd_g, "sd_p": sd_p}).to_csv(f"{V}/cells.tsv", sep="\t", index=False)
    pd.Series(fr.id[gi].to_numpy()).to_csv(f"{V}/genes.txt", index=False, header=False)
    pd.Series(fa.id[pi].to_numpy()).to_csv(f"{V}/peaks.txt", index=False, header=False)
    sys.exit()
pg = irls(R[:, gi].tocsc(), sd_g, b, len(donors), chunk=37)
pp = irls(A[:, pi].tocsc(), sd_p, b, len(donors), chunk=400)
est, ts, pv = wls_block(R[:, gi].toarray().astype(float), sd_g, pg, A[:, pi].toarray().astype(float), sd_p, pp, b)
py = pd.DataFrame({"gene": np.repeat(fr.id[gi].to_numpy(), len(pi)), "peak": np.tile(fa.id[pi].to_numpy(), len(gi)),
                   "pval_py": pv.ravel(), "ts_py": ts.ravel(), "covar_py": est.ravel()})
r = pd.read_csv(f"{V}/R_wls.tsv", sep="\t")
m = r.merge(py, on=["gene", "peak"])
rel = lambda a, c: np.nanmax(np.abs(a - c) / np.maximum(np.abs(a), 1e-12))
print(f"{len(m)} pairs; max rel diff test_stat {rel(m.test_stat, m.ts_py):.2e}, covar {rel(m.covar, m.covar_py):.2e}, "
      f"pval {rel(m.pval, m.pval_py):.2e}; corr(log p) {np.corrcoef(np.log10(m.pval), np.log10(m.pval_py))[0, 1]:.8f}")
rg = pd.read_csv(f"{V}/R_irls_gene.tsv", sep="\t")
print("IRLS gene sigma_sq max rel diff", rel(rg.sigma_sq.to_numpy(), pg["sigma_sq"]), "theta R", open(f"{V}/R_theta.txt").read().split()[:3], "py", pg["theta"][:3])
m.to_csv(f"{V}/compare.tsv.gz", sep="\t", index=False)
