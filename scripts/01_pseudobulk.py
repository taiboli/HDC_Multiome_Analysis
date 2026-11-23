"""
Pseudobulk HDC multiome counts per (highres_celltype, donor) in one streaming pass.

For every cell type and modality (RNA = 'Gene Expression', ATAC = 'Peaks') writes to PseudoBulk/:
  {CellType}_{RNA|ATAC}_SumCount.tsv.gz   feature x donor summed raw counts
  {CellType}_{RNA|ATAC}_NCellExpr.tsv.gz  feature x donor number of cells with count > 0
  {CellType}_NCell.tsv                    number of cells per donor
plus features_RNA.tsv.gz / features_ATAC.tsv.gz (id, chr, start, end).

These are exactly the quantities used by the EB_CellTypeData pipeline
(rowSums of counts per donor, and cells.count / cells.perc filters), so the
downstream R step can reproduce the same filtering + TMM + invnorm logic.

Only labeled 'paired' multiome cells with a donor assignment are used.
"""
import os
import re
import sys
import time

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

H5AD = os.environ["HDC_DIR"] + "/MultiVI/HDC_multiome_RemoveHighRibo_posttraining_adata_germlayer.h5ad"
OUT = os.environ["HDC_DIR"] + "/PseudoBulk/"
CHUNK = 4000


def read_col(group, key):
    x = group[key]
    if isinstance(x, h5py.Group):  # categorical
        cats = x["categories"][:].astype(str)
        codes = x["codes"][:]
        out = np.where(codes >= 0, cats[np.clip(codes, 0, None)], "")
        return out.astype(object)
    v = x[:]
    return v.astype(str) if v.dtype.kind in "SO" else v


def clean(name):
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")


t0 = time.time()
f = h5py.File(H5AD, "r")
obs, var = f["obs"], f["var"]

celltype = read_col(obs, "highres_celltype")
donor = read_col(obs, "donor")
modality = read_col(obs, "modality")
keep = (celltype != "") & (donor != "") & (modality == "paired")
print(f"cells kept: {keep.sum()} / {len(keep)}", flush=True)

groups = pd.DataFrame({"celltype": celltype[keep], "donor": donor[keep]})
groups["g"] = groups.celltype + "|" + groups.donor
glevels = np.array(sorted(groups.g.unique()))
gcode = np.full(len(keep), -1, dtype=np.int64)
gcode[np.where(keep)[0]] = pd.Categorical(groups.g, categories=glevels).codes
print(f"{len(glevels)} celltype x donor groups", flush=True)

# features
fid = read_col(var, "_index")
fmod = read_col(var, "modality")
fchr = read_col(var, "chr")
fstart = var["start"][:]
fend = var["end"][:]
nfeat = len(fid)

# accumulators: groups x features
sum_mat = np.zeros((len(glevels), nfeat), dtype=np.float64)
nz_mat = np.zeros((len(glevels), nfeat), dtype=np.float64)

X = f["X"]
indptr = X["indptr"][:].astype(np.int64)
ncell = len(indptr) - 1
for s in range(0, ncell, CHUNK):
    e = min(s + CHUNK, ncell)
    gc = gcode[s:e]
    if (gc >= 0).sum() == 0:
        continue
    a, b = indptr[s], indptr[e]
    data = X["data"][a:b].astype(np.float64)
    idx = X["indices"][a:b]
    m = sp.csr_matrix((data, idx, indptr[s:e + 1] - a), shape=(e - s, nfeat))
    rows = np.where(gc >= 0)[0]
    G = sp.csr_matrix((np.ones(len(rows)), (gc[rows], rows)), shape=(len(glevels), e - s))
    sum_mat += (G @ m).toarray()
    m.data[:] = 1.0
    nz_mat += (G @ m).toarray()
    if (s // CHUNK) % 5 == 0:
        print(f"  {e}/{ncell} cells, {time.time() - t0:.0f}s", flush=True)

ncells = groups.g.value_counts()

for mod, label in [("Gene Expression", "RNA"), ("Peaks", "ATAC")]:
    fsel = np.where(fmod == mod)[0]
    pd.DataFrame({"id": fid[fsel], "chr": fchr[fsel], "start": fstart[fsel], "end": fend[fsel]}).to_csv(
        f"{OUT}/features_{label}.tsv.gz", sep="\t", index=False)

for ct in sorted(groups.celltype.unique()):
    gi = [i for i, g in enumerate(glevels) if g.split("|")[0] == ct]
    donors = [glevels[i].split("|")[1] for i in gi]
    ctn = clean(ct)
    pd.DataFrame({"donor": donors, "ncell": [ncells[glevels[i]] for i in gi]}).to_csv(
        f"{OUT}/{ctn}_NCell.tsv", sep="\t", index=False)
    for mod, label in [("Gene Expression", "RNA"), ("Peaks", "ATAC")]:
        fsel = np.where(fmod == mod)[0]
        for mat, suffix in [(sum_mat, "SumCount"), (nz_mat, "NCellExpr")]:
            df = pd.DataFrame(mat[np.ix_(gi, fsel)].T.astype(np.int64), columns=donors)
            df.insert(0, "id", fid[fsel])
            df.to_csv(f"{OUT}/{ctn}_{label}_{suffix}.tsv.gz", sep="\t", index=False)
    print(f"wrote {ctn} ({len(donors)} donors), {time.time() - t0:.0f}s", flush=True)

print("done", flush=True)
