"""
Per-cell-type single-cell count matrices for peak-gene linking (scMultiMap, 23_scmultimap.py) and peak
characterization, in one streaming pass over the HDC multiome h5ad (same cells as 01_pseudobulk.py:
labeled 'paired' cells with a donor; cell type = highres_celltype).

Output PeakAnnot/sc/<CT>/:
  RNA.npz, ATAC.npz   scipy CSC (cells x features), float32 counts; RNA = all 'Gene Expression' features,
                      ATAC = all 'Peaks' features (column order = features_{RNA,ATAC}.tsv.gz)
  cells.tsv.gz        barcode, donor, depth_RNA, depth_ATAC (total counts over all features of the modality)
PeakAnnot/sc/features_{RNA,ATAC}.tsv.gz   id, gene_id (RNA), chr, start, end
"""
import os
import re
import time

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

H5AD = os.environ["HDC_DIR"] + "/MultiVI/HDC_multiome_RemoveHighRibo_posttraining_adata_germlayer.h5ad"
OUT = os.environ["HDC_DIR"] + "/PeakAnnot/sc/"
CHUNK = 5000


def read_col(group, key):
    x = group[key]
    if isinstance(x, h5py.Group):  # categorical
        cats = x["categories"][:].astype(str)
        codes = x["codes"][:]
        return np.where(codes >= 0, cats[np.clip(codes, 0, None)], "").astype(object)
    v = x[:]
    return v.astype(str) if v.dtype.kind in "SO" else v


def clean(name):
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")


t0 = time.time()
f = h5py.File(H5AD, "r")
obs, var = f["obs"], f["var"]
celltype = np.array([clean(c) if c else "" for c in read_col(obs, "highres_celltype")], dtype=object)
donor = read_col(obs, "donor")
modality = read_col(obs, "modality")
barcode = read_col(obs, obs.attrs["_index"])
keep = (celltype != "") & (donor != "") & (modality == "paired")
cts = sorted(set(celltype[keep]))
print(f"cells kept {keep.sum()}; cell types {cts}", flush=True)

fmod = read_col(var, "modality")
is_rna = fmod == "Gene Expression"
is_atac = fmod == "Peaks"
feat = pd.DataFrame({"id": read_col(var, "_index"), "gene_id": read_col(var, "ID"), "chr": read_col(var, "chr"),
                     "start": var["start"][:], "end": var["end"][:]})
os.makedirs(OUT, exist_ok=True)
feat[is_rna].to_csv(f"{OUT}/features_RNA.tsv.gz", sep="\t", index=False)
feat[is_atac].drop(columns="gene_id").to_csv(f"{OUT}/features_ATAC.tsv.gz", sep="\t", index=False)
irna, iatac = np.where(is_rna)[0], np.where(is_atac)[0]

X = f["X"]
indptr = X["indptr"][:].astype(np.int64)
ncell, nfeat = len(indptr) - 1, len(fmod)
parts = {ct: {"RNA": [], "ATAC": [], "cells": []} for ct in cts}
for s in range(0, ncell, CHUNK):
    e = min(s + CHUNK, ncell)
    k = keep[s:e]
    if not k.any():
        continue
    a, b = indptr[s], indptr[e]
    m = sp.csr_matrix((X["data"][a:b].astype(np.float32), X["indices"][a:b], indptr[s:e + 1] - a), shape=(e - s, nfeat))
    mr, ma = m[:, irna], m[:, iatac]
    dr, da = np.asarray(mr.sum(1)).ravel(), np.asarray(ma.sum(1)).ravel()
    cts_chunk = celltype[s:e]
    for ct in set(cts_chunk[k]):
        sel = np.where(k & (cts_chunk == ct))[0]
        parts[ct]["RNA"].append(mr[sel])
        parts[ct]["ATAC"].append(ma[sel])
        parts[ct]["cells"].append(pd.DataFrame({"barcode": barcode[s:e][sel], "donor": donor[s:e][sel],
                                                "depth_RNA": dr[sel], "depth_ATAC": da[sel]}))
    print(f"{e}/{ncell} cells, {time.time() - t0:.0f}s", flush=True)

for ct in cts:
    d = f"{OUT}/{ct}/"
    os.makedirs(d, exist_ok=True)
    for mod in ("RNA", "ATAC"):
        sp.save_npz(f"{d}/{mod}.npz", sp.vstack(parts[ct][mod]).tocsc(), compressed=False)
        parts[ct][mod] = None
    pd.concat(parts[ct]["cells"]).to_csv(f"{d}/cells.tsv.gz", sep="\t", index=False)
    print(f"{ct}: {sum(len(x) for x in parts[ct]['cells'])} cells written, {time.time() - t0:.0f}s", flush=True)
