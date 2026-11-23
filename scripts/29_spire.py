"""
SPIRE (sparse OMP + bootstrap stability selection; $SPIRE_DIR, numpy backend) TF -> peak networks
on HDC pseudocells. Same cells as 01 / 22 (PeakAnnot/sc/<CT>/{RNA,ATAC}.npz); TFs = Lambert et al. 2018 human TFs
(TFRegulation/OMP/human_tf.txt, matched by Ensembl id).

Designs (pseudocells are always drawn within one stratum, so none mixes strata):
  ct_<CT>        within-state: one stratum = all cells of the cell type; 300 pseudocells x 50 cells
                 (Vascular_endothelial: 25 cells). Network = covariation inside the cell type.
  nn_<CT>        within-state, neighbourhood pseudocells: 300 random seed cells of the cell type, each pseudocell = the
                 seed + its 49 nearest neighbours (same cell type) in the MultiVI latent space (TFNet/multivi_latent.tsv.gz),
                 so pseudocells keep the cell type's internal heterogeneity (random pooling averages it away).
                 Unpaired control: RNA from 50 other cells among the seed's 150 nearest neighbours.
  all            joint: strata = the 11 cell types, 60 pseudocells x 50 cells each (25 for Vascular_endothelial);
                 network includes between-cell-type variation (as in the SPIRE brain analysis) -> cell-type programs.
  traj_<lin>     trajectory: slingshot lineage ($PSEUDOTIME_FILE, highres_slingshot_pseudotimes.tsv),
                 cells = iPSC + the lineage's cell types with a pseudotime on it (Unknown/stress cells excluded);
                 strata = pseudotime bins (2 equal-count bins of iPSC + 8 of the derived cells), 40 pseudocells x 30
                 cells per bin; network = covariation along differentiation.
  <run>_s2      independent replicate (different pseudocells and bootstrap seed)
  <design>_unpaired  control: RNA of each pseudocell comes from an independent draw of cells of the same stratum
                 (= separately profiled RNA and ATAC, matched by cell state only).
Feature filter: detected in > 1% of the cells of at least one stratum and > 5 counts in total (SPIRE default).
Normalization TMM + log-CPM of pseudocells per modality; OMP (spire_omp_fast.omp_fast = spire.omp.omp) L = 10, eps = 0 on centred, unit-norm TF columns and centred
accessibility; B = 50 replicates, each with a freshly drawn pseudocell set (SPIRE resample = "cells": a fixed pseudocell
set resampled by rows keeps its own noise correlations, which then look stable); consensus = selection frequency >= 0.5.

Output TFNet/runs/<run>/:
  edges.tsv.gz        consensus edges: tf, peak, coef (mean over bootstraps where selected), frequency
  scope.tsv           tf, n_peaks, fraction, mean / max stratum expression (log-CPM)
  freq.npz            sparse TF x peak selection frequency (>= 0.1), tf / peak names
  strata.npz          per-stratum mean log-CPM of TFs and peaks, detection fraction of peaks, stratum labels,
                      pseudotime medians and cell-type composition (traj)
  pseudocells_tf.npz  normalized pseudocell TF matrix + stratum label (for TF-TF collinearity)
  log.txt
Usage: python 29_spire.py <run> [B]
"""
import os
import sys
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, os.environ["SPIRE_SRC"])
import spire.consensus  # noqa: E402
from spire.consensus import bootstrap_omp  # noqa: E402
from spire.normalize import tmm_log_cpm  # noqa: E402
from spire.scope import regulatory_scope  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from spire_omp_fast import omp_fast  # noqa: E402

spire.consensus.omp = omp_fast  # identical results to spire.omp.omp, ~3x faster with ~1,300 TFs (spire_omp_fast.py)

BASE = os.environ["HDC_DIR"] + "/"
SC = f"{BASE}/PeakAnnot/sc/"
OUT = os.environ.get("SPIRE_OUT", f"{BASE}/TFNet/runs/")
SLING = os.environ["PSEUDOTIME_FILE"]
CTS = ["iPSC", "Ectoderm_1", "Ectoderm_2", "Ectoderm_3", "Neural_crest", "Neuronal", "Retinal", "Mesoderm", "Endoderm",
       "Vascular_endothelial", "Unknown"]
LINEAGE = {"neur": ["Neuronal"], "neur_crest": ["Neural_crest"], "ret": ["Ectoderm_1", "Ectoderm_2", "Retinal"],
           "meso": ["Ectoderm_3", "Mesoderm"], "endo": ["Endoderm"]}
L, B, TAU = 10, int(sys.argv[2]) if len(sys.argv) > 2 else 50, 0.5
run = sys.argv[1]
rep2 = run.endswith("_s2")          # independent replicate (new pseudocells + bootstraps) for the reproducibility ceiling
unpaired = "_unpaired" in run
design = run.replace("_s2", "").replace("_unpaired", "")
od = f"{OUT}/{run}/"
os.makedirs(od, exist_ok=True)
log = open(f"{od}/log.txt", "w")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    log.write(s + "\n")
    log.flush()


t0 = time.time()
SEED = 2027 if rep2 else 2026
rng = np.random.default_rng(SEED)
frna = pd.read_csv(f"{SC}/features_RNA.tsv.gz", sep="\t")
fatac = pd.read_csv(f"{SC}/features_ATAC.tsv.gz", sep="\t")
tfl = pd.read_csv(os.environ["TF_LIST"], sep="\t")
tf_idx = np.where(frna.gene_id.isin(tfl.Ensembl) | frna.id.isin(tfl.Symbol))[0]
tf_names = frna.id.to_numpy()[tf_idx]

# ---- cells and strata
cells = []
for ct in CTS:
    d = pd.read_csv(f"{SC}/{ct}/cells.tsv.gz", sep="\t")
    cells.append(d.assign(ct=ct, row=np.arange(len(d))))
cells = pd.concat(cells, ignore_index=True)
k_cells = 30 if design.startswith("traj_") else 50
if design.startswith("ct_") or design.startswith("nn_"):
    ct = design[3:]
    cells = cells[cells.ct == ct].assign(stratum=ct)
    if design.startswith("nn_"):
        lat = pd.read_csv(f"{BASE}/TFNet/multivi_latent.tsv.gz", sep="\t", index_col=0)
        Z = lat.reindex(cells.barcode)[[f"z{i}" for i in range(11)]].to_numpy()
        assert np.isfinite(Z).all()
    n_per = {ct: 300}
    k_cells = 25 if ct == "Vascular_endothelial" else 50
elif design == "all":
    cells = cells.assign(stratum=cells.ct)
    n_per = {ct: 60 for ct in CTS}
elif design.startswith("traj_"):
    lin = design[5:]
    pt = pd.read_csv(SLING, sep="\t", index_col=0)[lin]
    key = cells.barcode.str.replace("_Lane", "_lane").str.replace("_paired", "-1")
    cells = cells.assign(pt=pt.reindex(key).to_numpy())
    cells = cells[cells.ct.isin(["iPSC"] + LINEAGE[lin]) & cells.pt.notna()].copy()
    ip = cells.ct == "iPSC"
    cells.loc[ip, "stratum"] = pd.qcut(cells.pt[ip], 2, labels=False).astype(int)
    cells.loc[~ip, "stratum"] = 2 + pd.qcut(cells.pt[~ip], 8, labels=False).astype(int)
    cells["stratum"] = [f"bin{int(s):02d}" for s in cells.stratum]
    n_per = {s: 40 for s in sorted(cells.stratum.unique())}
else:
    raise SystemExit(f"unknown run {run}")
strata = list(n_per)
say(run, "cells", len(cells), "strata", len(strata), "k", k_cells, "B", B, "unpaired", unpaired)
say(cells.groupby("stratum").size().reindex(strata).to_string())

# ---- pseudocell membership (global cell index = position in `cells`)
cells = cells.reset_index(drop=True)
stratum_arr = cells.stratum.to_numpy()
pools = {s: np.where(stratum_arr == s)[0] for s in strata}
if design.startswith("nn_"):
    from scipy.spatial import cKDTree
    tree = cKDTree(Z)


def draw(g):
    """One pseudocell set: lists of ATAC / RNA cell indices and stratum labels (RNA = ATAC cells unless unpaired)."""
    ra, rr, lb = [], [], []
    if design.startswith("nn_"):
        _, nn = tree.query(Z[g.choice(len(Z), n_per[ct], replace=False)], k=3 * k_cells)
        for row in nn:
            ra.append(row[:k_cells])
            rr.append(g.choice(row[k_cells:], k_cells, replace=False) if unpaired else row[:k_cells])
            lb.append(ct)
        return ra, rr, np.array(lb)
    for s in strata:
        kk = k_cells if design != "all" or s != "Vascular_endothelial" else 25
        for _ in range(n_per[s]):
            x = g.choice(pools[s], kk, replace=False)
            ra.append(x)
            rr.append(g.choice(pools[s], kk, replace=False) if unpaired else x)
            lb.append(s)
    return ra, rr, np.array(lb)


S = sp.csr_matrix((np.ones(len(cells), np.float32), (pd.Categorical(cells.stratum, strata).codes, np.arange(len(cells)))),
                  shape=(len(strata), len(cells)))
nstr = np.asarray(S.sum(1)).ravel()

# pass 1: detection per stratum -> feature filter
det_a = np.zeros((len(strata), len(fatac)), np.float32)
det_r = np.zeros((len(strata), len(tf_idx)), np.float32)
tot_a = np.zeros(len(fatac))
tot_r = np.zeros(len(tf_idx))
for c_, g in cells.groupby("ct", sort=False):
    gi, rr_ = g.index.to_numpy(), g.row.to_numpy()
    R = sp.load_npz(f"{SC}/{c_}/RNA.npz").tocsc()[:, tf_idx].tocsr()[rr_]
    tot_r += np.asarray(R.sum(0)).ravel()
    R.data[:] = 1
    det_r += (S[:, gi] @ R).toarray()
    del R
    A = sp.load_npz(f"{SC}/{c_}/ATAC.npz").tocsr()[rr_]
    tot_a += np.asarray(A.sum(0)).ravel()
    A.data[:] = 1
    det_a += (S[:, gi] @ A).toarray()
    del A
    say("pass 1", c_, len(gi), f"{time.time() - t0:.0f}s")
det_a /= nstr[:, None]
det_r /= nstr[:, None]
keep_a = (det_a.max(0) > 0.01) & (tot_a > 5)
keep_r = (det_r.max(0) > 0.01) & (tot_r > 5)
peaks = fatac.id.to_numpy()[keep_a]
tfs = tf_names[keep_r]
ka, kr = np.where(keep_a)[0], tf_idx[keep_r]

# pass 2: filtered single-cell matrices in memory (rows = `cells` order), to rebuild pseudocells in every replicate
parts_a, parts_r, order = [], [], []
for c_, g in cells.groupby("ct", sort=False):
    rr_ = g.row.to_numpy()
    parts_r.append(sp.load_npz(f"{SC}/{c_}/RNA.npz").tocsc()[:, kr].tocsr()[rr_])
    parts_a.append(sp.load_npz(f"{SC}/{c_}/ATAC.npz").tocsc()[:, ka].tocsr()[rr_])
    order.append(g.index.to_numpy())
    say("pass 2", c_, f"{time.time() - t0:.0f}s")
inv = np.argsort(np.concatenate(order))
Acells = sp.vstack(parts_a).tocsr()[inv]
Rcells = sp.vstack(parts_r).tocsr()[inv]
del parts_a, parts_r
say("cells x peaks nnz", Acells.nnz, "| TFs", len(tfs), "| peaks", len(peaks), f"{time.time() - t0:.0f}s")


def build(b):
    """Pseudocells for replicate b (b = -1: reference set used for the stratum summaries), TMM + log-CPM."""
    ra, rr, lb = draw(np.random.default_rng(SEED * 1000 + b + 1))
    pa = np.vstack([np.asarray(Acells[x].sum(0)).ravel() for x in ra])
    pr = np.vstack([np.asarray(Rcells[x].sum(0)).ravel() for x in rr])
    return tmm_log_cpm(pr.astype(np.float64)), tmm_log_cpm(pa.astype(np.float64)), lb


X, Y, lab = build(-1)
npc = len(lab)
say("pseudocells", npc, "TFs", len(tfs), "peaks", len(peaks), f"{time.time() - t0:.0f}s")

# per-stratum summaries (reference pseudocell set)
codes = pd.Categorical(lab, strata).codes
mean_tf = np.vstack([X[codes == i].mean(0) for i in range(len(strata))]).astype(np.float32)
mean_pk = np.vstack([Y[codes == i].mean(0) for i in range(len(strata))]).astype(np.float32)
extra = {}
if design.startswith("traj_"):
    extra["pt_median"] = cells.groupby("stratum").pt.median().reindex(strata).to_numpy()
    comp = cells.groupby(["stratum", "ct"]).size().unstack(fill_value=0).reindex(strata)
    extra["comp"], extra["comp_ct"] = comp.to_numpy(), comp.columns.to_numpy().astype(str)
np.savez_compressed(f"{od}/strata.npz", strata=np.array(strata), n_cells=nstr, tf=tfs, peak=peaks, mean_tf=mean_tf,
                    mean_peak=mean_pk, det_peak=det_a[:, keep_a].astype(np.float16), **extra)
np.savez_compressed(f"{od}/pseudocells_tf.npz", X=X.astype(np.float32), stratum=lab, tf=tfs)

# ---- SPIRE: stability selection over independently rebuilt pseudocell sets (resample = "cells")
res = bootstrap_omp(X, Y, L=L, eps=0.0, n_bootstrap=B, tau=TAU, resample="cells",
                    pseudocell_factory=lambda b: build(b)[:2], backend="numpy", seed=SEED, verbose=True)
say("bootstrap done", f"{time.time() - t0:.0f}s")
e = res.edges(tfs, peaks)
e.to_csv(f"{od}/edges.tsv.gz", sep="\t", index=False, float_format="%.4g")
cnt, frac, _ = regulatory_scope(res.consensus)
sc = pd.DataFrame({"tf": tfs, "n_peaks": cnt, "fraction": frac, "mean_expr": X.mean(0), "max_stratum_expr": mean_tf.max(0)})
sc.sort_values("n_peaks", ascending=False).to_csv(f"{od}/scope.tsv", sep="\t", index=False, float_format="%.5g")
F = sp.csr_matrix(np.where(res.frequency >= 0.1, res.frequency, 0).astype(np.float32))
sp.save_npz(f"{od}/freq.npz", F)
ever = int((res.frequency > 0).sum())
say("edges selected >= once", ever, "| consensus", len(e), f"({len(e) / max(ever, 1):.4f})",
    "| peaks with >= 1 TF", e.peak.nunique(), f"of {len(peaks)}", "| TFs with >= 1 target", e.tf.nunique(),
    "| frac positive", round((e.coef > 0).mean(), 3), "| mean TFs per peak", round(len(e) / len(peaks), 3))
say("top scope:", ", ".join(sc.sort_values("n_peaks", ascending=False).tf.head(25)))
say("done", f"{time.time() - t0:.0f}s")
