"""
Cell-type-specific peak-gene links with scMultiMap (NumPy port in scmultimap_py.py) on single-cell multiome counts
(PeakAnnot/sc/<CT>/ from 22_extract_ct_matrices.py).

Per cell type:
  genes   protein_coding / lncRNA (GENCODE v32 = cellranger-arc 2020-A annotation, matched by Ensembl id),
          chr1-22, X, detected (count > 0) in >= MINFRAC of the cell type's cells
  peaks   chr1-22, X, detected in >= MINFRAC of the cell type's cells
  pairs   peak center within +-500 kb of the gene TSS (scMultiMap get_top_peak_gene_pairs distance)
  bsample donor (17 YRI lines; per-donor means, so links are within-donor cell-to-cell covariation)
  seq depth = total counts over all features of the modality
FDR = BH over all pairs of the cell type (scMultiMap default).
--perm: RNA cells permuted within donor (null calibration).

Output: PeakAnnot/scMultiMap/<CT>{_perm}.pairs.parquet (all tested pairs), <CT>{_perm}.irls_{RNA,ATAC}.tsv.gz,
        <CT>{_perm}.log via SLURM
Usage: python 23_scmultimap.py <CT or AllCells> [--perm]
"""
import os
import sys
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scmultimap_py import irls, wls_block  # noqa: E402

BASE = os.environ["HDC_DIR"] + "/"
SC = f"{BASE}/PeakAnnot/sc/"
OUT = f"{BASE}/PeakAnnot/scMultiMap/"
GTF = f"{BASE}/PeakAnnot/ref/gencode.v32.basic.annotation.gtf.gz"
GENES = f"{BASE}/PeakAnnot/ref/genes_gencode32.tsv.gz"
MINFRAC = float(os.environ.get("MINFRAC", 0.01))
# optional stricter pre-filter (used for AllCells to keep IRLS tractable; equals the post-hoc filter of 26_link_analysis.py)
MINFRAC_GENE = float(os.environ.get("MINFRAC_GENE", MINFRAC))
MINFRAC_PEAK = float(os.environ.get("MINFRAC_PEAK", MINFRAC))
MIN_CELLS = int(os.environ.get("MIN_CELLS", 0))
DIST = 500_000
GBLOCK = 64
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]


def bh(p):
    """BH adjusted p over non-missing p (NaN stays NaN)"""
    q = np.full(len(p), np.nan)
    ok = np.where(~np.isnan(p))[0]
    o = ok[np.argsort(p[ok])]
    m = len(ok)
    q[o] = np.minimum(np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1], 1)
    return q


def gene_table():
    if os.path.exists(GENES):
        return pd.read_csv(GENES, sep="\t")
    g = pd.read_csv(GTF, sep="\t", comment="#", header=None, usecols=[0, 2, 3, 4, 6, 8],
                    names=["chr", "type", "start", "end", "strand", "attr"])
    g = g[g.type == "gene"]
    a = g.attr.str.extract(r'gene_id "([^".]+)[^"]*";.*gene_type "([^"]+)";.*gene_name "([^"]+)"')
    g = pd.DataFrame({"gene_id": a[0], "gene_type": a[1], "gene_name": a[2], "chr": g.chr, "start": g.start,
                      "end": g.end, "strand": g.strand})
    g["tss"] = np.where(g.strand == "+", g.start, g.end)
    g.to_csv(GENES, sep="\t", index=False)
    return g


def main():
    ct = sys.argv[1]
    perm = "--perm" in sys.argv
    tag = ct + ("_perm" if perm else "")
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    # AllCells: all cell types pooled (bulk-like, cell-type composition included), for comparison
    dirs = sorted(d for d in os.listdir(SC) if os.path.isdir(f"{SC}/{d}")) if ct == "AllCells" else [ct]
    cells = pd.concat([pd.read_csv(f"{SC}/{d}/cells.tsv.gz", sep="\t") for d in dirs], ignore_index=True)
    R = sp.vstack([sp.load_npz(f"{SC}/{d}/RNA.npz") for d in dirs]).tocsc()
    A = sp.vstack([sp.load_npz(f"{SC}/{d}/ATAC.npz") for d in dirs]).tocsc()
    fr = pd.read_csv(f"{SC}/features_RNA.tsv.gz", sep="\t")
    fa = pd.read_csv(f"{SC}/features_ATAC.tsv.gz", sep="\t")
    n = len(cells)
    donors = pd.Index(sorted(cells.donor.unique()))
    bcode = donors.get_indexer(cells.donor)
    sd_g = cells.depth_RNA.to_numpy(float)
    sd_p = cells.depth_ATAC.to_numpy(float)
    if perm:
        rng = np.random.default_rng(2025)
        idx = np.arange(n)
        for b in range(len(donors)):
            k = np.where(bcode == b)[0]
            idx[k] = rng.permutation(k)
        R = R[idx].tocsc()
        sd_g = sd_g[idx]

    # feature filters
    gt = gene_table()
    gt = gt[gt.gene_type.isin(["protein_coding", "lncRNA"]) & gt.chr.isin(CHROMS)].drop_duplicates("gene_id")
    fr["col"] = np.arange(len(fr))
    fr = fr.merge(gt[["gene_id", "gene_type", "chr", "tss", "strand"]].rename(columns={"chr": "gchr"}), on="gene_id")
    fr = fr[fr.gchr == fr.chr]
    det_g = np.diff(R.indptr) / n
    fr = fr[(det_g[fr.col] >= MINFRAC_GENE) & (det_g[fr.col] * n >= MIN_CELLS)].reset_index(drop=True)
    fa["col"] = np.arange(len(fa))
    fa = fa[fa.chr.isin(CHROMS)]
    det_p = np.diff(A.indptr) / n
    fa = fa[(det_p[fa.col] >= MINFRAC_PEAK) & (det_p[fa.col] * n >= MIN_CELLS)].reset_index(drop=True)
    fa["center"] = (fa.start + fa.end) // 2
    print(f"{tag}: {n} cells, {len(donors)} donors; {len(fr)} genes (>= {MINFRAC_GENE}), {len(fa)} peaks (>= {MINFRAC_PEAK} of cells)", flush=True)

    Rg = R[:, fr.col.to_numpy()].tocsc()
    Ap = A[:, fa.col.to_numpy()].tocsc()
    del R, A
    pg = irls(Rg, sd_g, bcode, len(donors), verbose=True)
    print(f"IRLS RNA done ({pg['iters']} iters), {time.time() - t0:.0f}s", flush=True)
    pp = irls(Ap, sd_p, bcode, len(donors), verbose=True)
    print(f"IRLS ATAC done ({pp['iters']} iters), {time.time() - t0:.0f}s", flush=True)
    for mod, p, f in (("RNA", pg, fr.id), ("ATAC", pp, fa.id)):
        d = pd.DataFrame(p["mu"], index=f, columns=[f"mu_{x}" for x in donors])
        d.insert(0, "sigma_sq", p["sigma_sq"])
        d.to_csv(f"{OUT}/{tag}.irls_{mod}.tsv.gz", sep="\t", float_format="%.5g")
        pd.Series(p["theta"], index=donors).to_csv(f"{OUT}/{tag}.theta_{mod}.tsv", sep="\t", header=["theta"])

    res = []
    for chrom in CHROMS:
        gi = np.where(fr.chr == chrom)[0]
        gi = gi[np.argsort(fr.tss.to_numpy()[gi])]
        pi = np.where(fa.chr == chrom)[0]
        pc = fa.center.to_numpy()[pi]
        for b0 in range(0, len(gi), GBLOCK):
            g = gi[b0:b0 + GBLOCK]
            tss = fr.tss.to_numpy()[g]
            k = pi[(pc >= tss.min() - DIST) & (pc <= tss.max() + DIST)]
            if len(k) == 0:
                continue
            xg = Rg[:, g].toarray().astype(np.float64)
            xp = Ap[:, k].toarray().astype(np.float64)
            est, ts, pv = wls_block(xg, sd_g, {"mu": pg["mu"][g], "sigma_sq": pg["sigma_sq"][g], "theta": pg["theta"]},
                                    xp, sd_p, {"mu": pp["mu"][k], "sigma_sq": pp["sigma_sq"][k], "theta": pp["theta"]}, bcode)
            dist = fa.center.to_numpy()[k][None, :] - tss[:, None]
            ii, jj = np.where(np.abs(dist) <= DIST)
            if len(ii) == 0:
                continue
            sgn = np.where(fr.strand.to_numpy()[g][ii] == "+", 1, -1)
            vg, vp = pg["sigma_sq"][g][ii], pp["sigma_sq"][k][jj]
            with np.errstate(invalid="ignore", divide="ignore"):
                cor = est[ii, jj] / np.sqrt(np.where(vg > 0, vg, np.nan) * np.where(vp > 0, vp, np.nan))
            res.append(pd.DataFrame({"gene": fr.id.to_numpy()[g][ii], "gene_id": fr.gene_id.to_numpy()[g][ii],
                                     "peak": fa.id.to_numpy()[k][jj], "chr": chrom,
                                     "dist_tss": (dist[ii, jj] * sgn).astype(np.int32),
                                     "pval": pv[ii, jj], "test_stat": ts[ii, jj].astype(np.float32),
                                     "covar": est[ii, jj].astype(np.float32),
                                     "cor": np.clip(np.nan_to_num(cor, nan=0.0), -1, 1).astype(np.float32)}))
        print(f"{chrom}: {sum(len(r) for r in res)} pairs so far, {time.time() - t0:.0f}s", flush=True)
    res = pd.concat(res, ignore_index=True)
    res["padj"] = bh(res.pval.to_numpy())
    res.to_parquet(f"{OUT}/{tag}.pairs.parquet", index=False)
    print(f"{tag}: {len(res)} pairs; FDR<0.1 {int((res.padj < 0.1).sum())}, FDR<0.05 {int((res.padj < 0.05).sum())} "
          f"(positive {int(((res.padj < 0.05) & (res.covar > 0)).sum())}); {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
