"""
QTL x GWAS colocalization (coloc.abf, vectorized) and SMR-type test, one cell type x modality x chromosome,
for every harmonized GWAS trait (06b_harmonize_gwas.py).

coloc.abf: identical to coloc::coloc.abf (approx.bf.estimates + combine.abf, priors p1=p2=1e-4, p12=1e-5)
  on the SNPs shared by the QTL feature's cis window and the GWAS.
  QTL: z from the nominal t-test p-value (z = qnorm(p/2), sign of slope) and varbeta = (slope/z)^2, so the
       small-sample (n=17, t-distributed) evidence is not overstated; sd.prior 0.15 (sdY = 1, invnorm phenotype).
  GWAS: beta, se; sd.prior 0.2 (case-control) or 0.15 * sdY (quantitative, sdY estimated in 06b).
  A feature x trait is reported when the GWAS has p < 1e-6 among the shared SNPs and >= 50 shared SNPs.
SMR (Zhu et al. 2016, no HEIDI): at each feature's lead QTL variant (permutation output),
  T = z_gwas^2 z_qtl^2 / (z_gwas^2 + z_qtl^2) ~ chi2_1, b_xy = b_gwas / b_qtl (effect per SD of feature).

Outputs in HDC_Filtered/GWAS/coloc/<MOD>/<CT>/chr<N>.*:
  coloc.parquet   feature x trait coloc results (PP.H0-H4, best causal SNP and its SNP.PP.H4, min p's)
  lead.parquet    per feature lead variant + GWAS z for all traits (wide; for enrichment / SMR)
Usage: python 07_coloc_smr.py <CellType> <RNA|ATAC> <chr>
"""
import os
import sys
import time

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import norm

BASE = os.environ["HDC_DIR"] + "/"
TAG = "noBootstrap_TMM_invnorm"
P1, P2, P12 = 1e-4, 1e-4, 1e-5
GWAS_Z2_MIN = 23.928  # chi2_1 for p = 1e-6
MIN_SNPS = 50

celltype, modality, chrom = sys.argv[1], sys.argv[2], int(sys.argv[3])
qdir = f"{BASE}/{'eQTL' if modality == 'RNA' else 'caQTL'}_Output_optCov/{celltype}/"
out_dir = f"{BASE}/GWAS/coloc/{modality}/{celltype}/"
os.makedirs(out_dir, exist_ok=True)
t0 = time.time()


def labf(z, V, W):
    r = W / (W + V)
    return 0.5 * (np.log1p(-r) + r * z * z)


def z_from_p(p, slope):
    z = norm.isf(np.clip(p, 1e-300, 1) / 2) * np.sign(slope)
    return z


# --- variant universe for this chromosome (vidx contiguous per chr)
u = pd.read_parquet(f"{BASE}/GWAS/variant_universe.parquet", columns=["vidx", "chr", "pos38", "id"],
                    filters=[("chr", "==", chrom)])
off, nvar = int(u.vidx.iloc[0]), len(u)
assert (u.vidx.to_numpy() == off + np.arange(nvar)).all()
pos = u.pos38.to_numpy(np.int64)
id_index = pd.Index(u.id)

# --- QTL nominal statistics
tab = pq.read_table(f"{qdir}/chr{chrom}_cis_nominal_{TAG}.cis_qtl_pairs.{chrom}.parquet",
                    columns=["phenotype_id", "variant_id", "pval_nominal", "slope"],
                    read_dictionary=["phenotype_id", "variant_id"])
ph = tab.column("phenotype_id").combine_chunks()
va = tab.column("variant_id").combine_chunks()
ph_codes, ph_names = ph.indices.to_numpy(), np.array(ph.dictionary.to_pylist())
v_map = id_index.get_indexer(va.dictionary.to_pylist())
v = v_map[va.indices.to_numpy()]
pv = tab.column("pval_nominal").to_numpy()
slope = tab.column("slope").to_numpy().astype(np.float64)
del tab, ph, va
keep = (v >= 0) & np.isfinite(pv) & np.isfinite(slope) & (slope != 0)
ph_codes, v, pv, slope = ph_codes[keep], v[keep], pv[keep], slope[keep]
# rows are written phenotype by phenotype; make segments explicit
brk = np.flatnonzero(np.diff(ph_codes)) + 1
seg_start = np.r_[0, brk]
seg_ph = ph_codes[seg_start]
assert len(np.unique(seg_ph)) == len(seg_ph), "nominal rows not grouped by phenotype"
seg_len = np.diff(np.r_[seg_start, len(v)])
zq = z_from_p(pv, slope)
l1 = labf(zq, (slope / zq) ** 2, 0.15 ** 2)
zq2 = zq * zq
vpos = pos[v]
seg_min = np.minimum.reduceat(vpos, seg_start)
seg_max = np.maximum.reduceat(vpos, seg_start)
print(f"{celltype} {modality} chr{chrom}: {len(v)} pairs, {len(seg_start)} features ({time.time() - t0:.0f}s)", flush=True)

# --- lead variants (permutation)
leads = pd.read_csv(f"{qdir}/chr{chrom}_cis_{TAG}.cis_qtl.txt.gz", sep="\t",
                    usecols=["phenotype_id", "variant_id", "start_distance", "af", "pval_nominal", "slope", "pval_beta"])
q = pd.read_csv(f"{qdir}/allchr_cis_{TAG}.cis_qtl.qval.txt.gz", sep="\t", usecols=["phenotype_id", "qval"])
leads = leads.merge(q, on="phenotype_id", how="left")
leads["lvidx"] = id_index.get_indexer(leads.variant_id)
leads = leads[leads.lvidx >= 0].reset_index(drop=True)
leads["z_qtl"] = z_from_p(leads.pval_nominal.to_numpy(), leads.slope.to_numpy())
lead_out = leads[["phenotype_id", "variant_id", "start_distance", "af", "pval_beta", "qval", "slope",
                  "z_qtl"]].copy()
lead_out["pos38"] = pos[leads.lvidx]

manifest = pd.read_csv(f"{BASE}/GWAS/gwas_manifest.tsv", sep="\t")
coloc_rows, smr_cols = [], {}
for _, m in manifest.iterrows():
    fn = f"{BASE}/GWAS/harmonized/{m.trait_id}.parquet"
    if not os.path.exists(fn):  # population not available for this trait (see harmonized/qc)
        continue
    g = pd.read_parquet(fn, filters=[("chr", "==", chrom)], columns=["vidx", "beta", "se"])
    if len(g) == 0:
        continue
    qc = pd.read_csv(f"{BASE}/GWAS/harmonized/qc/{m.trait_id}.tsv", sep="\t").iloc[0]
    sdy = qc.sdY if (m.type == "quant" and np.isfinite(qc.sdY)) else 1.0
    W2 = 0.2 ** 2 if m.type == "cc" else (0.15 * sdy) ** 2
    li = g.vidx.to_numpy() - off
    gb, gse = np.full(nvar, np.nan), np.full(nvar, np.nan)
    gb[li], gse[li] = g.beta.to_numpy(np.float64), g.se.to_numpy(np.float64)
    gz = gb / gse
    gl = labf(gz, gse * gse, W2)

    # lead-variant GWAS z (for enrichment and SMR)
    smr_cols[m.trait_id] = gz[leads.lvidx].astype(np.float32)

    # coloc: features whose cis window contains a GWAS p < 1e-6 variant
    hit = np.sort(pos[np.isfinite(gz) & (gz * gz > GWAS_Z2_MIN)])
    if len(hit) == 0:
        continue
    j = np.searchsorted(hit, seg_min)
    fsel = np.flatnonzero((j < len(hit)) & (hit[np.minimum(j, len(hit) - 1)] <= seg_max))
    if len(fsel) == 0:
        continue
    starts, lens = seg_start[fsel], seg_len[fsel]
    rows = np.concatenate([np.arange(s, s + n) for s, n in zip(starts, lens)])
    new_start = np.r_[0, np.cumsum(lens)[:-1]]
    a1 = l1[rows]
    a2 = gl[v[rows]]
    ok = np.isfinite(a2)
    a1 = np.where(ok, a1, -np.inf)
    a2 = np.where(ok, a2, -np.inf)
    a12 = a1 + a2
    nsnp = np.add.reduceat(ok.astype(np.int64), new_start)

    def lse(x):
        mx = np.maximum.reduceat(x, new_start)
        mxs = np.where(np.isfinite(mx), mx, 0.0)
        return mxs + np.log(np.add.reduceat(np.exp(x - np.repeat(mxs, lens)), new_start))

    s1, s2, s12 = lse(a1), lse(a2), lse(a12)
    gz2 = np.where(ok, gz[v[rows]] ** 2, 0.0)
    gmax = np.maximum.reduceat(gz2, new_start)
    qmax = np.maximum.reduceat(np.where(ok, zq2[rows], 0.0), new_start)
    good = (nsnp >= MIN_SNPS) & (gmax > GWAS_Z2_MIN)
    if not good.any():
        continue
    # combine.abf
    with np.errstate(divide="ignore", invalid="ignore"):
        lsum = s1 + s2
        h3 = np.log(P1) + np.log(P2) + lsum + np.log1p(-np.exp(s12 - lsum))
    h = np.vstack([np.zeros_like(s1), np.log(P1) + s1, np.log(P2) + s2, np.where(np.isfinite(h3), h3, -np.inf),
                   np.log(P12) + s12])
    hm = h.max(axis=0)
    pp = np.exp(h - hm)
    pp /= pp.sum(axis=0)
    # best causal SNP for H4
    seg_id = np.repeat(np.arange(len(fsel)), lens)
    mx12 = np.repeat(np.maximum.reduceat(a12, new_start), lens)
    first = pd.Series(np.flatnonzero(a12 == mx12)).groupby(seg_id[a12 == mx12]).first()
    best_row = np.full(len(fsel), -1)
    best_row[first.index.to_numpy()] = first.to_numpy()
    br = rows[np.maximum(best_row, 0)]
    gi = good & (best_row >= 0)
    coloc_rows.append(pd.DataFrame({
        "trait_id": m.trait_id, "phenotype_id": ph_names[seg_ph[fsel[gi]]], "nsnps": nsnp[gi],
        "PP.H0": pp[0, gi], "PP.H1": pp[1, gi], "PP.H2": pp[2, gi], "PP.H3": pp[3, gi], "PP.H4": pp[4, gi],
        "best_snp": u.id.to_numpy()[v[br[gi]]], "best_pos38": pos[v[br[gi]]],
        "best_snp_ppH4": np.exp(a12[best_row[gi]] - s12[gi]),
        "z_qtl_best": zq[br[gi]], "z_gwas_best": gz[v[br[gi]]],
        "qtl_minp": 2 * norm.sf(np.sqrt(qmax[gi])), "gwas_minp": 2 * norm.sf(np.sqrt(gmax[gi]))}))

if coloc_rows:
    res = pd.concat(coloc_rows, ignore_index=True)
    res.insert(0, "chr", chrom)
    for c in res.select_dtypes("float64").columns:
        if c not in ("qtl_minp", "gwas_minp"):
            res[c] = res[c].astype(np.float32)
    res.to_parquet(f"{out_dir}/chr{chrom}.coloc.parquet.tmp", index=False)
    os.replace(f"{out_dir}/chr{chrom}.coloc.parquet.tmp", f"{out_dir}/chr{chrom}.coloc.parquet")
    print(f"coloc rows {len(res)}; PP.H4>0.5: {(res['PP.H4'] > 0.5).sum()}; >0.8: {(res['PP.H4'] > 0.8).sum()}")
lead_out = pd.concat([lead_out, pd.DataFrame(smr_cols)], axis=1)
lead_out.to_parquet(f"{out_dir}/chr{chrom}.lead.parquet.tmp", index=False)
os.replace(f"{out_dir}/chr{chrom}.lead.parquet.tmp", f"{out_dir}/chr{chrom}.lead.parquet")
print(f"done ({time.time() - t0:.0f}s)", flush=True)
