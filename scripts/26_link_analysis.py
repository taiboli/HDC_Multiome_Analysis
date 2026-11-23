"""
Characterize scMultiMap peak-gene links (23_scmultimap.py) and test their biological support.

Feature filter (post hoc, from the within-donor permutation nulls; 26c_scmultimap_calibration.py ->
PeakAnnot/results/scmultimap_calibration.tsv): gene detected in >= GENE_DET and peak in >= PEAK_DET of the cell type's
cells, both in >= MIN_CELLS cells; BH re-computed on the retained pairs. (Sparse features give a positive heavy-tailed
null; with these thresholds permuted / real FDR-10% links = 1.5% (Neuronal), 4.8% (iPSC).) Link = padj < FDR (both signs kept; positive links are the regulatory candidates).

  summary        tested / linked pairs, peaks, genes per cell type; distance; nearest-gene agreement
  distance       link rate by TSS distance bin
  peak_enrich    linked vs unlinked tested peaks: logistic regression of link status on each annotation (ENCODE cCRE,
                 enhancer state in the cell type's best-matched Roadmap epigenome and in H1, phastCons, CT-specific
                 peak, not-in-registry), adjusting for log CPM in the cell type, GC, log width, log10 TSS distance
                 and log # candidate genes
  nearest_gene   for distal linked peaks whose link skips the nearest gene: PCHi-C / GTEx support of linked vs nearest gene
  pair_valid     GTEx v8 (eVariant of the gene inside the peak, any tissue) and PCHi-C (Jung 2019; peak in a fragment
                 interacting with the gene promoter in H1 / H1-derived ME, MSC, NPC, TB, or any cell): linked vs
                 unlinked pairs, Mantel-Haenszel OR over TSS-distance bins
  qtl_valid      HDC caQTL lead SNP (pval_beta < 0.05) of the linked peak: eQTL of the linked gene at that SNP in the
                 same cell type (fraction p < 0.05) and sign(caQTL slope x eQTL slope) = sign(link)
  sharing        links of cell type A replicated (p < 0.05, same sign; and pi1) in cell type B; AllCells pooled recovery
  gene_spec      linked genes / peaks vs cell-type-specific genes / peaks
  gwas_targets   scMultiMap target genes of trait-colocalized caQTL peaks (GWAS/results/coloc_H4gt0.5.tsv, PP.H4 > 0.5)

Output: PeakAnnot/results/links_fdr<FDR>.tsv.gz, link_*.tsv; figures in 27_link_figures.py
"""
import glob
import os

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import scipy.sparse as sp
from scipy.stats import binomtest, chi2, norm

BASE = os.environ["HDC_DIR"] + "/"
PA = f"{BASE}/PeakAnnot/"
SM, RES = f"{PA}/scMultiMap/", f"{PA}/results/"
TAG = "noBootstrap_TMM_invnorm"
GENE_DET = float(os.environ.get("GENE_DET", 0.10))
PEAK_DET = float(os.environ.get("PEAK_DET", 0.02))
MIN_CELLS = int(os.environ.get("MIN_CELLS", 100))
FDR = float(os.environ.get("FDR", 0.1))
CTS = ["iPSC", "Ectoderm_1", "Ectoderm_2", "Ectoderm_3", "Neural_crest", "Neuronal", "Retinal", "Mesoderm",
       "Vascular_endothelial", "Endoderm", "Unknown"]
ENH = [7, 8, 9, 10, 11]
ES_PCHIC = ["H1", "ME", "MSC", "NPC", "TB"]
DBINS = [0, 2_000, 10_000, 50_000, 100_000, 250_000, 500_001]
DLAB = ["0-2k", "2-10k", "10-50k", "50-100k", "100-250k", "250-500k"]
os.makedirs(RES, exist_ok=True)


def bh(p):
    q = np.full(len(p), np.nan)
    ok = np.where(~np.isnan(p))[0]
    o = ok[np.argsort(p[ok])]
    m = len(ok)
    q[o] = np.minimum(np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1], 1)
    return q


def logit(y, X, it=50):
    b = np.zeros(X.shape[1])
    for _ in range(it):
        mu = 1 / (1 + np.exp(-np.clip(X @ b, -30, 30)))
        W = mu * (1 - mu)
        H = X.T @ (X * W[:, None])
        step = np.linalg.solve(H, X.T @ (y - mu))
        b += step
        if np.max(np.abs(step)) < 1e-8:
            break
    return b, np.sqrt(np.diag(np.linalg.inv(H)))


def mh(y, x, strata):
    """Mantel-Haenszel OR of binary y on binary x over strata, with CMH chi2 p"""
    num = den = 0.0
    a_sum = e_sum = v_sum = 0.0
    for s in np.unique(strata):
        k = strata == s
        a = np.sum(y[k] & x[k]); b = np.sum(~y[k] & x[k]); c = np.sum(y[k] & ~x[k]); d = np.sum(~y[k] & ~x[k])
        n = a + b + c + d
        if n < 2:
            continue
        num += a * d / n
        den += b * c / n
        r1, c1 = a + b, a + c
        a_sum += a
        e_sum += r1 * c1 / n
        v_sum += r1 * (n - r1) * c1 * (n - c1) / (n * n * (n - 1))
    orr = num / den if den > 0 else np.nan
    p = chi2.sf((abs(a_sum - e_sum) - 0.5) ** 2 / v_sum, 1) if v_sum > 0 else np.nan
    return orr, p


def detection(ct):
    if ct == "AllCells":
        dirs = [d for d in CTS]
        n = sum(len(pd.read_csv(f"{PA}/sc/{d}/cells.tsv.gz", sep="\t")) for d in dirs)
        out = {}
        for mod in ("RNA", "ATAC"):
            out[mod] = sum(np.diff(sp.load_npz(f"{PA}/sc/{d}/{mod}.npz").indptr) for d in dirs)
        return n, out
    n = len(pd.read_csv(f"{PA}/sc/{ct}/cells.tsv.gz", sep="\t"))
    return n, {mod: np.diff(sp.load_npz(f"{PA}/sc/{ct}/{mod}.npz").indptr) for mod in ("RNA", "ATAC")}


fr = pd.read_csv(f"{PA}/sc/features_RNA.tsv.gz", sep="\t")
fa = pd.read_csv(f"{PA}/sc/features_ATAC.tsv.gz", sep="\t")


def load_pairs(ct):
    d = pd.read_parquet(f"{SM}/{ct}.pairs.parquet")
    n, det = detection(ct.replace("_perm", ""))
    gdet = pd.Series(det["RNA"], index=fr.id)
    pdet = pd.Series(det["ATAC"], index=fa.id)
    d["gene_ncell"] = gdet.reindex(d.gene).to_numpy()
    d["peak_ncell"] = pdet.reindex(d.peak).to_numpy()
    d = d[(d.gene_ncell >= max(GENE_DET * n, MIN_CELLS)) & (d.peak_ncell >= max(PEAK_DET * n, MIN_CELLS))].copy()
    d["padj"] = bh(d.pval.to_numpy())
    d["celltype"] = ct
    return d, n


# ---------------- load ----------------
pk = pd.read_parquet(f"{PA}/peak_annotation.parquet")
st = pd.read_parquet(f"{PA}/peak_chromhmm18.parquet").set_index("peak")
best_e = pd.read_csv(f"{RES}/ctspec_peaks_top_epigenomes.tsv", sep="\t").groupby("celltype").eid.first()
pki = pk.set_index("peak")
cts = [c for c in CTS if os.path.exists(f"{SM}/{c}.pairs.parquet")]
P, ncell = {}, {}
for ct in cts + ["AllCells"]:
    if os.path.exists(f"{SM}/{ct}.pairs.parquet"):
        P[ct], ncell[ct] = load_pairs(ct)
        print(ct, len(P[ct]), "pairs after filter", int((P[ct].padj < FDR).sum()), "links", flush=True)

# ---------------- summary ----------------
rows, links = [], []
for ct in P:
    d = P[ct]
    L = d[d.padj < FDR]
    links.append(L)
    lp = L[L.covar > 0]
    near = pki.nearest_gene.reindex(lp.peak).to_numpy() == lp.gene.to_numpy()
    distal = (pki.category.reindex(lp.peak) == "distal").to_numpy()
    rows.append({"celltype": ct, "n_cells": ncell[ct], "tested_pairs": len(d), "tested_genes": d.gene.nunique(),
                 "tested_peaks": d.peak.nunique(), "links": len(L), "links_pos": len(lp), "frac_pos": len(lp) / max(len(L), 1),
                 "linked_genes": lp.gene.nunique(), "linked_peaks": lp.peak.nunique(),
                 "median_abs_dist": lp.dist_tss.abs().median(), "frac_link_nearest_gene": near.mean() if len(lp) else np.nan,
                 "frac_linked_peaks_distal": distal.mean() if len(lp) else np.nan,
                 "median_peaks_per_gene": lp.groupby("gene").size().median() if len(lp) else np.nan,
                 "median_genes_per_peak": lp.groupby("peak").size().median() if len(lp) else np.nan})
summ = pd.DataFrame(rows)
summ.to_csv(f"{RES}/link_summary.tsv", sep="\t", index=False, float_format="%.4g")
print(summ.to_string(), flush=True)
links = pd.concat(links, ignore_index=True)
links.to_csv(f"{RES}/links_fdr{FDR}.tsv.gz", sep="\t", index=False, float_format="%.4g")

# ---------------- calibration (permuted nulls) ----------------
cal = []
for ct in ["Neuronal_perm", "iPSC_perm"]:
    if not os.path.exists(f"{SM}/{ct}.pairs.parquet"):
        continue
    d, _ = load_pairs(ct)
    p = d.pval.dropna().to_numpy()
    cal.append({"set": ct, "pairs": len(p), "lambda": np.median(chi2.isf(p, 1)) / chi2.ppf(0.5, 1),
                "frac_p_lt_1e-3": (p < 1e-3).mean(), "n_padj_lt_fdr": int((d.padj < FDR).sum())})
pd.DataFrame(cal).to_csv(f"{RES}/scmultimap_perm_after_filter.tsv", sep="\t", index=False, float_format="%.4g")
print(pd.DataFrame(cal).to_string(), flush=True)

# ---------------- distance decay ----------------
rows = []
for ct in cts:
    d = P[ct]
    b = pd.cut(d.dist_tss.abs(), DBINS, right=False, labels=DLAB)
    g = d.assign(bin=b, link=(d.padj < FDR) & (d.covar > 0)).groupby("bin", observed=True).agg(n=("link", "size"), rate=("link", "mean"))
    rows.append(g.reset_index().assign(celltype=ct))
pd.concat(rows).to_csv(f"{RES}/link_distance.tsv", sep="\t", index=False, float_format="%.4g")

# ---------------- peak-level enrichment ----------------
rows = []
for ct in cts:
    if ct == "Unknown":
        continue
    d = P[ct]
    tp = d.groupby("peak").agg(linked=("padj", lambda x: bool(((x < FDR)).any())), ncand=("gene", "size"))
    pos = d[(d.padj < FDR) & (d.covar > 0)].peak.unique()
    tp["linked"] = tp.index.isin(pos)
    a = pki.reindex(tp.index)
    e = best_e.get(ct, "E003")
    ann = {"cCRE_any": (a.ccre != "none"), "cCRE_dELS": a.ccre == "dELS", "cCRE_pELS": a.ccre == "pELS",
           "cCRE_PLS": a.ccre == "PLS", "not_in_registry": a.ccre == "none",
           f"enhancer_{e}_matched": st[e].reindex(tp.index).isin(ENH), "enhancer_E003_H1": st["E003"].reindex(tp.index).isin(ENH),
           "phastCons": a.phastcons_bp > 0, "CT_specific_peak": a[f"spec_{ct}"].astype(bool),
           "promoter": a.category == "promoter", "distal": a.category == "distal"}
    cov = np.column_stack([np.ones(len(tp)), np.log1p(a[f"cpm_{ct}"].to_numpy()), a.gc.to_numpy(), np.log(a.width.to_numpy()),
                           np.log10(np.abs(a.dist_tss.to_numpy()) + 1), np.log(tp.ncand.to_numpy())])
    y = tp.linked.to_numpy(float)
    for name, x in ann.items():
        x = np.asarray(x, float)
        if x.std() == 0:
            continue
        cv = cov if name not in ("promoter", "distal") else np.delete(cov, 4, axis=1)
        b, se = logit(y, np.column_stack([cv, x]))
        rows.append({"celltype": ct, "annotation": name.replace(f"_{e}_", "_"), "epigenome": e if "matched" in name else "",
                     "frac_linked": x[y == 1].mean(), "frac_unlinked": x[y == 0].mean(), "n_linked": int(y.sum()),
                     "adj_OR": np.exp(b[-1]), "lo": np.exp(b[-1] - 1.96 * se[-1]), "hi": np.exp(b[-1] + 1.96 * se[-1]),
                     "p": 2 * norm.sf(abs(b[-1] / se[-1]))})
pe = pd.DataFrame(rows)
pe.to_csv(f"{RES}/link_peak_enrichment.tsv", sep="\t", index=False, float_format="%.4g")
print(pe.to_string(), flush=True)

# ---------------- pair-level validation: GTEx, PCHi-C ----------------
gt = pd.read_parquet(f"{PA}/ref/gtex_peak_gene.parquet")
gt_key = set(gt.peak + "|" + gt.gene_id)
pc = pd.read_csv(f"{PA}/ref/pchic/pchic_po_hg38.tsv.gz", sep="\t")
pcs = []
fs = fa[fa.chr.str.match(r"^chr[0-9X]+$")]
for c, x in pc.groupby("chr"):
    p = fs[fs.chr == c].sort_values("start")
    j0 = np.searchsorted(p.end.to_numpy(), x.start.to_numpy(), side="right")
    j1 = np.searchsorted(p.start.to_numpy(), x.end.to_numpy(), side="left")
    k = j1 - j0
    rep = np.repeat(np.arange(len(x)), np.maximum(k, 0))
    idx = np.concatenate([np.arange(a, b) for a, b in zip(j0, j1) if b > a]) if (k > 0).any() else np.array([], int)
    pcs.append(pd.DataFrame({"peak": p.id.to_numpy()[idx], "gene": x.gene.to_numpy()[rep], "cell": x.cell.to_numpy()[rep]}))
pcs = pd.concat(pcs)
pc_any = set(pcs.peak + "|" + pcs.gene)
pc_es = set(pcs[pcs.cell.isin(ES_PCHIC)].peak + "|" + pcs[pcs.cell.isin(ES_PCHIC)].gene)
pcs_genes = set(pcs.gene)
rows = []
for ct in list(cts) + ["AllCells"]:
    if ct not in P or ct == "Unknown":
        continue
    d = P[ct]
    d = d[d.dist_tss.abs() >= 2000]  # promoter-proximal peaks excluded (trivially near the TSS)
    strata = pd.cut(d.dist_tss.abs(), DBINS, right=False, labels=False).to_numpy()
    link = ((d.padj < FDR) & (d.covar > 0)).to_numpy()
    null = (d.padj > 0.5).to_numpy()
    keys_g = (d.peak + "|" + d.gene_id).to_numpy()
    keys_n = (d.peak + "|" + d.gene).to_numpy()
    for name, y, sub in (("GTEx_eQTL_in_peak", pd.Series(keys_g).isin(gt_key).to_numpy(), np.ones(len(d), bool)),
                         ("PCHiC_ES_derived", pd.Series(keys_n).isin(pc_es).to_numpy(), d.gene.isin(pcs_genes).to_numpy()),
                         ("PCHiC_any", pd.Series(keys_n).isin(pc_any).to_numpy(), d.gene.isin(pcs_genes).to_numpy())):
        k = sub & (link | null)
        orr, p = mh(y[k], link[k], strata[k])
        rows.append({"celltype": ct, "support": name, "n_link": int((link & sub).sum()), "rate_link": y[link & sub].mean(),
                     "rate_null": y[null & sub].mean(), "MH_OR": orr, "p": p})
pv = pd.DataFrame(rows)
pv.to_csv(f"{RES}/link_pair_validation.tsv", sep="\t", index=False, float_format="%.4g")
print(pv.to_string(), flush=True)

# nearest gene vs scMultiMap-linked gene for the same distal peaks: GTEx support
rows = []
for ct in cts:
    if ct == "Unknown":
        continue
    d = P[ct]
    lp = d[(d.padj < FDR) & (d.covar > 0) & (d.dist_tss.abs() >= 2000)]
    lp = lp[lp.peak.map(pki.category) != "promoter"]
    ng = pki.nearest_gene.reindex(lp.peak).to_numpy()
    ngid = fr.drop_duplicates("id").set_index("id").gene_id.reindex(ng).to_numpy()
    other = lp[lp.gene.to_numpy() != ng]
    k_link = (other.peak + "|" + other.gene_id).isin(gt_key).to_numpy()
    ng_o = ngid[lp.gene.to_numpy() != ng]
    k_near = np.array([f"{p}|{g}" in gt_key for p, g in zip(other.peak, ng_o)])
    # PCHi-C (H1 / H1-derived) contact of the peak with the linked gene vs with its nearest gene, peaks where both genes
    # have captured promoters
    both = other.gene.isin(pcs_genes).to_numpy() & pd.Series(ng[lp.gene.to_numpy() != ng]).isin(pcs_genes).to_numpy()
    pl = (other.peak + "|" + other.gene).isin(pc_es).to_numpy()
    pn = pd.Series([f"{p}|{g}" for p, g in zip(other.peak, ng[lp.gene.to_numpy() != ng])]).isin(pc_es).to_numpy()
    rows.append({"celltype": ct, "distal_links_not_nearest": len(other), "n_both_captured": int(both.sum()),
                 "pchic_rate_linked_gene": pl[both].mean() if both.any() else np.nan,
                 "pchic_rate_nearest_gene": pn[both].mean() if both.any() else np.nan,
                 "gtex_rate_linked_gene": k_link.mean() if len(other) else np.nan,
                 "gtex_rate_nearest_gene": k_near.mean() if len(other) else np.nan,
                 "frac_links_nearest": np.mean(lp.gene.to_numpy() == ng) if len(lp) else np.nan})
pd.DataFrame(rows).to_csv(f"{RES}/link_vs_nearest_gene_gtex.tsv", sep="\t", index=False, float_format="%.4g")
print(pd.DataFrame(rows).to_string(), flush=True)

# ---------------- genetic validation with HDC caQTL / eQTL ----------------
rows, qv_pairs = [], []
for ct in cts:
    qf = f"{BASE}/caQTL_Output_optCov/{ct}/allchr_cis_{TAG}.cis_qtl.qval.txt.gz"
    if not os.path.exists(qf) or ct == "Unknown":
        continue
    q = pd.read_csv(qf, sep="\t", usecols=["phenotype_id", "variant_id", "pval_beta", "slope"])
    q = q[q.pval_beta < 0.05].rename(columns={"phenotype_id": "peak", "slope": "ca_slope", "variant_id": "variant_id"})
    d = P[ct].merge(q[["peak", "variant_id", "ca_slope"]], on="peak")
    d["cls"] = np.select([(d.padj < FDR) & (d.covar > 0), (d.padj < FDR) & (d.covar < 0), d.padj > 0.5],
                         ["link_pos", "link_neg", "unlinked"], "other")
    d = d[d.cls != "other"]
    ev = []
    for c, x in d.groupby("chr"):
        f = glob.glob(f"{BASE}/eQTL_Output_optCov/{ct}/{c}_cis_nominal_{TAG}.cis_qtl_pairs.*.parquet")
        if not f:
            continue
        t = pq.read_table(f[0], columns=["phenotype_id", "variant_id", "pval_nominal", "slope"],
                          filters=[("variant_id", "in", list(x.variant_id.unique()))]).to_pandas()
        ev.append(x.merge(t.rename(columns={"phenotype_id": "gene", "pval_nominal": "e_p", "slope": "e_slope"}),
                          on=["gene", "variant_id"]))
    if not ev:
        continue
    ev = pd.concat(ev).dropna(subset=["e_p"])
    ev["concordant"] = np.sign(ev.ca_slope * ev.e_slope) == np.sign(ev.test_stat)
    qv_pairs.append(ev.assign(celltype=ct))
    for cls, x in ev.groupby("cls"):
        s = x[x.e_p < 0.05]
        rows.append({"celltype": ct, "class": cls, "pairs": len(x), "frac_eqtl_p05": (x.e_p < 0.05).mean(),
                     "n_eqtl_p05": len(s), "concordance": s.concordant.mean() if len(s) else np.nan,
                     "binom_p": binomtest(int(s.concordant.sum()), len(s), 0.5).pvalue if len(s) else np.nan})
qv = pd.DataFrame(rows)
if len(qv_pairs):
    allq = pd.concat(qv_pairs)
    allq.to_csv(f"{RES}/link_qtl_pairs.tsv.gz", sep="\t", index=False, float_format="%.4g")
    for cls, x in allq.groupby("cls"):
        s = x[x.e_p < 0.05]
        qv = pd.concat([qv, pd.DataFrame([{"celltype": "pooled", "class": cls, "pairs": len(x), "frac_eqtl_p05": (x.e_p < 0.05).mean(),
                                           "n_eqtl_p05": len(s), "concordance": s.concordant.mean(),
                                           "binom_p": binomtest(int(s.concordant.sum()), len(s), 0.5).pvalue if len(s) else np.nan}])])
qv.to_csv(f"{RES}/link_qtl_validation.tsv", sep="\t", index=False, float_format="%.4g")
print(qv.to_string(), flush=True)

# ---------------- cross-cell-type sharing and pooled (AllCells) recovery ----------------
key = {ct: P[ct].assign(k=P[ct].gene + "|" + P[ct].peak).set_index("k")[["pval", "test_stat", "padj"]] for ct in P}
rows = []
named = [c for c in cts if c != "Unknown"]
for a in named:
    la = key[a][(key[a].padj < FDR) & (key[a].test_stat > 0)]
    for b in named + (["AllCells"] if "AllCells" in key else []):
        x = key[b].reindex(la.index).dropna()
        rows.append({"discovery": a, "replication": b, "n_tested": len(x),
                     "pi1": 1 - min(1, (x.pval > 0.5).mean() / 0.5) if len(x) >= 20 else np.nan,
                     "rep_rate": ((x.pval < 0.05) & (x.test_stat > 0)).mean() if len(x) else np.nan,
                     "rep_fdr": ((x.padj < FDR) & (x.test_stat > 0)).mean() if len(x) else np.nan})
sh = pd.DataFrame(rows)
sh.to_csv(f"{RES}/link_sharing.tsv", sep="\t", index=False, float_format="%.4g")
# number of cell types in which each positive link is significant; recovery by AllCells
pos = links[(links.celltype != "AllCells") & (links.celltype != "Unknown") & (links.covar > 0)].assign(k=lambda x: x.gene + "|" + x.peak)
nct = pos.groupby("k").celltype.nunique()
if "AllCells" in key:
    ac = key["AllCells"].reindex(nct.index)
    rec = pd.DataFrame({"n_ct": nct, "allcells_tested": ac.padj.notna(), "allcells_link": (ac.padj < FDR) & (ac.test_stat > 0)})
    r = rec.groupby(rec.n_ct.clip(upper=4)).agg(n=("allcells_link", "size"), tested=("allcells_tested", "mean"),
                                                 recovered=("allcells_link", "mean")).reset_index()
    r.to_csv(f"{RES}/link_allcells_recovery.tsv", sep="\t", index=False, float_format="%.4g")
    print(r.to_string(), flush=True)
    al = links[(links.celltype == "AllCells") & (links.covar > 0)].assign(k=lambda x: x.gene + "|" + x.peak)
    print("AllCells positive links:", len(al), "found in >= 1 cell type:", al.k.isin(nct.index).mean(), flush=True)

# ---------------- gene / peak specificity of links ----------------
rows = []
cpm_g = {}
for ct in named:
    x = pd.read_csv(f"{BASE}/PseudoBulk/{ct}_RNA_SumCount.tsv.gz", sep="\t", index_col=0).sum(axis=1)
    cpm_g[ct] = x / x.sum() * 1e6
cpm_g = pd.DataFrame(cpm_g).fillna(0)
for ct in named:
    other = cpm_g.drop(columns=ct).mean(axis=1)
    s = ((cpm_g[ct] + 0.1) / (other + 0.1))[cpm_g[ct] >= 1]
    gs = set(s.sort_values(ascending=False).index[: int(0.1 * len(s))])
    d = P[ct]
    tg = pd.Series(d.gene.unique())
    lg = d[(d.padj < FDR) & (d.covar > 0)].gene.unique()
    tpk = pd.Series(d.peak.unique())
    lpk = d[(d.padj < FDR) & (d.covar > 0)].peak.unique()
    spk = pki[f"spec_{ct}"]
    rows.append({"celltype": ct, "frac_linked_genes_ctspec": np.isin(lg, list(gs)).mean() if len(lg) else np.nan,
                 "frac_tested_genes_ctspec": tg.isin(gs).mean(),
                 "frac_linked_peaks_ctspec": spk.reindex(lpk).mean() if len(lpk) else np.nan,
                 "frac_tested_peaks_ctspec": spk.reindex(tpk).mean()})
gsp = pd.DataFrame(rows)
gsp.to_csv(f"{RES}/link_specificity.tsv", sep="\t", index=False, float_format="%.4g")
print(gsp.to_string(), flush=True)

# ---------------- target genes of trait-colocalized caQTL peaks ----------------
h = pd.read_csv(f"{BASE}/GWAS/results/coloc_H4gt0.5.tsv", sep="\t")
h = h[h.modality == "ATAC"].rename(columns={"phenotype_id": "peak", "celltype": "coloc_celltype"})
man = pd.read_csv(f"{BASE}/GWAS/gwas_manifest.tsv", sep="\t")[["trait_id", "label"]]
if "label" not in h.columns:
    h = h.merge(man, on="trait_id", how="left")
lp = links[(links.covar > 0) & (links.celltype != "Unknown")]
rows = []
for r in h.itertuples():
    x = lp[lp.peak == r.peak]
    same = x[x.celltype == r.coloc_celltype]
    rows.append({"trait_id": r.trait_id, "label": r.label, "coloc_celltype": r.coloc_celltype, "peak": r.peak,
                 "nearest_gene": pki.nearest_gene.get(r.peak),
                 "linked_genes_same_ct": ",".join(sorted(set(same.gene))),
                 "linked_genes_any_ct": ";".join(f"{g}({c})" for g, c in sorted(set(zip(x.gene, x.celltype)))),
                 "pchic_es_linked": ",".join(sorted(g for g in set(x.gene) if f"{r.peak}|{g}" in pc_es)),
                 "tested_in_ct": r.coloc_celltype in P and (P[r.coloc_celltype].peak == r.peak).any()})
gt_ = pd.DataFrame(rows)
gt_["PP.H4"] = h["PP.H4"].to_numpy()
gt_.to_csv(f"{RES}/link_gwas_coloc_peak_targets.tsv", sep="\t", index=False, float_format="%.4g")
u = gt_.drop_duplicates("peak")
print(f"trait-colocalized caQTL peaks: {len(u)}; tested in coloc CT {int(u.tested_in_ct.sum())}; "
      f"linked in coloc CT {int((u.linked_genes_same_ct != '').sum())}; linked in any CT {int((u.linked_genes_any_ct != '').sum())}", flush=True)
print(gt_[gt_.linked_genes_any_ct != ""].to_string(), flush=True)
print("done")
