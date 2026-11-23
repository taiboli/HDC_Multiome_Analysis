"""
Multiome follow-up and novelty annotation of trait colocalizations (08_summarize_coloc_smr.py, PP.H4 > 0.5).

1. Peak <-> gene links: for each trait-colocalized feature, coloc.abf between that feature and every
   other-modality feature of the same cell type whose cis window contains the trait coloc best SNP
   (caQTL x eQTL, both n = 17 donors; same z-from-p / varbeta conventions as 07_coloc_smr.py).
   -> trait - peak - gene chains when the QTL-QTL coloc PP.H4 > 0.5.
2. GTEx v8 (49 tissues; $GTEX_DIR, read-only):
   - eGene in any GTEx tissue (q < 0.05)?  (RNA features; gene symbol match)
   - is the trait-coloc best SNP (or the linked gene's) a significant GTEx eQTL variant for any gene / for that gene?
   Loci whose best SNP is not a GTEx eQTL for any gene in any tissue are "not explained by adult-tissue eQTLs".

Outputs: GWAS/results/coloc_H4gt0.5_annotated.tsv, GWAS/results/peak_gene_links.tsv
"""
import os
import glob

import numpy as np
import pandas as pd
from scipy.stats import norm

BASE = os.environ["HDC_DIR"] + "/"
TAG = "noBootstrap_TMM_invnorm"
R = f"{BASE}/GWAS/results/"
GTEX = os.environ["GTEX_DIR"] + "/"
P1, P2, P12 = 1e-4, 1e-4, 1e-5


def labf(z, V, W=0.15 ** 2):
    r = W / (W + V)
    return 0.5 * (np.log1p(-r) + r * z * z)


def lse(x):
    m = np.max(x)
    return m + np.log(np.sum(np.exp(x - m)))


def coloc_qq(a, b):
    """coloc.abf between two QTL nominal tables (variant_id, pval_nominal, slope)"""
    d = a.merge(b, on="variant_id")
    d = d[(d.slope_x != 0) & (d.slope_y != 0)]
    if len(d) < 50:
        return None
    out = []
    for s in ("x", "y"):
        z = norm.isf(np.clip(d[f"pval_nominal_{s}"], 1e-300, 1) / 2) * np.sign(d[f"slope_{s}"])
        out.append(labf(z, (d[f"slope_{s}"] / z) ** 2))
    l1, l2 = out
    s1, s2, s12 = lse(l1), lse(l2), lse(l1 + l2)
    with np.errstate(divide="ignore"):
        h3 = np.log(P1) + np.log(P2) + s1 + s2 + np.log1p(-np.exp(s12 - s1 - s2))
    h = np.array([0, np.log(P1) + s1, np.log(P2) + s2, h3, np.log(P12) + s12])
    pp = np.exp(h - h.max())
    pp /= pp.sum()
    best = d.variant_id.to_numpy()[np.argmax(l1 + l2)]
    return dict(nsnps_qq=len(d), qq_PP_H3=pp[3], qq_PP_H4=pp[4], qq_best_snp=best)


hits = pd.read_csv(f"{R}/coloc_H4gt0.5.tsv", sep="\t")
if len(hits) == 0:
    raise SystemExit("no PP.H4 > 0.5 colocalizations; nothing to annotate")
u = pd.read_parquet(f"{BASE}/GWAS/variant_universe.parquet", columns=["chr", "pos38", "id", "ref", "alt"])
u = u[u.id.isin(set(hits.best_snp))]

# ---------- 1. peak <-> gene links
links = []
for (ct, mod, chrom), h in hits.groupby(["celltype", "modality", "chr"]):
    other = "ATAC" if mod == "RNA" else "RNA"
    od = f"{BASE}/{'eQTL' if other == 'RNA' else 'caQTL'}_Output_optCov/{ct}/"
    sd = f"{BASE}/{'eQTL' if mod == 'RNA' else 'caQTL'}_Output_optCov/{ct}/"
    fo = f"{od}/chr{chrom}_cis_nominal_{TAG}.cis_qtl_pairs.{chrom}.parquet"
    if not glob.glob(fo):
        continue
    # other-modality features whose window contains a best SNP
    cand = pd.read_parquet(fo, columns=["phenotype_id", "variant_id"],
                           filters=[("variant_id", "in", list(h.best_snp.unique()))])
    if len(cand) == 0:
        continue
    feats = list(cand.phenotype_id.unique())
    oth = pd.read_parquet(fo, columns=["phenotype_id", "variant_id", "pval_nominal", "slope"],
                          filters=[("phenotype_id", "in", feats)])
    slf = pd.read_parquet(f"{sd}/chr{chrom}_cis_nominal_{TAG}.cis_qtl_pairs.{chrom}.parquet",
                          columns=["phenotype_id", "variant_id", "pval_nominal", "slope"],
                          filters=[("phenotype_id", "in", list(h.phenotype_id.unique()))])
    for f, snps in h.groupby("phenotype_id").best_snp:
        a = slf[slf.phenotype_id == f].drop(columns="phenotype_id")
        for g in cand[cand.variant_id.isin(set(snps))].phenotype_id.unique():
            r = coloc_qq(a, oth[oth.phenotype_id == g].drop(columns="phenotype_id"))
            if r is None:
                continue
            links.append(dict(celltype=ct, modality=mod, chr=chrom, phenotype_id=f, linked_modality=other,
                              linked_feature=g, **r))
    print(ct, mod, chrom, flush=True)
links = pd.DataFrame(links)
links.to_csv(f"{R}/peak_gene_links.tsv", sep="\t", index=False, float_format="%.4g")

# ---------- 2. GTEx v8
egenes, pairs = set(), []
id2name = {}
for f in sorted(glob.glob(f"{GTEX}/*.v8.egenes.txt.gz")):
    e = pd.read_csv(f, sep="\t", usecols=["gene_id", "gene_name", "qval"])
    id2name.update(dict(zip(e.gene_id, e.gene_name)))
    egenes |= set(e.gene_name[e.qval < 0.05])
key = set(u.chr.astype(str) + "_" + u.pos38.astype(str))
for f in sorted(glob.glob(f"{GTEX}/*.v8.signif_variant_gene_pairs.txt.gz")):
    tissue = f.split("/")[-1].split(".v8")[0]
    p = pd.read_csv(f, sep="\t", usecols=["variant_id", "gene_id"])
    k = p.variant_id.str.replace("chr", "").str.split("_").str[:2].str.join("_")
    p = p[k.isin(key)].assign(key=k[k.isin(key)], tissue=tissue)
    pairs.append(p)
pairs = pd.concat(pairs)
pairs["gene_name"] = pairs.gene_id.map(id2name)
u["key"] = u.chr.astype(str) + "_" + u.pos38.astype(str)
snp2key = dict(zip(u.id, u.key))
by_snp = pairs.groupby("key").agg(gtex_eqtl_genes=("gene_name", lambda x: ",".join(sorted(set(map(str, x))))),
                                  gtex_n_tissues=("tissue", "nunique"))

hits["gtex_key"] = hits.best_snp.map(snp2key)
hits = hits.merge(by_snp, left_on="gtex_key", right_index=True, how="left")
hits["best_snp_is_GTEx_eQTL"] = hits.gtex_eqtl_genes.notna()
is_rna = hits.modality == "RNA"
hits["gene_is_GTEx_eGene"] = np.where(is_rna, hits.phenotype_id.isin(egenes), np.nan)
hits["best_snp_is_GTEx_eQTL_for_this_gene"] = [
    (isinstance(g, str) and f in g.split(",")) if r else np.nan
    for f, g, r in zip(hits.phenotype_id, hits.gtex_eqtl_genes, is_rna)]
if len(links):
    l5 = links[links.qq_PP_H4 > 0.5]
    lk = l5.groupby(["celltype", "modality", "phenotype_id"]).linked_feature.agg(lambda x: ",".join(sorted(set(x))))
    hits = hits.merge(lk.rename("linked_features_qqH4gt0.5"), left_on=["celltype", "modality", "phenotype_id"],
                      right_index=True, how="left")
hits.drop(columns="gtex_key").to_csv(f"{R}/coloc_H4gt0.5_annotated.tsv", sep="\t", index=False, float_format="%.4g")
print(hits.groupby("modality")[["best_snp_is_GTEx_eQTL"]].mean())
print("RNA hits gene not GTEx eGene:", int((is_rna & ~hits.phenotype_id.isin(egenes)).sum()), "/", int(is_rna.sum()))
