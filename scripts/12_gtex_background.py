"""
Background for the GTEx novelty annotation of 10_multiome_link_gtex.py.

The fraction of HDC trait-colocalization SNPs that are GTEx v8 eQTLs is only interpretable against how often a
GWAS-significant SNP near an HDC feature is a GTEx eQTL anyway. Background = coloc tests with PP.H4 < 0.1 and
GWAS p < 5e-8 in the window (best SNP then mostly tracks the GWAS signal), one row per (cell type, modality,
feature, best SNP), sampled up to 20,000. Same GTEx definitions as 10_multiome_link_gtex.py.

Output: GWAS/results/gtex_novelty_vs_background.tsv
"""
import os
import glob

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

BASE = os.environ["HDC_DIR"] + "/"
R = f"{BASE}/GWAS/results/"
GTEX = os.environ["GTEX_DIR"] + "/"

hits = pd.read_csv(f"{R}/coloc_H4gt0.5.tsv", sep="\t")
cols = ["celltype", "modality", "phenotype_id", "best_snp", "PP.H4", "gwas_minp"]
bg = pd.concat([d[(d["PP.H4"] < 0.1) & (d.gwas_minp < 5e-8)]
                for d in pd.read_csv(f"{R}/coloc_all.tsv.gz", sep="\t", usecols=cols, chunksize=5_000_000)])
key = ["celltype", "modality", "phenotype_id", "best_snp"]
bg = bg.drop_duplicates(key)
bg = bg.sample(min(len(bg), 20_000), random_state=1)
hits = hits.drop_duplicates(key)
sets = pd.concat([hits[key].assign(set="coloc_H4gt0.5"), bg[key].assign(set="background")], ignore_index=True)

u = pd.read_parquet(f"{BASE}/GWAS/variant_universe.parquet", columns=["chr", "pos38", "id"])
u = u[u.id.isin(set(sets.best_snp))]
snp2key = dict(zip(u.id, u.chr.astype(str) + "_" + u.pos38.astype(str)))
sets["key"] = sets.best_snp.map(snp2key)

id2name, egenes, pairs = {}, set(), []
for f in sorted(glob.glob(f"{GTEX}/*.v8.egenes.txt.gz")):
    e = pd.read_csv(f, sep="\t", usecols=["gene_id", "gene_name", "qval"])
    id2name.update(dict(zip(e.gene_id, e.gene_name)))
    egenes |= set(e.gene_name[e.qval < 0.05])
keys = set(sets.key.dropna())
for f in sorted(glob.glob(f"{GTEX}/*.v8.signif_variant_gene_pairs.txt.gz")):
    p = pd.read_csv(f, sep="\t", usecols=["variant_id", "gene_id"])
    k = p.variant_id.str.replace("chr", "").str.split("_").str[:2].str.join("_")
    pairs.append(pd.DataFrame({"key": k[k.isin(keys)], "gene": p.gene_id[k.isin(keys)].map(id2name)}))
pairs = pd.concat(pairs).drop_duplicates()
is_eqtl = set(pairs.key)
pair_set = set(zip(pairs.key, pairs.gene))

sets["snp_GTEx_eQTL"] = sets.key.isin(is_eqtl)
rna = sets.modality == "RNA"
sets["gene_GTEx_eGene"] = np.where(rna, sets.phenotype_id.isin(egenes), np.nan)
sets["snp_GTEx_eQTL_same_gene"] = np.where(rna, [(k, g) in pair_set for k, g in zip(sets.key, sets.phenotype_id)],
                                           np.nan)

rows = []
for mod in ("RNA", "ATAC"):
    for col in ["snp_GTEx_eQTL", "gene_GTEx_eGene", "snp_GTEx_eQTL_same_gene"]:
        if mod == "ATAC" and col != "snp_GTEx_eQTL":
            continue
        d = sets[sets.modality == mod]
        a = d[d.set == "coloc_H4gt0.5"][col].astype(bool)
        b = d[d.set == "background"][col].astype(bool)
        orr, p = fisher_exact([[a.sum(), (~a).sum()], [b.sum(), (~b).sum()]])
        rows.append(dict(modality=mod, annotation=col, n_coloc=len(a), frac_coloc=a.mean(), n_background=len(b),
                         frac_background=b.mean(), OR=orr, p_fisher=p))
res = pd.DataFrame(rows)
res.to_csv(f"{R}/gtex_novelty_vs_background.tsv", sep="\t", index=False, float_format="%.4g")
print(res.to_string())
