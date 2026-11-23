"""
Feature sets for S-LDSC annotations (hg19 intervals from 15a_liftover_features.R).

batch "spec" (donor-number independent; pseudobulk over all donors, all 11 cell types):
  HDC_peaks              all HDC multiome peaks (control annotation for the peak sets)
  peakspec_<CT>          cell-type-specific peaks (same definition as 11_celltype_peak_enrichment.py)
  HDC_genes              genes with CPM >= 1 in at least one cell type, +/-100 kb (control for gene sets)
  genespec_<CT>          top 10% most cell-type-specific genes among those with CPM >= 1 in the cell type, +/-100 kb
                         (Finucane et al. 2018 style)
batch "qtl" (QTL-derived; peak / gene level, so the annotation does not depend on YRI lead-SNP LD):
  caQTLpeak_<CT>         peaks with mashr lfsr < 0.05 in the cell type (14_mash.R)
  eGene_<CT>             genes with mashr lfsr < 0.05 in the cell type, +/-100 kb
  caQTLpeak_AllCells / eGene_AllCells   pooled pseudobulk map, permutation pval_beta < 0.01
  caQTLpeak_anyCT / eGene_anyCT         lfsr < 0.05 in any cell type
  control annotations for this batch: HDC_peaks_tested / HDC_genes_tested (features tested in any map, +/-100 kb genes)
batch "link" (scMultiMap peak-gene links, PeakAnnot/results/links_fdr0.1.tsv.gz from 26_link_analysis.py; positive links):
  linked_<CT>            peaks with >= 1 positive link in the cell type           (model: baselineLD + HDC_peaks + set)
  speclinked_<CT>        cell-type-specific peaks (peakspec_<CT>) with a positive link in that cell type
                         (model: baselineLD + HDC_peaks + peakspec_<CT> + set: h2 of linked beyond CT-specific peaks)
  LD scores of HDC_peaks / peakspec_<CT> are reused from batch "spec".

Output: GWAS/ldsc/sets_<batch>.tsv.gz (annot, chr, start, end)   [hg19, chr without "chr"]
Usage: python 15_ldsc_sets.py <spec|qtl|link>
"""
import os
import glob
import sys

import numpy as np
import pandas as pd

BASE = os.environ["HDC_DIR"] + "/"
L = f"{BASE}/GWAS/ldsc/"
WIN = 100_000
batch = sys.argv[1]
fa = pd.read_csv(f"{L}/features_hg19_ATAC.tsv.gz", sep="\t", dtype={"chr": str}).set_index("id")
fr = pd.read_csv(f"{L}/features_hg19_RNA.tsv.gz", sep="\t", dtype={"chr": str}).set_index("id")
sets = []


def add(name, ids, gene=False):
    f = fr if gene else fa
    d = f.loc[f.index.intersection(list(ids))].copy()
    if gene:
        d["start"] = np.maximum(d.start - WIN, 1)
        d["end"] = d.end + WIN
    sets.append(d[["chr", "start", "end"]].assign(annot=name))
    print(name, len(d), flush=True)


cts = sorted(f.split("/")[-1].replace("_ATAC_SumCount.tsv.gz", "")
             for f in glob.glob(f"{BASE}/PseudoBulk/*_ATAC_SumCount.tsv.gz")
             if not f.split("/")[-1].startswith("AllCells"))

if batch == "spec":
    add("HDC_peaks", fa.index)
    sp = pd.read_csv(f"{BASE}/GWAS/celltype_specific_peaks.tsv.gz", sep="\t")
    for ct in cts:
        add(f"peakspec_{ct}", sp.peak[sp[ct]])
    cpm = {}
    for ct in cts:
        x = pd.read_csv(f"{BASE}/PseudoBulk/{ct}_RNA_SumCount.tsv.gz", sep="\t", index_col=0).sum(axis=1)
        cpm[ct] = x / x.sum() * 1e6
    cpm = pd.DataFrame(cpm).fillna(0)
    add("HDC_genes", cpm.index[(cpm >= 1).any(axis=1)], gene=True)
    for ct in cts:
        other = cpm.drop(columns=ct).mean(axis=1)
        s = ((cpm[ct] + 0.1) / (other + 0.1))[cpm[ct] >= 1]
        add(f"genespec_{ct}", s.sort_values(ascending=False).index[: int(0.1 * len(s))], gene=True)

elif batch == "qtl":
    tested = {"ATAC": set(), "RNA": set()}
    for mod in ("ATAC", "RNA"):
        gene = mod == "RNA"
        tag = "eGene" if gene else "caQTLpeak"
        m = pd.read_csv(f"{BASE}/mash/{mod}_mash_strong.tsv.gz", sep="\t")
        lf = m.filter(like="lfsr_")
        tested[mod] |= set(m.phenotype_id)
        for c in lf.columns:
            add(f"{tag}_{c[5:]}", m.phenotype_id[lf[c] < 0.05], gene=gene)
        add(f"{tag}_anyCT", m.phenotype_id[lf.min(axis=1) < 0.05], gene=gene)
        qd = "eQTL" if gene else "caQTL"
        a = pd.read_csv(f"{BASE}/{qd}_Output_optCov/AllCells/allchr_cis_noBootstrap_TMM_invnorm.cis_qtl.qval.txt.gz",
                        sep="\t", usecols=["phenotype_id", "pval_beta"])
        tested[mod] |= set(a.phenotype_id)
        add(f"{tag}_AllCells", a.phenotype_id[a.pval_beta < 0.01], gene=gene)
    add("HDC_peaks_tested", tested["ATAC"])
    add("HDC_genes_tested", tested["RNA"], gene=True)

elif batch == "link":
    lk = pd.read_csv(f"{BASE}/PeakAnnot/results/links_fdr0.1.tsv.gz", sep="\t", usecols=["celltype", "peak", "covar"])
    lk = lk[(lk.covar > 0) & ~lk.celltype.isin(["AllCells", "Unknown"])]
    sp = pd.read_csv(f"{BASE}/GWAS/celltype_specific_peaks.tsv.gz", sep="\t")
    for ct, d in lk.groupby("celltype"):
        if d.peak.nunique() < 200:
            continue
        add(f"linked_{ct}", d.peak.unique())
        add(f"speclinked_{ct}", set(d.peak) & set(sp.peak[sp[ct]]))

pd.concat(sets)[["annot", "chr", "start", "end"]].to_csv(f"{L}/sets_{batch}.tsv.gz", sep="\t", index=False)
print("done")
