"""
Audit of the coloc summary (08_summarize_coloc_smr.py) before it goes into the paper.

  1. integrity of coloc_all: PP.H0-H4 sum to 1, no duplicate feature x trait rows, H4>0.5 table = filter of coloc_all
  2. counts: pairs, features, trait-loci per cell type / modality; pooled (AllCells) vs per-CT overlap at three levels
       trait-locus  same trait, best SNPs < 500 kb apart
       phenotype-locus  same trait after merging the same phenotype measured in several cohorts
                        (UKB_EUR / UKB_meta_hq / FinnGen / BBJ / PAGE labels mapped to one name, PHENO below)
       region       any trait, best SNPs < 500 kb apart (independent genomic regions)
  3. prior sensitivity: PP.H4 recomputed analytically for p12 = 5e-6 (other priors unchanged)
  4. QTL strength at hits: feature permutation p (pval_beta) and min nominal p
  5. exports every PP.H4 > 0.5 pair for an independent rerun with R coloc::coloc.abf (07v_validate_coloc.R)
Output: GWAS/results/coloc_audit/{hits_audited.tsv, loci.tsv, overlap.tsv, summary.txt, R/pair*.tsv, R/meta.tsv}
"""
import os
import re

import numpy as np
import pandas as pd
from scipy.stats import norm

BASE = os.environ["HDC_DIR"] + "/"
R = f"{BASE}/GWAS/results/"
OUT = f"{R}/coloc_audit/"
os.makedirs(f"{OUT}/R", exist_ok=True)
TAG = "noBootstrap_TMM_invnorm"
log = open(f"{OUT}/summary.txt", "w")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    log.write(s + "\n")


# phenotype harmonization across cohorts (same trait measured in several GWAS)
PHENO = [(r"body.mass|^bmi$", "BMI"), (r"^height$|comparative height", "height"), (r"platelet|^plt$", "platelet count"),
         (r"haemoglobin conc|hemoglobin conc|^hb$", "haemoglobin"), (r"^mchc$|corpuscular hemoglobin conc", "MCHC"),
         (r"ldl", "LDL"), (r"hdl", "HDL"), (r"^cholesterol|total cholesterol|^tc$", "total cholesterol"),
         (r"triglycer|^tg$", "triglycerides"), (r"c-reactive|^crp$", "CRP"), (r"systolic|^sbp$", "SBP"),
         (r"diastolic|^dbp$", "DBP"), (r"hba1c|glycated", "HbA1c"), (r"egfr|glomerular", "eGFR"),
         (r"white blood|leukocyte|^wbc$", "WBC"), (r"type 2 diabetes|non-insulin|^t2d$", "T2D"),
         (r"hypertension", "hypertension"), (r"obesity", "obesity"),
         (r"ischaemic heart|coronary|angina", "CAD")]


def pheno(label):
    l = str(label).lower()
    for pat, name in PHENO:
        if re.search(pat, l):
            return name
    return l


def cluster(pos, gap=500_000):
    o = np.argsort(pos)
    cid = np.empty(len(pos), int)
    cid[o] = np.concatenate([[0], np.cumsum(np.diff(pos[o]) > gap)])
    return cid


c = pd.read_csv(f"{R}/coloc_all.tsv.gz", sep="\t")
pp = c[[f"PP.H{i}" for i in range(5)]].to_numpy()
say("coloc_all rows", len(c), "| max |sum PP - 1|", np.nanmax(np.abs(pp.sum(1) - 1)).round(6),
    "| NaN PP rows", int(np.isnan(pp).any(1).sum()))
say("duplicate celltype x modality x feature x trait rows:",
    int(c.duplicated(["celltype", "modality", "phenotype_id", "trait_id"]).sum()))
h = c[c["PP.H4"] > 0.5].copy()
old = pd.read_csv(f"{R}/coloc_H4gt0.5.tsv", sep="\t")
say("PP.H4>0.5 pairs recomputed", len(h), "| in coloc_H4gt0.5.tsv", len(old))
summ = pd.read_csv(f"{R}/coloc_summary.tsv", sep="\t")
say("coloc_summary n_H4>0.5 total", int(summ["n_H4>0.5"].sum()), "| n_tests total", int(summ.n_tests.sum()),
    "vs coloc_all rows", len(c))

# --- prior sensitivity (p12 = 5e-6: H4 prior halves, H0-H3 unchanged)
k = 0.5
h["PP.H4_p12_5e-6"] = h["PP.H4"] * k / (h[["PP.H0", "PP.H1", "PP.H2", "PP.H3"]].sum(1) + h["PP.H4"] * k)
c["PP.H4_p12_5e-6"] = c["PP.H4"] * k / (c[["PP.H0", "PP.H1", "PP.H2", "PP.H3"]].sum(1) + c["PP.H4"] * k)

# --- QTL strength: feature-level permutation p
qv = []
for mod, d in (("RNA", "eQTL"), ("ATAC", "caQTL")):
    for ct in h[h.modality == mod].celltype.unique():
        q = pd.read_csv(f"{BASE}/{d}_Output_optCov/{ct}/allchr_cis_{TAG}.cis_qtl.qval.txt.gz", sep="\t",
                        usecols=["phenotype_id", "pval_beta", "qval"])
        qv.append(q.assign(celltype=ct, modality=mod))
h = h.merge(pd.concat(qv), on=["celltype", "modality", "phenotype_id"], how="left")
h["pheno"] = h.label.map(pheno)
h["design"] = np.where(h.celltype == "AllCells", "pooled", "per-CT")

# --- loci
rows = []
for design, d in [("pooled", h[h.design == "pooled"]), ("per-CT", h[h.design == "per-CT"]), ("all", h)]:
    for level, keys in [("trait-locus", ["trait_id", "chr"]), ("phenotype-locus", ["pheno", "chr"]), ("region", ["chr"])]:
        n = sum(cluster(g.best_pos38.to_numpy()).max() + 1 for _, g in d.groupby(keys))
        rows.append({"design": design, "level": level, "n_loci": n})
loci = pd.DataFrame(rows).pivot(index="level", columns="design", values="n_loci")
# shared: cluster pooled + per-CT jointly; a locus is shared if it holds hits of both designs
ov = []
for level, keys in [("trait-locus", ["trait_id", "chr"]), ("phenotype-locus", ["pheno", "chr"]), ("region", ["chr"])]:
    sh = po = pc = 0
    for _, g in h.groupby(keys):
        cid = cluster(g.best_pos38.to_numpy())
        for i in np.unique(cid):
            dz = set(g.design.to_numpy()[cid == i])
            sh += dz == {"pooled", "per-CT"}
            po += dz == {"pooled"}
            pc += dz == {"per-CT"}
    ov.append({"level": level, "pooled_only": po, "perCT_only": pc, "shared": sh, "total": po + pc + sh})
ov = pd.DataFrame(ov)
loci.to_csv(f"{OUT}/loci.tsv", sep="\t")
ov.to_csv(f"{OUT}/overlap.tsv", sep="\t", index=False)
say("\nloci (joint clustering, < 500 kb):\n" + ov.to_string(index=False))
say("\npairs by design x modality:\n" + h.groupby(["design", "modality"]).size().to_string())
say("\nunique features by design x modality:\n" + h.groupby(["design", "modality"]).phenotype_id.nunique().to_string())
say("\nPP.H4 > 0.8:", int((h["PP.H4"] > 0.8).sum()), "| per design:", h[h["PP.H4"] > 0.8].design.value_counts().to_dict())
s = h["PP.H4_p12_5e-6"] > 0.5
say("survive p12=5e-6 (PP.H4>0.5): all", int(s.sum()), "of", len(h), "| per-CT", int(s[h.design == "per-CT"].sum()),
    "of", int((h.design == "per-CT").sum()), "| pooled", int(s[h.design == "pooled"].sum()), "of", int((h.design == "pooled").sum()))
say("QTL feature permutation p (pval_beta) at hits: <0.01", int((h.pval_beta < 0.01).sum()), "| <0.05",
    int((h.pval_beta < 0.05).sum()), "| >=0.05", int((h.pval_beta >= 0.05).sum()), "| median", round(h.pval_beta.median(), 4))
say("min nominal QTL p at hits: median", f"{h.qtl_minp.median():.2g}", "| >1e-3:", int((h.qtl_minp > 1e-3).sum()))
say("GWAS min p at hits: median", f"{h.gwas_minp.median():.2g}", "| <5e-8:", int((h.gwas_minp < 5e-8).sum()))
say("best SNP SNP.PP.H4 median", round(h.best_snp_ppH4.median(), 3))
h.sort_values("PP.H4", ascending=False).to_csv(f"{OUT}/hits_audited.tsv", sep="\t", index=False, float_format="%.4g")

# --- export all hits for R coloc
man = pd.read_csv(f"{BASE}/GWAS/gwas_manifest.tsv", sep="\t").set_index("trait_id")
i = 0
meta = []
for (ct, mod, chrom), g0 in h.groupby(["celltype", "modality", "chr"]):
    qdir = f"{BASE}/{'eQTL' if mod == 'RNA' else 'caQTL'}_Output_optCov/{ct}/"
    u = pd.read_parquet(f"{BASE}/GWAS/variant_universe.parquet", columns=["vidx", "id"], filters=[("chr", "==", chrom)])
    nom = pd.read_parquet(f"{qdir}/chr{chrom}_cis_nominal_{TAG}.cis_qtl_pairs.{chrom}.parquet",
                          columns=["phenotype_id", "variant_id", "pval_nominal", "slope"],
                          filters=[("phenotype_id", "in", list(g0.phenotype_id.unique()))]).merge(u, left_on="variant_id", right_on="id")
    for _, r in g0.iterrows():
        g = pd.read_parquet(f"{BASE}/GWAS/harmonized/{r.trait_id}.parquet", filters=[("chr", "==", chrom)])
        qc = pd.read_csv(f"{BASE}/GWAS/harmonized/qc/{r.trait_id}.tsv", sep="\t").iloc[0]
        d = nom[nom.phenotype_id == r.phenotype_id].merge(g, on="vidx")
        d = d[np.isfinite(d.slope) & (d.slope != 0)]
        z = norm.isf(np.clip(d.pval_nominal, 1e-300, 1) / 2) * np.sign(d.slope)
        d = d.assign(varbeta_qtl=(d.slope / z) ** 2, varbeta_gwas=d.se.astype(float) ** 2)
        d[["variant_id", "slope", "varbeta_qtl", "beta", "varbeta_gwas"]].to_csv(f"{OUT}/R/pair{i}.tsv", sep="\t", index=False)
        sdy = qc.sdY if man.loc[r.trait_id, "type"] == "quant" and np.isfinite(qc.sdY) else 1.0
        meta.append(dict(pair=i, trait=r.trait_id, feature=f"{ct}|{r.phenotype_id}", type=man.loc[r.trait_id, "type"],
                         sdY=sdy, nsnps=r.nsnps, H3=r["PP.H3"], H4=r["PP.H4"]))
        i += 1
pd.DataFrame(meta).to_csv(f"{OUT}/R/meta.tsv", sep="\t", index=False)
say("exported", i, "pairs for R coloc")
