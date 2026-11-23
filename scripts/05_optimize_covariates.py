"""
Covariate / parameter optimization for QTL discovery power on a subset of chromosomes.

For one cell type, for RNA and ATAC (existing inputs in eQTL/<CT>/, caQTL/<CT>/):
  Stage A1: PhenoPC {0..6,8,10} with no other covariates
  Stage A2: with best PhenoPC count, Sex {0,1} x GenoPC {0,1,2,3,5}
           (maf 0, window 1 Mb; skip if residual dof = n - 2 - n_cov < 2)
  Stage B: with best Stage-A covariates, grid over maf_threshold x cis-window
Permutation mapping only (nperm=1000, fixed seed so all configs share permutations),
q-values computed across the optimization chromosomes.

PhenoPCs reproduce R princomp(features x donors)$loadings (eigenvectors of the donor covariance).
Sex / GenoPCs are taken from the existing covariate file.

Output: CovOptimization/<CT>_<MOD>/summary.tsv and per-config permutation results.
Usage: python 05_optimize_covariates.py <CellType> <RNA,ATAC> [chroms, default 2,10,17,20]
"""
import os
import sys
import time

import numpy as np
import pandas as pd
from tensorqtl import cis, genotypeio, read_phenotype_bed, rfunc
from tensorqtl.core import SimpleLogger

GENO = os.environ["GENOTYPE_DIR"] + "/Genotype/plink/chr{}"
BASE = os.environ["HDC_DIR"] + "/"
TAG = "noBootstrap_TMM_invnorm"
NPERM, SEED = 1000, 2025

celltype = sys.argv[1]
modalities = sys.argv[2].split(",")
chroms = sys.argv[3].split(",") if len(sys.argv) > 3 else ["2", "10", "17", "20"]


def pheno_pcs(phenotype_df, k):
    x = phenotype_df.values
    w, v = np.linalg.eigh(np.cov(x, rowvar=False))
    v = v[:, np.argsort(w)[::-1][:k]]
    return pd.DataFrame(v, index=phenotype_df.columns, columns=[f"PhenoPC_{i + 1}" for i in range(k)])


def run_config(geno, p_df, p_pos_df, cov_df, maf, window):
    res = []
    for chrom, (genotype_df, variant_df) in geno.items():
        sel = p_pos_df["chr"] == chrom
        if sel.sum() == 0:
            continue
        res.append(cis.map_cis(genotype_df, variant_df, p_df.loc[sel], p_pos_df.loc[sel],
                               covariates_df=cov_df, nperm=NPERM, window=window, maf_threshold=maf,
                               seed=SEED, verbose=False, warn_monomorphic=False,
                               logger=SimpleLogger(verbose=False)))
    res = pd.concat(res)
    res = res[res["pval_beta"].notna()]
    q, pi0 = rfunc.qvalue(res["pval_beta"])
    res["qval"] = q
    return res, {"nTested": len(res), "pi1": 1 - pi0,
                 "nSigFDR05": int((res["qval"] <= 0.05).sum()), "nSigFDR10": int((res["qval"] <= 0.10).sum())}


for modality in modalities:
    in_sub, bed, cov = ("eQTL", "geneExpr", "geneCov") if modality == "RNA" else ("caQTL", "atacExpr", "atacCov")
    pheno_file = f"{BASE}/{in_sub}/{celltype}/tensorqtl.{bed}_{TAG}.bed.gz"
    if not os.path.exists(pheno_file):
        continue
    out_dir = f"{BASE}/CovOptimization/{celltype}_{modality}/"
    os.makedirs(out_dir, exist_ok=True)
    t0 = time.time()

    phenotype_df, phenotype_pos_df = read_phenotype_bed(pheno_file)
    base_cov = pd.read_csv(f"{BASE}/{in_sub}/{celltype}/tensorqtl.{cov}_{TAG}.txt", sep="\t", index_col=0).T
    n = phenotype_df.shape[1]
    pcs = pheno_pcs(phenotype_df, 10)

    sel = phenotype_pos_df["chr"].isin(chroms)
    p_df, p_pos_df = phenotype_df[sel], phenotype_pos_df[sel]

    geno = {}
    for chrom in chroms:
        g, v = genotypeio.load_genotypes(GENO.format(chrom), select_samples=phenotype_df.columns)
        poly = g.nunique(axis=1) > 1
        geno[chrom] = (g[poly], v[poly.values])
    print(f"{celltype} {modality}: n={n}, {sel.sum()} phenotypes on chr{','.join(chroms)}", flush=True)

    rows = []

    def record(stage, sex, ngeno, npheno, maf, window):
        cols = (["Sex"] if sex else []) + [f"GenoPC_{i + 1}" for i in range(ngeno)]
        cov_df = pd.concat([base_cov[cols], pcs.iloc[:, :npheno]], axis=1)
        cov_df = cov_df if cov_df.shape[1] > 0 else None
        res, stats = run_config(geno, p_df, p_pos_df, cov_df, maf, window)
        name = f"sex{sex}_geno{ngeno}_pheno{npheno}_maf{maf}_win{window}"
        res.to_csv(f"{out_dir}/{name}.cis_qtl.txt.gz", sep="\t", float_format="%.6g")
        row = {"Stage": stage, "Sex": sex, "nGenoPC": ngeno, "nPhenoPC": npheno, "maf": maf, "window": window,
               "nCov": len(cols) + npheno, "dof": n - 2 - len(cols) - npheno, **stats}
        rows.append(row)
        pd.DataFrame(rows).to_csv(f"{out_dir}/summary.tsv", sep="\t", index=False)
        print(f"  {name}: {stats} ({time.time() - t0:.0f}s)", flush=True)

    def best():  # most FDR 10% hits, then FDR 5%, then pi1
        return pd.DataFrame(rows).sort_values(["nSigFDR10", "nSigFDR05", "pi1"], ascending=False).iloc[0]

    # Stage A1: number of PhenoPCs
    for npheno in [0, 1, 2, 3, 4, 5, 6, 8, 10]:
        if n - 2 - npheno >= 2:
            record("A1", 0, 0, npheno, 0, 1000000)
    k = int(best().nPhenoPC)

    # Stage A2: Sex and GenoPCs on top of the best PhenoPC count
    for sex in [0, 1]:
        for ngeno in [0, 1, 2, 3, 5]:
            if (sex, ngeno) != (0, 0) and n - 2 - sex - ngeno - k >= 2:
                record("A2", sex, ngeno, k, 0, 1000000)

    # Stage B: maf x window with best covariates
    a = best()
    windows = [1000000, 500000, 250000, 100000] if modality == "RNA" else [1000000, 250000, 100000, 50000, 25000]
    for maf in [0, 0.05, 0.1, 0.15, 0.2]:
        for window in windows:
            if maf == 0 and window == 1000000:
                continue
            record("B", int(a.Sex), int(a.nGenoPC), int(a.nPhenoPC), maf, window)
    print(f"{celltype} {modality} done in {time.time() - t0:.0f}s", flush=True)
