"""
Write the chosen optimized covariate file (tensorqtl.{gene,atac}Cov_optCov.txt) for one cell type x modality,
using the same PhenoPC construction as 05_optimize_covariates.py (eigenvectors of the donor covariance,
= R princomp loadings) plus optional Sex / GenoPC rows from the base covariate file.
Usage: python 05b_write_optcov.py <CellType> <RNA|ATAC> <nPhenoPC> <Sex 0|1> <nGenoPC>
"""
import os
import sys

import numpy as np
import pandas as pd
from tensorqtl import read_phenotype_bed

BASE = os.environ["HDC_DIR"] + "/"
TAG = "noBootstrap_TMM_invnorm"
ct, mod, npheno, sex, ngeno = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
sub, bed, cov = ("eQTL", "geneExpr", "geneCov") if mod == "RNA" else ("caQTL", "atacExpr", "atacCov")
p, _ = read_phenotype_bed(f"{BASE}/{sub}/{ct}/tensorqtl.{bed}_{TAG}.bed.gz")
w, v = np.linalg.eigh(np.cov(p.values, rowvar=False))
v = v[:, np.argsort(w)[::-1][:npheno]]
pcs = pd.DataFrame(v.T, columns=p.columns, index=[f"PhenoPC_{i + 1}" for i in range(npheno)])
base = pd.read_csv(f"{BASE}/{sub}/{ct}/tensorqtl.{cov}_{TAG}.txt", sep="\t", index_col=0)[p.columns]
rows = (["Sex"] if sex else []) + [f"GenoPC_{i + 1}" for i in range(ngeno)]
out = pd.concat([base.loc[rows], pcs])
out.index.name = "id"
out.to_csv(f"{BASE}/{sub}/{ct}/tensorqtl.{cov}_optCov.txt", sep="\t")
print(ct, mod, list(out.index))
