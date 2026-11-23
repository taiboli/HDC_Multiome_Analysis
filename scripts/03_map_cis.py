"""
tensorQTL cis mapping (permutation + nominal) for one HDC cell type x modality, per chromosome.

Mirrors the EB_CellTypeData submit scripts (python -m tensorqtl ... --mode cis / cis_nominal,
default 10,000 permutations, +/-1 Mb window, maf_threshold 0, eigenMT in nominal mode), but loads
each chromosome's genotypes once for both modes, and skips chromosomes whose outputs already exist.

Output per chromosome in {eQTL,caQTL}_Output/<CellType>/:
  chr{N}_cis_noBootstrap_TMM_invnorm.cis_qtl.txt.gz                      permutation (per-phenotype lead)
  chr{N}_cis_noBootstrap_TMM_invnorm.tensorQTL.cis.log
  chr{N}_cis_nominal_noBootstrap_TMM_invnorm.cis_qtl_pairs.{N}.parquet    all cis variant-phenotype pairs
  chr{N}_cis_nominal_noBootstrap_TMM_invnorm.tensorQTL.cis_nominal.log
q-values are computed genome-wide afterwards (04_merge_qvalues.py), not per chromosome.

Usage: python 03_map_cis.py <CellType> <RNA|ATAC> [chr list, default 1-22]
"""
import os
import sys
import time
from datetime import datetime

import pandas as pd
import torch
from tensorqtl import cis, genotypeio, read_phenotype_bed
from tensorqtl.core import SimpleLogger

GENO = os.environ["GENOTYPE_DIR"] + "/Genotype/plink/chr{}"
BASE = os.environ["HDC_DIR"] + "/"
TAG = "noBootstrap_TMM_invnorm"
SEED = 2025
# COVSET=optCov (env) -> optimized covariates (05_optimize_covariates.py), outputs in *_Output_optCov/
COVSET = os.environ.get("COVSET", "")
COV_TAG = COVSET if COVSET else TAG
OUT_SUFFIX = f"_{COVSET}" if COVSET else ""

celltype, modality = sys.argv[1], sys.argv[2]
chroms = sys.argv[3].split(",") if len(sys.argv) > 3 else [str(i) for i in range(1, 23)]

if modality == "RNA":
    in_dir = f"{BASE}/eQTL/{celltype}/"
    pheno_file = f"{in_dir}/tensorqtl.geneExpr_{TAG}.bed.gz"
    cov_file = f"{in_dir}/tensorqtl.geneCov_{COV_TAG}.txt"
    out_dir = f"{BASE}/eQTL_Output{OUT_SUFFIX}/{celltype}/"
else:
    in_dir = f"{BASE}/caQTL/{celltype}/"
    pheno_file = f"{in_dir}/tensorqtl.atacExpr_{TAG}.bed.gz"
    cov_file = f"{in_dir}/tensorqtl.atacCov_{COV_TAG}.txt"
    out_dir = f"{BASE}/caQTL_Output{OUT_SUFFIX}/{celltype}/"

if not os.path.exists(pheno_file):
    print(f"No input for {celltype} {modality} (skipped at input step); nothing to do.")
    sys.exit(0)
os.makedirs(out_dir, exist_ok=True)

device = "GPU " + torch.cuda.get_device_name(0) if torch.cuda.is_available() else "CPU"
print(f"{celltype} {modality} on {device}", flush=True)

phenotype_df, phenotype_pos_df = read_phenotype_bed(pheno_file)
covariates_df = pd.read_csv(cov_file, sep="\t", index_col=0).T
assert phenotype_df.columns.equals(covariates_df.index)

for chrom in chroms:
    perm_prefix = f"chr{chrom}_cis_{TAG}"
    nom_prefix = f"chr{chrom}_cis_nominal_{TAG}"
    perm_out = f"{out_dir}/{perm_prefix}.cis_qtl.txt.gz"
    nom_out = f"{out_dir}/{nom_prefix}.cis_qtl_pairs.{chrom}.parquet"
    if os.path.exists(perm_out) and os.path.exists(nom_out):
        print(f"chr{chrom}: done already, skipping", flush=True)
        continue

    sel = phenotype_pos_df["chr"] == chrom
    if sel.sum() == 0:
        print(f"chr{chrom}: no phenotypes", flush=True)
        continue
    p_df, p_pos_df = phenotype_df.loc[sel], phenotype_pos_df.loc[sel]

    t0 = time.time()
    genotype_df, variant_df = genotypeio.load_genotypes(GENO.format(chrom), select_samples=phenotype_df.columns)
    # drop variants monomorphic in these donors: they only give NaN correlations, and a cis-window
    # containing only such variants makes map_cis fail ("All correlations resulted in NaN")
    poly = genotype_df.nunique(axis=1) > 1
    genotype_df, variant_df = genotype_df[poly], variant_df[poly.values]

    if not os.path.exists(perm_out):
        logger = SimpleLogger(f"{out_dir}/{perm_prefix}.tensorQTL.cis.log")
        logger.write(f'[{datetime.now().strftime("%b %d %H:%M:%S")}] tensorQTL cis-QTL mapping on {device}')
        logger.write(f"  * phenotypes: {pheno_file}\n  * covariates: {cov_file}\n  * seed {SEED}")
        res_df = cis.map_cis(genotype_df, variant_df, p_df, p_pos_df, covariates_df=covariates_df,
                             nperm=10000, window=1000000, maf_threshold=0, logger=logger,
                             seed=SEED, verbose=True)
        res_df.to_csv(perm_out + ".tmp", sep="\t", float_format="%.6g", compression="gzip")
        os.replace(perm_out + ".tmp", perm_out)
        logger.write(f'[{datetime.now().strftime("%b %d %H:%M:%S")}] Finished mapping')

    if not os.path.exists(nom_out):
        logger = SimpleLogger(f"{out_dir}/{nom_prefix}.tensorQTL.cis_nominal.log")
        logger.write(f'[{datetime.now().strftime("%b %d %H:%M:%S")}] tensorQTL cis_nominal mapping on {device}')
        cis.map_nominal(genotype_df, variant_df, p_df, p_pos_df, nom_prefix, covariates_df=covariates_df,
                        window=1000000, maf_threshold=0, run_eigenmt=True, output_dir=out_dir,
                        write_top=True, write_stats=True, logger=logger, verbose=True)
        logger.write(f'[{datetime.now().strftime("%b %d %H:%M:%S")}] Finished mapping')

    print(f"chr{chrom}: {sel.sum()} phenotypes, {time.time() - t0:.0f}s", flush=True)
    del genotype_df, variant_df
