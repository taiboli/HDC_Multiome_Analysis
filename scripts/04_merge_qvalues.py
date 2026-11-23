"""
Combine per-chromosome permutation results for each cell type x modality, compute genome-wide
q-values (Storey qvalue via tensorqtl.post.calculate_qvalues; needs `module load R/4.2.0`), and
check that every per-chromosome nominal parquet is present and readable.

Writes:
  {eQTL,caQTL}_Output/<CellType>/allchr_cis_noBootstrap_TMM_invnorm.cis_qtl.qval.txt.gz
  HDC_Filtered/QTL_summary.tsv
"""
import glob
import os

import pandas as pd
import pyarrow.parquet as pq
from tensorqtl import post, read_phenotype_bed
from tensorqtl.core import SimpleLogger

BASE = os.environ["HDC_DIR"] + "/"
TAG = "noBootstrap_TMM_invnorm"
COVSET = os.environ.get("COVSET", "")
SUFFIX = f"_{COVSET}" if COVSET else ""

rows = []
for modality, in_sub, out_sub, bed in [("RNA", "eQTL", f"eQTL_Output{SUFFIX}", "geneExpr"),
                                        ("ATAC", "caQTL", f"caQTL_Output{SUFFIX}", "atacExpr")]:
    for pheno_file in sorted(glob.glob(f"{BASE}/{in_sub}/*/tensorqtl.{bed}_{TAG}.bed.gz")):
        celltype = pheno_file.split("/")[-2]
        out_dir = f"{BASE}/{out_sub}/{celltype}/"
        _, pos_df = read_phenotype_bed(pheno_file)
        chroms = [c for c in map(str, range(1, 23)) if (pos_df["chr"] == c).any()]

        missing, perm = [], []
        for c in chroms:
            f_perm = f"{out_dir}/chr{c}_cis_{TAG}.cis_qtl.txt.gz"
            f_nom = f"{out_dir}/chr{c}_cis_nominal_{TAG}.cis_qtl_pairs.{c}.parquet"
            if not os.path.exists(f_perm):
                missing.append(f"chr{c}:perm")
            else:
                perm.append(pd.read_csv(f_perm, sep="\t", index_col=0))
            try:
                pq.read_metadata(f_nom)
            except Exception:
                missing.append(f"chr{c}:nominal")

        row = {"CellType": celltype, "Modality": modality,
               "nPhenotypeInput": int(pos_df["chr"].isin(chroms).sum()), "Missing": ",".join(missing)}
        if perm and not missing:
            res_df = pd.concat(perm)
            logger = SimpleLogger(f"{out_dir}/allchr_cis_{TAG}.qvalue.log")
            post.calculate_qvalues(res_df, fdr=0.05, logger=logger)
            res_df.to_csv(f"{out_dir}/allchr_cis_{TAG}.cis_qtl.qval.txt.gz", sep="\t", float_format="%.6g")
            row.update(nTested=len(res_df), nSigFDR05=int((res_df["qval"] <= 0.05).sum()),
                       nSigFDR10=int((res_df["qval"] <= 0.10).sum()))
        rows.append(row)
        print(row, flush=True)

pd.DataFrame(rows).to_csv(f"{BASE}/QTL_summary{SUFFIX}.tsv", sep="\t", index=False)
