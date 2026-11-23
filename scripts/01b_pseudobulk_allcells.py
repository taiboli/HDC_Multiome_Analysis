"""
All-cell ("AllCells") pseudobulk per donor, for a pooled QTL map across the heterogeneous culture.

Cell types partition the labeled paired cells, so summing the per-cell-type outputs of 01_pseudobulk.py
(all 11 highres_celltype levels, incl. Ectoderm_2) is exactly the per-donor pseudobulk over all cells:
  AllCells_{RNA|ATAC}_SumCount.tsv.gz, AllCells_{RNA|ATAC}_NCellExpr.tsv.gz, AllCells_NCell.tsv  in PseudoBulk/
Also writes AllCells_composition.tsv (donor x cell-type cell fractions) for interpreting pooled PCs.
"""
import glob
import os

import pandas as pd

PB = os.environ["HDC_DIR"] + "/PseudoBulk/"
cts = sorted(os.path.basename(f)[:-len("_NCell.tsv")] for f in glob.glob(f"{PB}/*_NCell.tsv")
             if not os.path.basename(f).startswith("AllCells"))
print(cts, flush=True)

nc = {ct: pd.read_csv(f"{PB}/{ct}_NCell.tsv", sep="\t", index_col="donor").ncell for ct in cts}
comp = pd.DataFrame(nc).fillna(0).astype(int)
tot = comp.sum(axis=1)
tot.rename("ncell").rename_axis("donor").reset_index().to_csv(f"{PB}/AllCells_NCell.tsv", sep="\t", index=False)
comp.div(tot, axis=0).round(4).rename_axis("donor").to_csv(f"{PB}/AllCells_composition.tsv", sep="\t")
donors = list(tot.index)

for mod in ("RNA", "ATAC"):
    for kind in ("SumCount", "NCellExpr"):
        acc = None
        for ct in cts:
            d = pd.read_csv(f"{PB}/{ct}_{mod}_{kind}.tsv.gz", sep="\t", index_col="id")
            d = d.reindex(columns=donors, fill_value=0)
            acc = d if acc is None else acc.add(d, fill_value=0)
        feats = pd.read_csv(f"{PB}/features_{mod}.tsv.gz", sep="\t").id
        assert acc.index.equals(pd.Index(feats)), "feature order differs from features file"
        acc.astype("int64").reset_index().to_csv(f"{PB}/AllCells_{mod}_{kind}.tsv.gz", sep="\t", index=False)
        print(mod, kind, acc.shape, flush=True)
print(tot.describe())
