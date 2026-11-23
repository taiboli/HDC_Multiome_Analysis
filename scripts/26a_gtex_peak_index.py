"""
GTEx v8 support index for HDC peak-gene pairs: for every (gene, HDC peak) the number of GTEx v8 tissues (49;
$GTEX_DIR, read-only) in which a significant eQTL variant of that gene
(signif_variant_gene_pairs, b38) lies inside the peak, and the best (min) nominal p across tissues.

Output: PeakAnnot/ref/gtex_peak_gene.parquet (peak, gene_id, n_tissues, min_p)
"""
import os
import glob

import numpy as np
import pandas as pd

GTEX = os.environ["GTEX_DIR"] + "/"
PA = os.environ["HDC_DIR"] + "/PeakAnnot/"
pk = pd.read_csv(os.environ["HDC_DIR"] + "/PseudoBulk/features_ATAC.tsv.gz", sep="\t")
pk = pk[pk.chr.str.match(r"^chr[0-9X]+$")]
idx = {c: d.sort_values("start") for c, d in pk.groupby("chr")}
parts = []
for f in sorted(glob.glob(f"{GTEX}/*.v8.signif_variant_gene_pairs.txt.gz")):
    t = f.split("/")[-1].split(".")[0]
    d = pd.read_csv(f, sep="\t", usecols=["variant_id", "gene_id", "pval_nominal"])
    v = d.variant_id.str.split("_", expand=True)
    d = pd.DataFrame({"chr": v[0], "pos": v[1].astype(int) - 1, "gene_id": d.gene_id.str.split(".").str[0],
                      "p": d.pval_nominal})
    out = []
    for c, x in d.groupby("chr"):
        if c not in idx:
            continue
        p = idx[c]
        j = np.searchsorted(p.start.to_numpy(), x.pos.to_numpy(), side="right") - 1
        ok = (j >= 0) & (x.pos.to_numpy() < p.end.to_numpy()[np.maximum(j, 0)])
        out.append(pd.DataFrame({"peak": p.id.to_numpy()[j[ok]], "gene_id": x.gene_id.to_numpy()[ok], "p": x.p.to_numpy()[ok]}))
    out = pd.concat(out).groupby(["peak", "gene_id"]).p.min().reset_index().assign(tissue=t)
    parts.append(out)
    print(t, len(out), flush=True)
g = pd.concat(parts)
g = g.groupby(["peak", "gene_id"]).agg(n_tissues=("tissue", "nunique"), min_p=("p", "min")).reset_index()
g.to_parquet(f"{PA}/ref/gtex_peak_gene.parquet", index=False)
print(len(g), "peak-gene pairs with GTEx support;", g.peak.nunique(), "peaks")
