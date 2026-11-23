"""
HOMER inputs for known-motif enrichment of cell-type-specific HDC peaks (25b_homer.sh).
Target: top 3000 peaks of each cell type's peakspec set ranked by specificity (CPM_ct + 0.1) / (mean CPM other + 0.1);
background: 30000 random HDC peaks not in that cell type's peakspec set (seed 1). 200 bp centered on the peak center
(HOMER -size given). Linked-peak sets (26_link_analysis.py) are written to the same folder when available.
Output: PeakAnnot/homer/input/<set>.target.bed, <set>.bg.bed
"""
import os

import numpy as np
import pandas as pd

PA = os.environ["HDC_DIR"] + "/PeakAnnot/"
D = f"{PA}/homer/input/"
os.makedirs(D, exist_ok=True)
pk = pd.read_parquet(f"{PA}/peak_annotation.parquet")
cts = [c[5:] for c in pk.columns if c.startswith("spec_")]
cpm = pk[[f"cpm_{c}" for c in cts]].to_numpy()
rng = np.random.default_rng(1)
for i, ct in enumerate(cts):
    if ct == "Unknown":
        continue
    other = np.delete(cpm, i, axis=1).mean(1)
    s = (cpm[:, i] + 0.1) / (other + 0.1)
    m = pk[f"spec_{ct}"].to_numpy()
    tgt = pk[m].assign(s=s[m]).sort_values("s", ascending=False).head(3000)
    bg = pk[~m].sample(30000, random_state=1)
    for name, d in (("target", tgt), ("bg", bg)):
        d.assign(strand="+")[["chr", "start", "end", "peak", "width", "strand"]].to_csv(
            f"{D}/peakspec_{ct}.{name}.bed", sep="\t", header=False, index=False)
print(sorted(os.listdir(D)))
