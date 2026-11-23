"""
Summary of S-LDSC for the SPIRE TF sets (BATCH=tf; 17_sldsc.sh -> 18_sldsc_summarize.py tf -> GWAS/ldsc/sldsc_tf.tsv).
Per set (primary traits only): peaks, traits, # with tau p < 0.05 (vs 5% expected; binomial), mean tau z (and one-sample
t-test of tau z > 0 across traits), median enrichment, best trait and its p / FDR.
tfprog_<CT> is conditioned on baselineLD + all HDC peaks + peakspec_<CT>; tfdyn_<lin> on baselineLD + all HDC peaks.
Output TFNet/results/sldsc_tf_summary.tsv
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import binomtest, ttest_1samp

BASE = os.environ["HDC_DIR"] + "/"
d = pd.read_csv(f"{BASE}/GWAS/ldsc/sldsc_tf.tsv", sep="\t")
d = d[d.primary]
sets = pd.read_csv(f"{BASE}/GWAS/ldsc/sets_tf.tsv.gz", sep="\t").annot.value_counts()
rows = []
for a, g in d.groupby("annot"):
    b = g.nsmallest(1, "p_tau").iloc[0]
    k = int((g.p_tau < 0.05).sum())
    rows.append(dict(set=a, conditioned_on="HDC + CT-specific peaks" if a.startswith("tfprog") else "HDC peaks",
                     peaks=int(sets.get(a, 0)), traits=len(g), n_p05=k,
                     p_excess=binomtest(k, len(g), 0.05, alternative="greater").pvalue,
                     mean_tau_z=g.tau_z.mean(), p_tau_z_gt0=ttest_1samp(g.tau_z, 0, alternative="greater").pvalue,
                     median_enrichment=g.enrichment.median(), best_trait=b.label, best_p=b.p_tau, best_FDR=b.FDR,
                     min_FDR=g.FDR.min()))
r = pd.DataFrame(rows)
CTS = ["iPSC", "Ectoderm_1", "Ectoderm_2", "Ectoderm_3", "Neural_crest", "Neuronal", "Retinal", "Mesoderm", "Endoderm",
       "Vascular_endothelial", "Unknown"]
want = [f"tfprog_{c}" for c in CTS] + [f"tfdyn_{x}" for x in ("neur", "neur_crest", "ret", "meso", "endo")]
r = r.set_index("set").reindex([w for w in want if w in set(r.set)]).reset_index()
r.to_csv(f"{BASE}/TFNet/results/sldsc_tf_summary.tsv", sep="\t", index=False, float_format="%.3g")
pd.set_option("display.width", 250)
print(r.to_string(index=False))
print("all sets: tests", int(r.traits.sum()), "p<0.05", int(r.n_p05.sum()), f"({r.n_p05.sum() / r.traits.sum():.3f})",
      "min FDR", round(r.min_FDR.min(), 3))
