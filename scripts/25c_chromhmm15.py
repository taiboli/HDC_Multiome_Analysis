"""
Roadmap 15-state core-marks ChromHMM (127 epigenomes, hg38 liftover) at HDC peak centers and at the shuffled null.
Complements the 18-state core_K27ac model (98 epigenomes) of 24/25, which lacks the neuronal epigenomes
(H9-derived neurons / NPC E009-E010, fetal brain E081-E082, neurospheres E053-E054 have no H3K27ac).

  all-peak enrichment   log2(frac HDC centers in state / frac shuffled centers) per epigenome x state
  CT-specific peaks     enhancer states (6_EnhG, 7_Enh) relative to all HDC peaks, Fisher test, per CT x epigenome

Output: PeakAnnot/peak_chromhmm15{,_shuffled}.parquet, results/chromhmm15_enrichment_allpeaks.tsv,
        results/ctspec_peaks_epigenome15_enhancer_enrichment.tsv, results/ctspec_peaks_top_epigenomes15.tsv,
        figures/fig_peaks_ctspec_epigenome15.pdf
"""
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402
from scipy.stats import fisher_exact  # noqa: E402

PA = os.environ["HDC_DIR"] + "/PeakAnnot/"
RES, FIG = f"{PA}/results/", f"{PA}/figures/"
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
STATES = ["TssA", "TssAFlnk", "TxFlnk", "Tx", "TxWk", "EnhG", "Enh", "ZNF/Rpts", "Het", "TssBiv", "BivFlnk", "EnhBiv",
          "ReprPC", "ReprPCWk", "Quies"]
ENH = [6, 7]
CTS = ["iPSC", "Ectoderm_1", "Ectoderm_2", "Ectoderm_3", "Neural_crest", "Neuronal", "Retinal", "Mesoderm",
       "Vascular_endothelial", "Endoderm"]
plt.rcParams.update({"font.size": 8, "pdf.fonttype": 42})


def annotate(pk, files):
    out = {}
    for f in files:
        eid = os.path.basename(f).split("_")[0]
        s = pd.read_csv(f, sep="\t", header=None, skiprows=1, usecols=[0, 1, 2, 3], names=["chr", "start", "end", "state"])
        s = s[s.chr.isin(CHROMS)].sort_values(["chr", "start"])
        s["code"] = s.state.str.split("_").str[0].astype(np.int8)
        v = np.zeros(len(pk), np.int8)
        for c, d in pk.groupby("chr"):
            sc = s[s.chr == c]
            j = np.searchsorted(sc.start.to_numpy(), d.center.to_numpy(), side="right") - 1
            ok = (j >= 0) & (d.center.to_numpy() < sc.end.to_numpy()[np.maximum(j, 0)])
            v[d.index[ok]] = sc.code.to_numpy()[j[ok]]
        out[eid] = v
    st = pd.DataFrame(out)
    st.insert(0, "peak", pk.peak.to_numpy())
    return st


files = sorted(glob.glob(f"{PA}/ref/chromhmm15/E*_15_coreMarks_hg38lift_dense.bed.gz"))
pk = pd.read_parquet(f"{PA}/peak_annotation.parquet")
sh = pd.read_parquet(f"{PA}/peak_annotation_shuffled.parquet")
st = annotate(pk, files)
sst = annotate(sh, files)
st.to_parquet(f"{PA}/peak_chromhmm15.parquet", index=False)
sst.to_parquet(f"{PA}/peak_chromhmm15_shuffled.parquet", index=False)
st, sst = st.set_index("peak"), sst.set_index("peak")
meta = pd.read_csv(f"{PA}/ref/roadmap_metadata.tsv", sep="\t").set_index("EID")
eids = [e for e in st.columns if e in meta.index]

rows = []
for e in eids:
    a = np.bincount(st[e].to_numpy(), minlength=16)[1:] / len(st)
    b = np.bincount(sst[e].to_numpy(), minlength=16)[1:] / len(sst)
    for k in range(15):
        rows.append({"eid": e, "group": meta.GROUP[e], "name": meta.STD_NAME[e], "state": STATES[k], "frac_hdc": a[k],
                     "frac_shuffled": b[k], "log2_enr": np.log2((a[k] + 1e-5) / (b[k] + 1e-5))})
pd.DataFrame(rows).to_csv(f"{RES}/chromhmm15_enrichment_allpeaks.tsv", sep="\t", index=False, float_format="%.4g")

in_enh = st[eids].isin(ENH).to_numpy()
base = in_enh.mean(0)
rows = []
for ct in CTS:
    m = pk[f"spec_{ct}"].to_numpy()
    a = in_enh[m].sum(0)
    for j, e in enumerate(eids):
        tab = [[a[j], m.sum() - a[j]], [in_enh[~m, j].sum(), (~m).sum() - in_enh[~m, j].sum()]]
        orr, p = fisher_exact(tab)
        rows.append({"celltype": ct, "eid": e, "group": meta.GROUP[e], "name": meta.STD_NAME[e], "frac_enh_spec": a[j] / m.sum(),
                     "frac_enh_all": base[j], "log2_enr": np.log2(a[j] / m.sum() / base[j]), "OR": orr, "p": p})
ce = pd.DataFrame(rows)
ce.to_csv(f"{RES}/ctspec_peaks_epigenome15_enhancer_enrichment.tsv", sep="\t", index=False, float_format="%.4g")
top = ce.sort_values("log2_enr", ascending=False).groupby("celltype").head(5).sort_values(["celltype", "log2_enr"], ascending=[True, False])
top.to_csv(f"{RES}/ctspec_peaks_top_epigenomes15.tsv", sep="\t", index=False, float_format="%.4g")
print(top[["celltype", "eid", "name", "log2_enr", "frac_enh_spec", "frac_enh_all"]].to_string(), flush=True)

grp_order = ["ESC", "iPSC", "ES-deriv", "Neurosph", "Brain", "Heart", "Muscle", "Smooth Muscle", "Digestive", "Other",
             "Epithelial", "Thymus", "HSC & B-cell", "Blood & T-cell", "ENCODE", "Sm. Muscle", "IMR90"]
eo = sorted(eids, key=lambda e: (grp_order.index(meta.GROUP[e]) if meta.GROUP[e] in grp_order else 99, e))
pal = dict(zip(grp_order, sns.color_palette("tab20", len(grp_order))))
H = ce.pivot(index="celltype", columns="eid", values="log2_enr").loc[CTS, eo]
fig, a = plt.subplots(figsize=(16, 3.4))
sns.heatmap(H, cmap="RdBu_r", center=0, vmin=-1.5, vmax=1.5, ax=a,
            cbar_kws={"label": "log2 enhancer-state enrichment\n(CT-specific vs all HDC peaks)", "shrink": 0.7})
a.set_xticks(np.arange(len(eo)) + 0.5, [f"{e} {meta.MNEMONIC[e]}" for e in eo], rotation=90, fontsize=4.5)
for j, e in enumerate(eo):
    a.add_patch(plt.Rectangle((j, len(CTS)), 1, 0.35, color=pal.get(meta.GROUP[e], "grey"), clip_on=False))
a.set(xlabel="", ylabel="", title="Roadmap 15-state (127 epigenomes): EnhG + Enh at cell-type-specific HDC peaks")
fig.savefig(f"{FIG}/fig_peaks_ctspec_epigenome15.pdf", bbox_inches="tight")

# per-CT: rank of lineage-matched epigenome groups
grp = ce.groupby(["celltype", "group"]).log2_enr.max().unstack()
grp.loc[CTS].to_csv(f"{RES}/ctspec_peaks_epigenome15_group_max.tsv", sep="\t", float_format="%.3g")
print(grp.loc[CTS].round(2).to_string())
