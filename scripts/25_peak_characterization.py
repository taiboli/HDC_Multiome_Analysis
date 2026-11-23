"""
Characterize the HDC ATAC peak set (figure panels + tables) from 24_peak_annotation.py (HDC peaks and the
width/chromosome-matched shuffled null, --shuffle).

  overview       genomic category, ENCODE cCRE v4 class, phastCons overlap: HDC vs shuffled
  chromhmm       per Roadmap epigenome x 18 state: log2( frac HDC peak centers in state / frac shuffled ), and the
                 fraction of peaks in active states (1-11) per epigenome
  specificity    tau deciles vs promoter / cCRE / ESC-active / conservation
  ct_epigenome   cell-type-specific peaks (peakspec, 11_celltype_peak_enrichment.py): log2 enrichment of enhancer
                 states (7-11) in each epigenome relative to all HDC peaks; Fisher test per CT x epigenome
  poised         CT-specific peaks in bivalent / Polycomb states (14-17) of H1 (E003) and H9 (E008)
  novel          CT-specific peaks absent from the ENCODE registry: conservation and matched-epigenome enhancer state

Output: PeakAnnot/results/*.tsv, PeakAnnot/figures/fig_peaks_*.pdf
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402
from scipy.stats import fisher_exact  # noqa: E402

BASE = os.environ["HDC_DIR"] + "/"
PA = f"{BASE}/PeakAnnot/"
RES, FIG = f"{PA}/results/", f"{PA}/figures/"
os.makedirs(RES, exist_ok=True)
os.makedirs(FIG, exist_ok=True)
STATES = ["TssA", "TssFlnk", "TssFlnkU", "TssFlnkD", "Tx", "TxWk", "EnhG1", "EnhG2", "EnhA1", "EnhA2", "EnhWk",
          "ZNF/Rpts", "Het", "TssBiv", "EnhBiv", "ReprPC", "ReprPCWk", "Quies"]
ENH = [7, 8, 9, 10, 11]
ACTIVE = list(range(1, 12))
POISED = [14, 15, 16, 17]
CTS = ["iPSC", "Ectoderm_1", "Ectoderm_2", "Ectoderm_3", "Neural_crest", "Neuronal", "Retinal", "Mesoderm",
       "Vascular_endothelial", "Endoderm"]
CATS = ["promoter", "exon", "intron", "distal"]
CCRE = ["PLS", "pELS", "dELS", "CA-H3K4me3", "CA-CTCF", "CA-TF", "CA", "TF", "none"]
plt.rcParams.update({"font.size": 8, "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False})

pk = pd.read_parquet(f"{PA}/peak_annotation.parquet")
sh = pd.read_parquet(f"{PA}/peak_annotation_shuffled.parquet")
st = pd.read_parquet(f"{PA}/peak_chromhmm18.parquet").set_index("peak").reindex(pk.peak)
sst = pd.read_parquet(f"{PA}/peak_chromhmm18_shuffled.parquet").set_index("peak").reindex(sh.peak)
meta = pd.read_csv(f"{PA}/ref/roadmap_metadata.tsv", sep="\t").set_index("EID")
eids = [e for e in st.columns if e in meta.index]
grp_order = ["ESC", "iPSC", "ES-deriv", "Neurosph", "Brain", "Heart", "Muscle", "Smooth Muscle", "Digestive", "Other",
             "Epithelial", "Thymus", "HSC & B-cell", "Blood & T-cell", "ENCODE", "Sm. Muscle", "IMR90"]
eids = sorted(eids, key=lambda e: (grp_order.index(meta.GROUP[e]) if meta.GROUP[e] in grp_order else 99, e))
lab = {e: f"{e} {meta.MNEMONIC[e]}" for e in eids}
pk["phast"] = pk.phastcons_bp > 0
sh["phast"] = sh.phastcons_bp > 0
pk["in_ccre"] = pk.ccre != "none"
sh["in_ccre"] = sh.ccre != "none"

# ---------------- 1. overview ----------------
ov = []
for name, d in (("HDC peaks", pk), ("shuffled", sh)):
    r = {"set": name, "n": len(d), "median_width": d.width.median(), "gc": d.gc.mean(),
         "frac_ccre": d.in_ccre.mean(), "frac_phastcons": d.phast.mean(),
         "median_abs_dist_tss": d.dist_tss.abs().median()}
    r.update({f"frac_{c}": (d.category == c).mean() for c in CATS})
    r.update({f"ccre_{c}": (d.ccre == c).mean() for c in CCRE})
    ov.append(r)
ov = pd.DataFrame(ov)
ov.to_csv(f"{RES}/peaks_overview.tsv", sep="\t", index=False, float_format="%.4g")
print(ov.T.to_string(), flush=True)

fig, ax = plt.subplots(1, 4, figsize=(11, 2.6))
ax[0].hist(pk.width.clip(upper=2000), bins=60, color="#4C72B0")
ax[0].set(xlabel="peak width (bp)", ylabel="peaks", title=f"{len(pk):,} HDC peaks")
x = np.arange(len(CATS))
for i, (name, d, c) in enumerate((("HDC", pk, "#4C72B0"), ("shuffled", sh, "#BBBBBB"))):
    ax[1].bar(x + (i - 0.5) * 0.4, [(d.category == k).mean() for k in CATS], 0.4, color=c, label=name)
ax[1].set_xticks(x, CATS, rotation=30)
ax[1].set(ylabel="fraction of peaks", title="genomic context")
ax[1].legend(frameon=False)
x = np.arange(len(CCRE))
for i, (name, d, c) in enumerate((("HDC", pk, "#4C72B0"), ("shuffled", sh, "#BBBBBB"))):
    ax[2].bar(x + (i - 0.5) * 0.4, [(d.ccre == k).mean() for k in CCRE], 0.4, color=c, label=name)
ax[2].set_xticks(x, CCRE, rotation=60)
ax[2].set(ylabel="fraction of peaks", title=f"ENCODE cCRE v4 ({pk.in_ccre.mean():.0%} vs {sh.in_ccre.mean():.0%})")
dd = np.sign(pk.dist_tss) * np.log10(pk.dist_tss.abs() + 1)
ax[3].hist(dd, bins=100, color="#4C72B0")
ax[3].set(xlabel="signed log10 distance to nearest TSS", ylabel="peaks", title="distance to TSS")
fig.tight_layout()
fig.savefig(f"{FIG}/fig_peaks_overview.pdf")

# ---------------- 2. ChromHMM enrichment over shuffled, all HDC peaks ----------------
rows = []
for e in eids:
    a = np.bincount(st[e].to_numpy(), minlength=19)[1:] / len(st)
    b = np.bincount(sst[e].to_numpy(), minlength=19)[1:] / len(sst)
    for k in range(18):
        rows.append({"eid": e, "group": meta.GROUP[e], "name": meta.STD_NAME[e], "state": STATES[k], "frac_hdc": a[k],
                     "frac_shuffled": b[k], "log2_enr": np.log2((a[k] + 1e-5) / (b[k] + 1e-5))})
ch = pd.DataFrame(rows)
ch.to_csv(f"{RES}/chromhmm18_enrichment_allpeaks.tsv", sep="\t", index=False, float_format="%.4g")
act = ch[ch.state.isin([STATES[i - 1] for i in ACTIVE])].groupby(["eid", "group"]).frac_hdc.sum().reset_index()
act = act.sort_values("frac_hdc", ascending=False)
act.to_csv(f"{RES}/chromhmm18_active_frac_by_epigenome.tsv", sep="\t", index=False, float_format="%.4g")
print("epigenomes with most HDC peaks in active states:\n", act.head(10).to_string(), flush=True)

M = ch.pivot(index="state", columns="eid", values="log2_enr").loc[STATES, eids]
fig = plt.figure(figsize=(14, 5.2))
gs = fig.add_gridspec(2, 1, height_ratios=[1, 3], hspace=0.05)
a0 = fig.add_subplot(gs[0])
pal = dict(zip(grp_order, sns.color_palette("tab20", len(grp_order))))
af = act.set_index("eid").frac_hdc.reindex(eids)
a0.bar(np.arange(len(eids)), af, color=[pal.get(meta.GROUP[e], "grey") for e in eids])
a0.set(xlim=(-0.5, len(eids) - 0.5), xticks=[], ylabel="frac peaks\nactive (1-11)")
a1 = fig.add_subplot(gs[1])
sns.heatmap(M, cmap="RdBu_r", center=0, vmin=-3, vmax=3, ax=a1, cbar_kws={"label": "log2(HDC / shuffled)", "shrink": 0.6})
a1.set_xticks(np.arange(len(eids)) + 0.5, [lab[e] for e in eids], rotation=90, fontsize=5)
a1.set(xlabel="", ylabel="")
handles = [plt.Rectangle((0, 0), 1, 1, color=pal[g]) for g in grp_order if g in set(meta.GROUP[eids])]
a0.legend(handles, [g for g in grp_order if g in set(meta.GROUP[eids])], ncol=9, fontsize=6, frameon=False,
          loc="lower center", bbox_to_anchor=(0.5, 1.0))
fig.savefig(f"{FIG}/fig_peaks_chromhmm_allpeaks.pdf", bbox_inches="tight")

# ---------------- 3. specificity deciles ----------------
pk["tau_bin"] = pd.qcut(pk.tau.rank(method="first"), 10, labels=[f"D{i}" for i in range(1, 11)])
esc = ["E003", "E008"]
pk["esc_active"] = st[esc].isin(ACTIVE).any(axis=1).to_numpy()
pk["esc_poised"] = st[esc].isin(POISED).any(axis=1).to_numpy() & ~pk.esc_active
sp = pk.groupby("tau_bin", observed=True).agg(tau=("tau", "median"), promoter=("category", lambda x: (x == "promoter").mean()),
                                              ccre=("in_ccre", "mean"), esc_active=("esc_active", "mean"),
                                              esc_poised=("esc_poised", "mean"), phastcons=("phast", "mean"),
                                              width=("width", "median")).reset_index()
sp.to_csv(f"{RES}/specificity_deciles.tsv", sep="\t", index=False, float_format="%.4g")
print(sp.to_string(), flush=True)

# ---------------- 4. CT-specific peaks x epigenome enhancer enrichment ----------------
in_enh = st.isin(ENH).to_numpy()  # peaks x epigenomes
cols = list(st.columns)
base_rate = in_enh.mean(0)
rows = []
for ct in CTS:
    m = pk[f"spec_{ct}"].to_numpy()
    a = in_enh[m].sum(0)
    for j, e in enumerate(cols):
        if e not in meta.index:
            continue
        tab = [[a[j], m.sum() - a[j]], [in_enh[~m, j].sum(), (~m).sum() - in_enh[~m, j].sum()]]
        orr, p = fisher_exact(tab)
        rows.append({"celltype": ct, "eid": e, "group": meta.GROUP[e], "name": meta.STD_NAME[e], "n_spec": int(m.sum()),
                     "frac_enh_spec": a[j] / m.sum(), "frac_enh_all": base_rate[j],
                     "log2_enr": np.log2(a[j] / m.sum() / base_rate[j]), "OR": orr, "p": p})
ce = pd.DataFrame(rows)
ce.to_csv(f"{RES}/ctspec_peaks_epigenome_enhancer_enrichment.tsv", sep="\t", index=False, float_format="%.4g")
top = ce.sort_values("log2_enr", ascending=False).groupby("celltype").head(5)
top.sort_values(["celltype", "log2_enr"], ascending=[True, False]).to_csv(f"{RES}/ctspec_peaks_top_epigenomes.tsv",
                                                                         sep="\t", index=False, float_format="%.4g")
print(top.sort_values(["celltype", "log2_enr"], ascending=[True, False])[["celltype", "eid", "name", "log2_enr", "p"]].to_string(), flush=True)

H = ce.pivot(index="celltype", columns="eid", values="log2_enr").loc[CTS, eids]
fig, a = plt.subplots(figsize=(14, 3.4))
sns.heatmap(H, cmap="RdBu_r", center=0, vmin=-1.5, vmax=1.5, ax=a, cbar_kws={"label": "log2 enhancer-state enrichment\n(CT-specific vs all HDC peaks)", "shrink": 0.7})
a.set_xticks(np.arange(len(eids)) + 0.5, [lab[e] for e in eids], rotation=90, fontsize=5)
for j, e in enumerate(eids):
    a.add_patch(plt.Rectangle((j, len(CTS)), 1, 0.35, color=pal.get(meta.GROUP[e], "grey"), clip_on=False))
a.set(xlabel="", ylabel="")
fig.savefig(f"{FIG}/fig_peaks_ctspec_epigenome.pdf", bbox_inches="tight")

# compact: selected developmental / matched epigenomes
SEL = ["E003", "E008", "E020", "E021", "E004", "E005", "E006", "E007", "E009", "E010", "E011", "E012", "E013", "E053",
       "E054", "E081", "E082", "E083", "E058", "E059", "E061", "E065", "E066", "E106", "E109", "E122", "E055", "E057", "E129"]
SEL = [e for e in SEL if e in H.columns]
fig, a = plt.subplots(figsize=(8, 3.4))
sns.heatmap(H[SEL], cmap="RdBu_r", center=0, vmin=-1.5, vmax=1.5, ax=a, cbar_kws={"label": "log2 enrichment", "shrink": 0.7})
a.set_xticks(np.arange(len(SEL)) + 0.5, [f"{meta.STD_NAME[e][:28]}" for e in SEL], rotation=90, fontsize=6)
a.set(xlabel="", ylabel="", title="Enhancer states (EnhG/EnhA/EnhWk) at cell-type-specific HDC peaks")
fig.savefig(f"{FIG}/fig_peaks_ctspec_epigenome_selected.pdf", bbox_inches="tight")

# ---------------- 5. poised in ESC ----------------
rows = []
for ct in CTS + ["all"]:
    m = np.ones(len(pk), bool) if ct == "all" else pk[f"spec_{ct}"].to_numpy()
    for e in esc:
        s = st[e].to_numpy()[m]
        rows.append({"celltype": ct, "eid": e, "active": np.isin(s, ACTIVE).mean(), "bivalent": np.isin(s, [14, 15]).mean(),
                     "polycomb": np.isin(s, [16, 17]).mean(), "quiescent_het": np.isin(s, [12, 13, 18]).mean()})
po = pd.DataFrame(rows)
po.to_csv(f"{RES}/ctspec_peaks_esc_states.tsv", sep="\t", index=False, float_format="%.4g")
print(po[po.eid == "E003"].to_string(), flush=True)
fig, a = plt.subplots(figsize=(5, 2.6))
p3 = po[po.eid == "E003"].set_index("celltype")[["active", "bivalent", "polycomb", "quiescent_het"]]
p3.plot.barh(stacked=True, ax=a, color=["#E41A1C", "#984EA3", "#377EB8", "#CCCCCC"], width=0.8)
a.set(xlabel="fraction of peaks (H1 ESC, E003 ChromHMM)", ylabel="")
a.legend(frameon=False, fontsize=6, ncol=4, loc="lower center", bbox_to_anchor=(0.5, 1))
fig.tight_layout()
fig.savefig(f"{FIG}/fig_peaks_esc_states.pdf")

# ---------------- 6. novel (not in ENCODE registry) CT-specific peaks ----------------
best_e = ce.sort_values("log2_enr", ascending=False).groupby("celltype").eid.first()
rows = []
for ct in CTS:
    m = pk[f"spec_{ct}"].to_numpy()
    nov = m & ~pk.in_ccre.to_numpy()
    e = best_e[ct]
    rows.append({"celltype": ct, "n_spec": int(m.sum()), "frac_novel": nov.sum() / m.sum(),
                 "phast_novel": pk.phast[nov].mean(), "phast_known": pk.phast[m & pk.in_ccre].mean(),
                 "phast_shuffled": sh.phast.mean(), "best_epigenome": e,
                 "enh_best_epigenome_novel": np.isin(st[e].to_numpy()[nov], ENH).mean(),
                 "enh_best_epigenome_shuffled": np.isin(sst[e].to_numpy(), ENH).mean()})
nv = pd.DataFrame(rows)
nv.to_csv(f"{RES}/ctspec_novel_peaks.tsv", sep="\t", index=False, float_format="%.4g")
print(nv.to_string(), flush=True)
pk[["peak", "chr", "start", "end", "width", "gc", "category", "nearest_gene", "dist_tss", "ccre", "phastcons_bp", "tau",
    "top_ct"] + [c for c in pk.columns if c.startswith("spec_") or c.startswith("cpm_")]].to_csv(
    f"{RES}/peak_annotation_table.tsv.gz", sep="\t", index=False, float_format="%.4g")
print("done")
