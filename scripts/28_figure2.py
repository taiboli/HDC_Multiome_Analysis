"""
Figure 2: developmental QTLs and regulatory elements in adult traits.

  a  genome-wide pi1 (1 - pi0 of permutation pval_beta) per caQTL / eQTL map
  b  mashR QTLs per cell type (lfsr < 0.05)
  c  replication in adult brain neuron / glia caQTLs (pi1, sign concordance) for mashR vs null features
  d  S-LDSC of cell-type-specific peaks: -log10 p_tau for traits with any FDR < 10% hit (dot = FDR < 10%)
  e  colocalized loci (PP.H4 > 0.5): per-cell-type vs pooled maps (08c_coloc_audit.py, phenotype-locus level)
  f  example: Neural crest caQTL x angina / coronary atherosclerosis at 11q13 (mirror plot)
  g  GTEx v8 eGene status of HDC-colocalized genes vs background (12_gtex_background.py)
Output: Figures/Main/Figure2.{pdf,png}, panels in PeakAnnot/figures/panels/fig2_*.pdf
"""
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402
import seaborn as sns  # noqa: E402

BASE = os.environ["HDC_DIR"] + "/"
F = f"{BASE}/PeakAnnot/figures/"
MAINF = f"{BASE}/Figures/Main/"
os.makedirs(MAINF, exist_ok=True)
TAG = "noBootstrap_TMM_invnorm"
CTS = ["iPSC", "Ectoderm_1", "Ectoderm_3", "Neural_crest", "Neuronal", "Retinal", "Mesoderm", "Vascular_endothelial",
       "Endoderm", "Unknown"]
ALLCT = ["iPSC", "Ectoderm_1", "Ectoderm_2", "Ectoderm_3", "Neural_crest", "Neuronal", "Retinal", "Mesoderm",
         "Vascular_endothelial", "Endoderm"]
COL = dict(zip(ALLCT, sns.color_palette("tab10", len(ALLCT))))
COL.update({"Unknown": "#999999", "AllCells": "#000000"})
MODCOL = {"ATAC": "#E6550D", "RNA": "#3B6FB6"}
SHORT = {"Vascular_endothelial": "Vasc_endo"}
EX = dict(celltype="Neural_crest", modality="ATAC", phenotype_id="chr11_65586716_65587229", chr=11,
          traits=["FinnGen__I9_ANGINA", "FinnGen__I9_CORATHER"], lo=65_450_000, hi=65_750_000)
plt.rcParams.update({"font.size": 7, "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.titlesize": 7, "axes.titleweight": "bold"})


def pi1(p):
    p = p[np.isfinite(p)]
    lam = np.arange(0.4, 0.81, 0.05)
    return 1 - min(1, np.median([(p > l).mean() / (1 - l) for l in lam]))


def p_pi1(ax):
    rows = []
    for mod, d in (("ATAC", "caQTL"), ("RNA", "eQTL")):
        for ct in CTS + ["AllCells"]:
            f = f"{BASE}/{d}_Output_optCov/{ct}/allchr_cis_{TAG}.cis_qtl.qval.txt.gz"
            if os.path.exists(f):
                rows.append((ct, mod, pi1(pd.read_csv(f, sep="\t", usecols=["pval_beta"]).pval_beta.to_numpy())))
    d = pd.DataFrame(rows, columns=["ct", "mod", "pi1"])
    d.to_csv(f"{BASE}/QTL_pi1_optCov.tsv", sep="\t", index=False, float_format="%.4g")
    order = CTS + ["AllCells"]
    x = np.arange(len(order))
    for i, mod in enumerate(("ATAC", "RNA")):
        v = d[d["mod"] == mod].set_index("ct").reindex(order).pi1
        ax.bar(x + (i - 0.5) * 0.4, v, 0.4, color=MODCOL[mod], label={"ATAC": "caQTL", "RNA": "eQTL"}[mod])
    ax.set_xticks(x, [SHORT.get(c, c) for c in order], rotation=60, ha="right", fontsize=5.5)
    ax.set(ylabel="pi1 (permutation p)", title="Genome-wide QTL signal\n(17 donors; vasc. endo. 12)")
    ax.legend(frameon=False, fontsize=5.5)


def p_mash(ax):
    x = np.arange(len(CTS))
    for i, mod in enumerate(("ATAC", "RNA")):
        m = pd.read_csv(f"{BASE}/mash/{mod}_mash_summary.tsv", sep="\t").set_index("celltype").reindex(CTS)
        ax.bar(x + (i - 0.5) * 0.4, m.mash_lfsr05, 0.4, color=MODCOL[mod], label={"ATAC": "caQTL peaks", "RNA": "eGenes"}[mod])
        a = pd.read_csv(f"{BASE}/mash/{mod}_mash_summary.tsv", sep="\t").set_index("celltype").loc["any"]
        ax.text(0.98, 0.95 - 0.09 * i, f"any CT: {int(a.mash_lfsr05):,} {'caQTL peaks' if mod == 'ATAC' else 'eGenes'}",
                transform=ax.transAxes, ha="right", fontsize=5.5, color=MODCOL[mod])
    ax.set_xticks(x, [SHORT.get(c, c) for c in CTS], rotation=60, ha="right", fontsize=5.5)
    ax.set(ylabel="QTLs (mashR lfsr < 0.05)", title="Joint mapping across cell types", ylim=(0, 1750))


def p_brain(ax):
    p = pd.read_csv(f"{BASE}/BrainRep/pi1.tsv", sep="\t")
    p = p[(p.snp == "best_shared_snp") & p.set.isin(["mash", "mash_null"])]
    p = p[p.celltype.isin(CTS)]
    for j, (metric, lab) in enumerate((("pi1", "pi1"), ("sign_conc_p05", "sign concordance"))):
        for k, bc in enumerate(("neuron", "glia")):
            for s, col in (("mash", "#7B3294"), ("mash_null", "#BBBBBB")):
                v = p[(p.brain_ct == bc) & (p.set == s)][metric].dropna()
                xpos = j * 2.6 + k * 1.1 + (0.22 if s == "mash_null" else -0.22)
                ax.scatter(np.full(len(v), xpos) + np.random.default_rng(1).uniform(-0.08, 0.08, len(v)), v, s=5,
                           color=col, zorder=3, label=("mashR caQTL" if s == "mash" else "null peaks") if j == k == 0 else None)
                ax.plot([xpos - 0.15, xpos + 0.15], [v.median()] * 2, color="k", lw=1)
    ax.axhline(0.5, ls="--", lw=0.5, color="grey")
    ax.set_xticks([0, 1.1, 2.6, 3.7], ["pi1\nneuron", "pi1\nglia", "sign\nneuron", "sign\nglia"], fontsize=5.5)
    ax.set(ylim=(0, 0.85), ylabel="replication", title="Adult brain replication\n(one dot per HDC cell type)")
    ax.legend(frameon=False, fontsize=5.5, loc="upper left")


def p_sldsc(ax):
    d = pd.read_csv(f"{BASE}/GWAS/ldsc/sldsc_spec.tsv", sep="\t")
    d = d[d.primary & (d.family == "peakspec")]
    d["ct"] = d.annot.str.replace("peakspec_", "")
    # one row per trait label: keep the cohort with the smallest p (same phenotype in UKB_EUR / meta etc.)
    keep = d[d.FDR < 0.1].label.unique()
    s = d[d.label.isin(keep)].copy()
    s["lab"] = s.label.str.replace(", combined previous and current smoking", "").str.replace(" (BMI)", "", regex=False)
    s["lab"] = s.lab.str.replace("Body-mass index, inverse-rank normalized", "Body mass index").str[:34]
    best = s.loc[s.groupby(["lab", "ct"]).p_tau.idxmin()]
    M = best.pivot(index="lab", columns="ct", values="p_tau")
    Q = best.pivot(index="lab", columns="ct", values="FDR")
    cols = [c for c in ALLCT if c in M.columns]
    M, Q = -np.log10(M[cols]), Q[cols]
    M = M.loc[M.max(1).sort_values(ascending=False).index]
    Q = Q.loc[M.index]
    sns.heatmap(M.clip(upper=6), cmap="Reds", vmin=0, vmax=6, ax=ax, cbar_kws={"label": "-log10 p (tau)", "shrink": 0.7},
                xticklabels=[SHORT.get(c, c) for c in cols], yticklabels=True)
    yy, xx = np.where(Q.to_numpy() < 0.1)
    ax.scatter(xx + 0.5, yy + 0.5, s=6, color="k", marker="o")
    ax.tick_params(axis="x", rotation=60, labelsize=5.5)
    ax.tick_params(axis="y", labelsize=5.5)
    ax.set(xlabel="", ylabel="", title="S-LDSC: cell-type-specific peaks\n(dot: FDR < 10%)")


def p_loci(ax):
    o = pd.read_csv(f"{BASE}/GWAS/results/coloc_audit/overlap.tsv", sep="\t").set_index("level")
    r = o.loc["region"]
    vals = [r.perCT_only, r.shared, r.pooled_only]
    ax.barh([0], [vals[0]], color="#E6550D")
    ax.barh([0], [vals[1]], left=vals[0], color="#7B3294")
    ax.barh([0], [vals[2]], left=vals[0] + vals[1], color="#555555")
    for v, l, c in zip(vals, np.cumsum([0] + vals[:-1]), ("per-cell-type only", "both", "pooled only")):
        ax.text(l + v / 2, 0, f"{int(v)}", ha="center", va="center", color="w", fontsize=7, fontweight="bold")
        ax.text(l + v / 2, 0.55, c, ha="center", va="bottom", fontsize=5.5)
    ax.set(xlim=(0, sum(vals)), ylim=(-0.6, 1.2), yticks=[], xlabel="independent colocalized regions (PP.H4 > 0.5; lead variants > 500 kb apart)",
           title="Cell-type maps add trait loci")
    ax.spines["left"].set_visible(False)


def p_locus(ax):
    e = EX
    qdir = f"{BASE}/{'eQTL' if e['modality'] == 'RNA' else 'caQTL'}_Output_optCov/{e['celltype']}/"
    q = pq.read_table(f"{qdir}/chr{e['chr']}_cis_nominal_{TAG}.cis_qtl_pairs.{e['chr']}.parquet",
                      columns=["phenotype_id", "variant_id", "pval_nominal"],
                      filters=[("phenotype_id", "==", e["phenotype_id"])]).to_pandas()
    u = pd.read_parquet(f"{BASE}/GWAS/variant_universe.parquet", columns=["vidx", "chr", "pos38", "id"],
                        filters=[("chr", "==", e["chr"])])
    q = q.merge(u, left_on="variant_id", right_on="id")
    h = pd.read_csv(f"{BASE}/GWAS/results/coloc_H4gt0.5.tsv", sep="\t")
    h = h[(h.celltype == e["celltype"]) & (h.phenotype_id == e["phenotype_id"]) & h.trait_id.isin(e["traits"])]
    best = h.sort_values("PP.H4", ascending=False).iloc[0]
    g = pd.read_parquet(f"{BASE}/GWAS/harmonized/{best.trait_id}.parquet", filters=[("chr", "==", e["chr"])])
    if "pval" not in g.columns:
        from scipy.stats import norm
        g["pval"] = 2 * norm.sf(np.abs(g.beta / g.se))
    g = g.merge(u[["vidx", "pos38"]], on="vidx")
    g = g[(g.pos38 > e["lo"]) & (g.pos38 < e["hi"])]
    q = q[(q.pos38 > e["lo"]) & (q.pos38 < e["hi"])]
    lead = best.best_pos38
    ax.scatter(g.pos38 / 1e6, -np.log10(g.pval.astype(float)), s=4, color="#3B6FB6", rasterized=True)
    ax.scatter(q.pos38 / 1e6, np.log10(q.pval_nominal), s=4, color="#E6550D", rasterized=True)
    for d, sgn in ((g.assign(lp=-np.log10(g.pval.astype(float))), 1), (q.assign(lp=-np.log10(q.pval_nominal)), -1)):
        r = d[d.pos38 == lead]
        if len(r):
            ax.scatter(lead / 1e6, sgn * r.lp.iloc[0], s=30, marker="D", facecolor="none", edgecolor="k", lw=0.8, zorder=4)
    ax.axhline(0, color="k", lw=0.5)
    ax.axvline(lead / 1e6, color="grey", lw=0.5, ls="--")
    pk = [int(v) for v in e["phenotype_id"].split("_")[1:]]
    ax.axvspan(pk[0] / 1e6, pk[1] / 1e6, color="#E6550D", alpha=0.3, lw=0)
    yl = ax.get_ylim()
    ax.text(0.01, 0.97, f"{best.label} GWAS (FinnGen)", transform=ax.transAxes, va="top", fontsize=5.5, color="#3B6FB6")
    ax.text(0.01, 0.03, f"Neural crest caQTL ({e['phenotype_id']})", transform=ax.transAxes, va="bottom", fontsize=5.5,
            color="#E6550D")
    ax.text(lead / 1e6, yl[1] * 0.85, f" {best.best_snp}", fontsize=5.5)
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{abs(v):.0f}"))
    pp = ", ".join(f"{r.label.split(' ')[0].lower()} {r['PP.H4']:.2f}" for _, r in h.iterrows())
    ax.set(xlabel=f"chr{e['chr']} (Mb)", ylabel="-log10 p", title=f"11q13 coronary disease locus: PP.H4 {pp}")


def p_gtex(ax):
    t = pd.read_csv(f"{BASE}/GWAS/results/gtex_novelty_vs_background.tsv", sep="\t")
    t = t[(t.modality == "RNA") & (t.annotation == "gene_GTEx_eGene")].iloc[0]
    ax.bar([0, 1], [t.frac_coloc, t.frac_background], color=["#E6550D", "#BBBBBB"], width=0.6)
    for i, (v, n) in enumerate(((t.frac_coloc, t.n_coloc), (t.frac_background, t.n_background))):
        ax.text(i, v + 0.02, f"{v:.0%}\nn={int(n):,}", ha="center", fontsize=5.5)
    ax.set_xticks([0, 1], ["HDC-colocalized\neGenes", "genes at other\nGWAS leads"], fontsize=5.5)
    ax.set(ylim=(0, 1.2), ylabel="GTEx v8 eGene (any tissue)", title=f"Beyond adult eQTLs\n(OR {t.OR:.2f}, p = {t.p_fisher:.2g})")


PANELS = [("a", p_pi1), ("b", p_mash), ("c", p_brain), ("d", p_sldsc), ("e", p_loci), ("f", p_locus), ("g", p_gtex)]
fig = plt.figure(figsize=(7.5, 9))
gs = fig.add_gridspec(4, 6, hspace=1.1, wspace=1.8, height_ratios=[1, 1.5, 0.45, 1])
pos = [gs[0, 0:2], gs[0, 2:4], gs[0, 4:6], gs[1, 1:5], gs[2, 0:6], gs[3, 0:4], gs[3, 4:6]]
XS = {"a": -0.35, "b": -0.48, "c": -0.35, "d": -0.5, "e": -0.04, "f": -0.12, "g": -0.45}
for (t, fn), p in zip(PANELS, pos):
    ax = fig.add_subplot(p)
    try:
        fn(ax)
    except Exception as ex:  # noqa: BLE001
        ax.set_title(f"{t}: missing ({type(ex).__name__}: {ex})"[:80])
    ax.text(XS.get(t, -0.25), 1.12 if t != "e" else 1.6, t, transform=ax.transAxes, fontsize=10, fontweight="bold", va="top")
fig.savefig(f"{MAINF}/Figure2.pdf", bbox_inches="tight")
fig.savefig(f"{MAINF}/Figure2.png", dpi=200, bbox_inches="tight")
plt.close(fig)
for t, fn in PANELS:
    fig, ax = plt.subplots(figsize=(4.5, 3.2) if fn is not p_sldsc else (5, 4.5))
    try:
        fn(ax)
        fig.savefig(f"{F}/panels/fig2_{fn.__name__[2:]}.pdf", bbox_inches="tight")
    except Exception as ex:  # noqa: BLE001
        print("panel", t, "skipped:", ex)
    plt.close(fig)
print("Figure 2 written")
