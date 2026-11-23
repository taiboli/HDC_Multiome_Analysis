"""
Main-figure panels for the HDC peak / peak-gene link characterization, from the tables of 25, 25b (HOMER), 25c, 26
(and S-LDSC link batch when present). Panels whose inputs are missing are skipped.

Page 1 (peaks):  a genomic context, b ENCODE cCRE class (HDC vs shuffled), c active ChromHMM fraction by epigenome group,
                 d specificity deciles, e CT-specific peaks x matched epigenomes (15-state enhancers), f HOMER motifs
Page 2 (links):  g links per cell type (+ permuted null), h link rate vs TSS distance, i linked-peak enrichment,
                 j PCHi-C support, k caQTL-eQTL concordance, l cross-cell-type sharing (pi1), m linked genes / peaks
                 cell-type specificity, n pooled (AllCells) recovery of per-CT links; supp_sldsc = S-LDSC of linked peaks
Output: Figures/Main/Figure1.{pdf,png} (main: cCRE, matched-epigenome enhancers, motifs, links + permuted null,
        PCHi-C, caQTL->eQTL, pooled recovery, GAP43 locus from 27b), PeakAnnot/figures/FigureS_peaks_links.pdf (remaining panels; S1-S2 via 34) and
        one PDF per panel in PeakAnnot/figures/panels/
"""
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

BASE = os.environ["HDC_DIR"] + "/"
PA = f"{BASE}/PeakAnnot/"
R, F = f"{PA}/results/", f"{PA}/figures/"
MAINF = f"{BASE}/Figures/Main/"
os.makedirs(f"{F}/panels", exist_ok=True)
os.makedirs(MAINF, exist_ok=True)
CTS = ["iPSC", "Ectoderm_1", "Ectoderm_2", "Ectoderm_3", "Neural_crest", "Neuronal", "Retinal", "Mesoderm",
       "Vascular_endothelial", "Endoderm"]
COL = dict(zip(CTS, sns.color_palette("tab10", len(CTS))))
plt.rcParams.update({"font.size": 7, "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.titlesize": 8, "axes.titleweight": "bold"})


def rd(name):
    f = f"{R}/{name}"
    return pd.read_csv(f, sep="\t", keep_default_na=False, na_values=["", "NaN", "nan"]) if os.path.exists(f) else None


TAGX = {"c": -0.3, "i": -0.55, "j": -0.3, "l": -0.35, "m": -0.28, "n": -0.28}


def tag(ax, s):
    ax.text(TAGX.get(s, -0.18), 1.08, s, transform=ax.transAxes, fontsize=10, fontweight="bold", va="top")


# ---------------- page 1 ----------------
def p_context(ax):
    ov = rd("peaks_overview.tsv").set_index("set")
    cats = ["promoter", "exon", "intron", "distal"]
    x = np.arange(len(cats))
    for i, (s, c) in enumerate((("HDC peaks", "#3B6FB6"), ("shuffled", "#BBBBBB"))):
        ax.bar(x + (i - 0.5) * 0.38, [ov.loc[s, f"frac_{k}"] for k in cats], 0.38, color=c, label=s)
    ax.set_xticks(x, cats)
    ax.set(ylabel="fraction of peaks", title=f"{int(ov.loc['HDC peaks', 'n']):,} HDC peaks")
    ax.legend(frameon=False)


def p_ccre(ax):
    ov = rd("peaks_overview.tsv").set_index("set")
    cl = ["PLS", "pELS", "dELS", "CA-H3K4me3", "CA-CTCF", "CA-TF", "CA", "TF", "none"]
    x = np.arange(len(cl))
    for i, (s, c) in enumerate((("HDC peaks", "#3B6FB6"), ("shuffled", "#BBBBBB"))):
        ax.bar(x + (i - 0.5) * 0.38, [ov.loc[s, f"ccre_{k}"] for k in cl], 0.38, color=c)
    ax.set_xticks(x, cl, rotation=60, ha="right")
    ax.set(ylabel="fraction of peaks",
           title=f"ENCODE cCRE v4\n{ov.loc['HDC peaks', 'frac_ccre']:.0%} vs {ov.loc['shuffled', 'frac_ccre']:.0%} (shuffled)")


def p_active(ax):
    a = rd("chromhmm18_active_frac_by_epigenome.tsv")
    order = a.groupby("group").frac_hdc.median().sort_values(ascending=False).index
    sns.boxplot(data=a, y="group", x="frac_hdc", order=order, ax=ax, color="#9ecae1", fliersize=0, linewidth=0.6)
    sns.stripplot(data=a, y="group", x="frac_hdc", order=order, ax=ax, color="k", size=2)
    ax.set(xlabel="fraction of HDC peaks in active states\n(Roadmap 18-state, per epigenome)", ylabel="",
           title="Active chromatin by tissue")
    ax.tick_params(axis="y", labelsize=5.5)


def p_deciles(ax):
    d = rd("specificity_deciles.tsv")
    for c, lab, col in (("promoter", "promoter", "#E41A1C"), ("esc_active", "active in ESC (H1/H9)", "#377EB8"),
                        ("esc_poised", "bivalent/Polycomb in ESC", "#984EA3"), ("phastcons", "phastCons element", "#4DAF4A")):
        ax.plot(np.arange(1, 11), d[c], "o-", ms=3, color=col, label=lab)
    ax.set(xlabel="cell-type specificity decile (tau)", ylabel="fraction of peaks", title="Specificity deciles")
    ax.legend(frameon=False, fontsize=5, loc="upper center", bbox_to_anchor=(0.5, -0.32), ncol=2)


def p_epi(ax):
    ce = rd("ctspec_peaks_epigenome15_enhancer_enrichment.tsv")
    sel = ["E003", "E008", "E021", "E004", "E007", "E009", "E010", "E011", "E012", "E013", "E053", "E081", "E082", "E070",
           "E059", "E061", "E092", "E084", "E085", "E066", "E013", "E122", "E057"]
    sel = list(dict.fromkeys(e for e in sel if e in set(ce.eid)))
    H = ce.pivot(index="celltype", columns="eid", values="log2_enr").reindex(index=[c for c in CTS if c in set(ce.celltype)], columns=sel)
    nm = ce.drop_duplicates("eid").set_index("eid").name
    sns.heatmap(H, cmap="RdBu_r", center=0, vmin=-1.5, vmax=1.5, ax=ax, yticklabels=True,
                cbar_kws={"label": "log2 enrichment\n(EnhG+Enh, 15-state)", "shrink": 0.7})
    ax.set_xticks(np.arange(len(sel)) + 0.5, [nm[e][:24] for e in sel], rotation=60, ha="right", fontsize=5)
    ax.tick_params(axis="y", labelsize=5.5)
    ax.set(xlabel="", ylabel="", title="CT-specific peaks: enhancers of matched epigenomes")


def p_homer(ax, wide=False):
    rows = []
    for f in sorted(glob.glob(f"{PA}/homer/out/peakspec_*/knownResults.txt")):
        ct = f.split("/")[-2].replace("peakspec_", "")
        d = pd.read_csv(f, sep="\t")
        d["motif"] = d.iloc[:, 0].str.split("/").str[0]
        d["mlog10p"] = -d.iloc[:, 3] / np.log(10)
        d["fold"] = d.iloc[:, 6].str.rstrip("%").astype(float) / d.iloc[:, 8].str.rstrip("%").astype(float).clip(lower=0.1)
        rows.append(d[["motif", "mlog10p", "fold"]].drop_duplicates("motif").assign(celltype=ct))
    if not rows:
        return
    d = pd.concat(rows)
    top = []
    for ct in CTS:
        x = d[d.celltype == ct].sort_values("mlog10p", ascending=False)
        top += [m for m in x.motif if m not in top][:3]
    M = d[d.motif.isin(top)].pivot(index="motif", columns="celltype", values="mlog10p").reindex(index=top, columns=[c for c in CTS if c in set(d.celltype)])
    if wide:  # cell types as rows, motifs as columns (main figure)
        sns.heatmap(M.T.clip(upper=200), cmap="magma_r", ax=ax, xticklabels=True, yticklabels=True,
                    cbar_kws={"label": "-log10 p\n(capped 200)", "shrink": 0.8})
        ax.set_xticklabels([m.get_text().split("(")[0].replace("-TCF", "") for m in ax.get_xticklabels()], rotation=60,
                           ha="right", fontsize=5)
    else:
        sns.heatmap(M.clip(upper=200), cmap="magma_r", ax=ax, yticklabels=True, cbar_kws={"label": "-log10 p (capped 200)", "shrink": 0.6})
    ax.set(xlabel="", ylabel="", title="Known motifs (HOMER) in CT-specific peaks")
    ax.tick_params(axis="y", labelsize=5.5)


def p_homer_wide(ax):
    p_homer(ax, wide=True)


def p_locus(ax, gene="GAP43"):
    """placeholder axis replaced by a native locus drawing (27b_locus.draw_locus) at the same grid position"""
    import importlib.util
    spec = importlib.util.spec_from_file_location("locus", f"{BASE}/scripts/27b_locus.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    fig, ss = ax.figure, ax.get_subplotspec()
    ax.remove()
    return m.draw_locus(gene, fig=fig, spec=ss, title_fs=6.5)


# ---------------- page 2 ----------------
def p_links(ax):
    s = rd("link_summary.tsv")
    s = s[s.celltype.isin(CTS)].set_index("celltype").reindex([c for c in CTS if c in set(s.celltype)])
    x = np.arange(len(s))
    ax.bar(x, s.links_pos, color=[COL[c] for c in s.index])
    for i, (lp, g) in enumerate(zip(s.links_pos, s.linked_genes)):
        ax.text(i, lp, f"{int(lp)}", ha="center", va="bottom", fontsize=5, rotation=90 if lp < 1000 else 0)
    pm = rd("scmultimap_perm_after_filter.tsv")
    t = "scMultiMap links (FDR 10%, positive)"
    if pm is not None:
        t += "\npermuted null: " + ", ".join(f"{r.set.replace('_perm', '')} {int(r.n_padj_lt_fdr)}" for r in pm.itertuples())
    ax.set_xticks(x, s.index, rotation=60, ha="right", fontsize=5.5)
    ax.set(ylabel="positive links", title=t)
    ax.set_ylim(0, s.links_pos.max() * 1.25)


def p_distance(ax):
    d = rd("link_distance.tsv")
    lab = ["0-2k", "2-10k", "10-50k", "50-100k", "100-250k", "250-500k"]
    for ct in CTS:
        x = d[d.celltype == ct].set_index("bin").reindex(lab)
        if x.rate.notna().any() and ct != "Vascular_endothelial":
            ax.plot(np.arange(len(lab)), x.rate * 1000, "o-", ms=2.5, color=COL[ct], label=ct)
    ax.set_xticks(np.arange(len(lab)), lab, rotation=45, ha="right")
    ax.set(xlabel="peak - TSS distance", ylabel="links per 1,000 tested pairs", title="Link rate decays with distance")
    ax.legend(frameon=False, fontsize=5, ncol=2)


def p_enrich(ax):
    e = rd("link_peak_enrichment.tsv")
    e = e[e.celltype.isin(CTS) & (e.celltype != "Vascular_endothelial")]
    ann = ["CT_specific_peak", "enhancer_matched", "cCRE_dELS", "phastCons", "promoter"]
    lab = {"CT_specific_peak": "CT-specific\npeak", "enhancer_matched": "enhancer,\nmatched epi.",
           "cCRE_dELS": "ENCODE\ndELS", "phastCons": "phastCons\nelement", "promoter": "promoter"}
    cts = [c for c in CTS if c in set(e.celltype)]
    for j, a in enumerate(ann):
        for i, ct in enumerate(cts):
            r = e[(e.celltype == ct) & (e.annotation == a)]
            if len(r):
                r = r.iloc[0]
                y = j + (i - len(cts) / 2) * 0.09
                ax.plot([r.lo, r.hi], [y, y], color=COL[ct], lw=0.8)
                ax.plot(r.adj_OR, y, "o", color=COL[ct], ms=2.5)
    ax.axvline(1, color="grey", lw=0.6, ls="--")
    ax.set_xscale("log")
    ax.set_yticks(range(len(ann)), [lab[a] for a in ann])
    ax.set(xlabel="adjusted OR (linked vs unlinked tested peaks)", title="Linked-peak annotation enrichment")
    ax.legend(handles=[plt.Line2D([], [], marker="o", ls="", ms=3, color=COL[c], label=c) for c in cts], frameon=False,
              fontsize=4.5, loc="upper center", bbox_to_anchor=(0.4, -0.4), ncol=5, handletextpad=0.1, columnspacing=0.5)
    ax.set_ylim(-0.6, len(ann) - 0.4)


def p_pchic(ax):
    v = rd("link_pair_validation.tsv")
    v = v[(v.support == "PCHiC_ES_derived") & (v.celltype != "Vascular_endothelial")]
    v = v.set_index("celltype").reindex([c for c in CTS + ["AllCells"] if c in set(v.celltype)])
    x = np.arange(len(v))
    ax.bar(x - 0.2, v.rate_link, 0.4, color=[COL.get(c, "k") for c in v.index], label="linked")
    ax.bar(x + 0.2, v.rate_null, 0.4, color="#CCCCCC", label="unlinked")
    for i, r in enumerate(v.itertuples()):
        ax.text(i, max(r.rate_link, r.rate_null) + 0.004, f"{r.MH_OR:.1f}", ha="center", fontsize=5)
    ax.set_xticks(x, v.index, rotation=45, ha="right")
    ax.set_ylim(0, v.rate_link.max() * 1.35)
    ax.set(ylabel="frac. pairs with contact",
           title=f"3D contact support (MH OR above bars;\nall p < {v.p.max():.0e})")
    ax.legend(frameon=False, fontsize=5.5, loc="upper right", ncol=2)


def p_qtl(ax):
    q = rd("link_qtl_validation.tsv")
    q = q[q.celltype == "pooled"].set_index("class").rename(index={"null": "unlinked"})
    cls = [c for c in ["link_pos", "unlinked"] if c in q.index]
    x = np.arange(len(cls))
    ax.bar(x - 0.2, q.loc[cls, "frac_eqtl_p05"], 0.4, color="#3B6FB6", label="eQTL p<0.05 at caQTL lead")
    ax.bar(x + 0.2, q.loc[cls, "concordance"], 0.4, color="#E6550D", label="direction concordant")
    for i, c in enumerate(cls):
        ax.text(i + 0.2, q.loc[c, "concordance"] + 0.02, f"{q.loc[c, 'concordance']:.0%}\nn={int(q.loc[c, 'n_eqtl_p05'])}",
                ha="center", fontsize=5)
        ax.text(i - 0.2, q.loc[c, "frac_eqtl_p05"] + 0.02, f"{q.loc[c, 'frac_eqtl_p05']:.0%}", ha="center", fontsize=5)
    ax.axhline(0.5, color="grey", lw=0.5, ls="--")
    ax.set_xticks(x, [{"link_pos": "positive\nlinks", "unlinked": "unlinked"}[c] for c in cls])
    p = q.loc["link_pos", "binom_p"] if "link_pos" in q.index else np.nan
    ax.set(ylim=(0, 1.1), ylabel="fraction", title=f"caQTL lead -> eQTL of linked gene\n(concordance p={p:.1g})")
    ax.legend(frameon=False, fontsize=5, loc="upper center", bbox_to_anchor=(0.5, -0.3), ncol=1)


def p_sharing(ax):
    s = rd("link_sharing.tsv")
    col = "pi1" if "pi1" in s.columns else "rep_rate"
    M = s.pivot(index="discovery", columns="replication", values=col)
    order = [c for c in CTS if c in M.index and c != "Vascular_endothelial"]
    M = M.reindex(index=order, columns=order + (["AllCells"] if "AllCells" in M.columns else []))
    sns.heatmap(M, cmap="viridis", vmin=0, vmax=1, annot=True, fmt=".2f", annot_kws={"fontsize": 4.5}, ax=ax, yticklabels=True,
                cbar_kws={"label": "pi1" if col == "pi1" else "replication rate", "shrink": 0.6})
    ax.set(xlabel="replication cell type", ylabel="", title="Link sharing across cell types (rows: discovery)")


def p_spec(ax):
    g = rd("link_specificity.tsv")
    g = g[g.celltype != "Vascular_endothelial"].set_index("celltype")
    g = g.reindex([c for c in CTS if c in g.index])
    x = np.arange(len(g))
    ax.bar(x - 0.3, g.frac_linked_genes_ctspec, 0.2, color="#3B6FB6", label="linked genes")
    ax.bar(x - 0.1, g.frac_tested_genes_ctspec, 0.2, color="#9ecae1", label="tested genes")
    ax.bar(x + 0.1, g.frac_linked_peaks_ctspec, 0.2, color="#E6550D", label="linked peaks")
    ax.bar(x + 0.3, g.frac_tested_peaks_ctspec, 0.2, color="#fdae6b", label="tested peaks")
    ax.set_xticks(x, g.index, rotation=45, ha="right")
    ax.set(ylabel="fraction cell-type specific", title="CT-specific peaks and genes among links")
    ax.legend(frameon=False, fontsize=5, ncol=2)


def p_sldsc(ax):
    f = f"{BASE}/GWAS/ldsc/sldsc_link.tsv"
    if not os.path.exists(f):
        ax.axis("off")
        return
    d = pd.read_csv(f, sep="\t")
    d = d[d.primary & (d.family == "speclinked")]
    d["ct"] = d.annot.str.replace("speclinked_", "")
    agg = d.groupby("ct").agg(med_tau_star=("tau_star", "median"), n_fdr=("FDR", lambda x: int((x < 0.1).sum())),
                              n_nom=("p_tau", lambda x: int((x < 0.05).sum())), n=("p_tau", "size"))
    agg = agg.reindex([c for c in CTS if c in agg.index])
    x = np.arange(len(agg))
    ax.bar(x, agg.n_nom / agg.n, color=[COL[c] for c in agg.index])
    ax.axhline(0.05, color="grey", ls="--", lw=0.6)
    for i, r in enumerate(agg.itertuples()):
        ax.text(i, r.n_nom / r.n, f"{r.n_fdr}", ha="center", va="bottom", fontsize=5)
    ax.set_xticks(x, agg.index, rotation=45, ha="right")
    ax.set(ylabel="fraction of traits with tau p<0.05\n(label: FDR<10%)",
           title="S-LDSC: linked CT-specific peaks\nbeyond all CT-specific peaks")


def p_pooled(ax):
    r = rd("link_allcells_recovery.tsv")
    s = rd("link_summary.tsv").set_index("celltype")
    r["lab"] = r.n_ct.astype(str)
    r = pd.concat([r[r.n_ct < 4], pd.DataFrame([{"lab": "4+", "n": r[r.n_ct >= 4].n.sum(),
        "tested": np.average(r[r.n_ct >= 4].tested, weights=r[r.n_ct >= 4].n),
        "recovered": np.average(r[r.n_ct >= 4].recovered, weights=r[r.n_ct >= 4].n)}])])
    x = np.arange(len(r))
    ax.bar(x - 0.2, r.tested, 0.4, color="#BBBBBB", label="tested in pooled")
    ax.bar(x + 0.2, r.recovered, 0.4, color="k", label="pooled FDR<10%")
    for i, q in enumerate(r.itertuples()):
        ax.text(i, max(q.tested, q.recovered) + 0.03, f"n={int(q.n):,}", ha="center", fontsize=5)
    ax.set_xticks(x, r.lab)
    a = s.loc["AllCells"]
    ax.set(ylim=(0, 1.45), xlabel="cell types in which the link is significant", ylabel="fraction of per-CT positive links",
           title=f"Pooled (all cells) analysis misses CT-specific links\n(pooled: {int(a.links):,} links, "
                 f"{a.frac_pos:.0%} positive)")
    ax.legend(frameon=False, fontsize=5.5, loc="upper left", ncol=2)


MAIN = [("a", p_ccre), ("b", p_epi), ("c", p_homer_wide), ("d", p_links), ("e", p_pchic), ("f", p_qtl), ("g", p_pooled),
        ("h", p_locus)]
SUPP_PEAKS = [("a", p_context), ("b", p_active), ("c", p_deciles), ("d", p_homer)]
SUPP_LINKS = [("a", p_distance), ("b", p_enrich), ("c", p_sharing), ("d", p_spec), ("e", p_sldsc)]
TAGX.clear()
TAGX.update({"b": -0.3, "c": -0.12, "d": -0.3, "e": -0.3, "f": -0.3, "g": -0.3, "h": -0.02})


def draw(fig, panels, pos, xs=None, titles=None, ys=None):
    for (t, fn), p in zip(panels, pos):
        ax = fig.add_subplot(p)
        try:
            ax = fn(ax) or ax
            if titles and t in titles:
                ax.set_title(titles[t], fontsize=7)
            ax.text((xs or {}).get(t, TAGX.get(t, -0.18)), (ys or {}).get(t, 1.08), t, transform=ax.transAxes, fontsize=10, fontweight="bold", va="top")
        except Exception as e:  # noqa: BLE001
            ax.set_title(f"{t}: missing ({type(e).__name__})")


# main Figure 1: the strongest positive results
fig = plt.figure(figsize=(7.5, 10.5))
gs = fig.add_gridspec(4, 6, hspace=1.25, wspace=1.9, height_ratios=[1.1, 0.95, 1, 1.5])
draw(fig, MAIN, [gs[0, 0:2], gs[0, 2:6], gs[1, 0:6], gs[2, 0:2], gs[2, 2:4], gs[2, 4:6], gs[3, 0:2], gs[3, 2:6]],
     titles={"d": "Positive links (FDR 10%)\npermuted: Neuronal 41, iPSC 165", "e": "PCHi-C contact\n(MH OR above bars)",
             "f": "caQTL lead -> eQTL of\nlinked gene (p=2e-4)", "g": "Pooled analysis misses\nCT-specific links"},
     xs={"e": -0.4, "g": -0.35, "h": -0.3}, ys={"c": 1.35, "h": 1.0})
fig.savefig(f"{MAINF}/Figure1.pdf", bbox_inches="tight")
fig.savefig(f"{MAINF}/Figure1.png", dpi=200, bbox_inches="tight")
plt.close(fig)

# supplementary: remaining peak and link panels
with PdfPages(f"{F}/FigureS_peaks_links.pdf") as pdf:
    fig = plt.figure(figsize=(7.5, 8))
    gs = fig.add_gridspec(2, 6, hspace=0.4, wspace=1.6, height_ratios=[1, 2.2])
    draw(fig, SUPP_PEAKS, [gs[0, 0:2], gs[0, 2:4], gs[0, 4:6], gs[1, 1:5]], {"b": -0.35})
    pdf.savefig(fig, bbox_inches="tight")
    fig.savefig(f"{F}/FigureS_peaks.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    fig = plt.figure(figsize=(7.5, 8))
    gs = fig.add_gridspec(3, 6, hspace=1.2, wspace=2.6)
    draw(fig, SUPP_LINKS, [gs[0, 0:3], gs[0, 3:6], gs[1, 0:6], gs[2, 0:3], gs[2, 3:6]], {"b": -0.35, "c": -0.2, "d": -0.28, "e": -0.28})
    pdf.savefig(fig, bbox_inches="tight")
    fig.savefig(f"{F}/FigureS_links.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
for t, fn in MAIN + SUPP_PEAKS + SUPP_LINKS:
    if fn is p_locus:
        continue
    fig, ax = plt.subplots(figsize=(4.5, 3.4) if fn not in (p_epi, p_homer, p_homer_wide, p_sharing) else (7, 4.5))
    try:
        fn(ax)
        fig.savefig(f"{F}/panels/{fn.__name__[2:]}.pdf", bbox_inches="tight")
    except Exception as e:  # noqa: BLE001
        print("panel", fn.__name__, "skipped:", e)
    plt.close(fig)
print("figures written")
