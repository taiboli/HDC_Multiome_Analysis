"""
Assemble the supplementary figures S1-S13 into Figures/Supplementary/FigureS<n>.{pdf,png}.

Drawn here from tables (no earlier figure):
  S4  scMultiMap calibration (PeakAnnot/results/scmultimap_calibration.tsv)
  S6  QTL covariate optimization (CovOptimization/all_summary.tsv)
  S7  mashR pairwise sharing (mash/<MOD>_mash_sharing.tsv)
  S8  adult brain replication per cell type (BrainRep/pi1.tsv, overlap_enrichment.tsv)
  S9  full S-LDSC: peak, gene, QTL and link sets (GWAS/ldsc/sldsc_{spec,qtl,link}.tsv)
  S10 colocalization audit (GWAS/results/coloc_audit/hits_audited.tsv, R_coloc_check.txt)
Taken from earlier scripts (run those first):
  S1, S2   pages 1-2 of PeakAnnot/figures/FigureS_peaks_links.pdf (27_figures.py), vector
  S3       fig_peaks_ctspec_epigenome{,15}.pdf (25, 25c), rasterized at 300 dpi and stacked as a / b
  S5       panels/locus_{TDGF1,GATA6}.png (27b_locus.py), stacked as a / b
  S11-S13  pages 1-3 of PeakAnnot/figures/FigureS_tfnet.pdf (32_figure3.py), vector
Needs ghostscript (gs) for PDF page extraction and PNG previews.
Legends: Figures/Figure_legends.md (written by hand; keep in sync with the panels).
"""
import os
import re
import shutil
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402

BASE = os.environ["HDC_DIR"] + "/"
PF = f"{BASE}/PeakAnnot/figures/"
OUT = f"{BASE}/Figures/Supplementary/"
TMP = f"{OUT}/.tmp/"
os.makedirs(TMP, exist_ok=True)
CTS = ["iPSC", "Ectoderm_1", "Ectoderm_2", "Ectoderm_3", "Neural_crest", "Neuronal", "Retinal", "Mesoderm", "Endoderm",
       "Vascular_endothelial", "Unknown"]
COL = dict(zip(CTS, sns.color_palette("tab10", 10).as_hex() + ["#555555"]))
plt.rcParams.update({"font.size": 7, "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.titlesize": 7, "axes.titleweight": "bold"})


def tag(ax, s, x=-0.18, y=1.1):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=10, fontweight="bold", va="top")


def save(fig, n):
    fig.savefig(f"{OUT}/FigureS{n}.pdf", bbox_inches="tight")
    fig.savefig(f"{OUT}/FigureS{n}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def gs(*args):
    subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-dSAFER", *args], check=True)


def page(src, p, n, png):
    """Copy page p of a vector PDF to FigureS<n>.pdf, with the matching PNG written by the source script."""
    gs("-sDEVICE=pdfwrite", f"-dFirstPage={p}", f"-dLastPage={p}", f"-sOutputFile={OUT}/FigureS{n}.pdf", src)
    shutil.copy(png, f"{OUT}/FigureS{n}.png")


def raster(pdf):
    png = f"{TMP}/{os.path.basename(pdf)[:-4]}.png"
    gs("-sDEVICE=png16m", "-r300", "-dTextAlphaBits=4", "-dGraphicsAlphaBits=4", f"-sOutputFile={png}", pdf)
    return png


def stack(pngs, n, width=7.5):
    """Stack images vertically as panels a, b, ... (heights follow each image's aspect)."""
    ims = [plt.imread(p) for p in pngs]
    hs = [width * i.shape[0] / i.shape[1] for i in ims]
    fig, axs = plt.subplots(len(ims), 1, figsize=(width, sum(hs) + 0.3 * len(ims)), gridspec_kw={"height_ratios": hs})
    for k, (ax, im) in enumerate(zip(np.atleast_1d(axs), ims)):
        ax.imshow(im, interpolation="lanczos")
        ax.axis("off")
        tag(ax, "ab"[k], x=-0.02, y=1.04)
    fig.savefig(f"{OUT}/FigureS{n}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{OUT}/FigureS{n}.png", bbox_inches="tight", dpi=200)
    plt.close(fig)


# ---------------- S4: scMultiMap calibration ----------------
def s4():
    d = pd.read_csv(f"{BASE}/PeakAnnot/results/scmultimap_calibration.tsv", sep="\t")
    fig, axs = plt.subplots(1, 3, figsize=(7.5, 2.4), gridspec_kw={"width_ratios": [1, 1, 1.1], "wspace": 0.75})
    for k, ct in enumerate(("iPSC", "Neuronal")):
        g = d[d.ct == ct]
        r = g.pivot(index="gene_det", columns="peak_det", values="perm/real_fdr10")
        n = g.pivot(index="gene_det", columns="peak_det", values="real_fdr10")
        ann = r.applymap(lambda v: f"{v:.2f}") + "\n" + n.applymap(lambda v: f"{v:,}")
        sns.heatmap(r, ax=axs[k], cmap="Reds", vmin=0, vmax=0.6, annot=ann, fmt="", annot_kws={"fontsize": 5},
                    cbar=k == 1, cbar_kws={"label": "permuted / real links", "shrink": 0.8})
        i, j = list(r.index).index(0.1), list(r.columns).index(0.02)
        axs[k].add_patch(plt.Rectangle((j, i), 1, 1, fill=False, lw=1.5, ec="k"))
        axs[k].set(xlabel="peak detection (fraction of cells)", ylabel="gene detection (fraction of cells)", title=ct)
        axs[k].invert_yaxis()
        tag(axs[k], "ab"[k], x=-0.3)
    ax = axs[2]
    for ct, m in (("iPSC", "o"), ("Neuronal", "s")):
        g = d[d.ct == ct]
        for gd, c in zip(sorted(g.gene_det.unique()), sns.color_palette("viridis", 4)):
            x = g[g.gene_det == gd].sort_values("peak_det")
            ax.plot(x.real_fdr10, x.perm_fdr10 + 0.5, m + "-", ms=3, color=c, lw=0.8,
                    label=f"gene {gd:g}" if ct == "iPSC" else None)
    ax.plot([300, 2e4], [300, 2e4], ls="--", color="grey", lw=0.6)
    ax.set(xscale="log", yscale="log", xlabel="real links (FDR 10%)", ylabel="permuted links + 0.5",
           title="Real vs permuted links\n(● iPSC, ■ Neuronal)")
    ax.legend(frameon=False, fontsize=5, title="detection", title_fontsize=5)
    tag(ax, "c", x=-0.3)
    save(fig, 4)


# ---------------- S6: covariate optimization ----------------
def s6():
    c = pd.read_csv(f"{BASE}/CovOptimization/all_summary.tsv", sep="\t")
    c["ct"] = c.Set.str.rsplit("_", n=1).str[0]
    c["mod"] = c.Set.str.rsplit("_", n=1).str[1]
    fig = plt.figure(figsize=(7.5, 4.6))
    g = fig.add_gridspec(2, 6, hspace=0.6, wspace=1.4)
    a1 = c[c.Stage == "A1"]
    for k, (mod, sel) in enumerate((("RNA", 2), ("ATAC", 3))):
        ax = fig.add_subplot(g[0, 3 * k:3 * k + 3])
        for ct, x in a1[a1["mod"] == mod].groupby("ct"):
            x = x.sort_values("nPhenoPC")
            ax.plot(x.nPhenoPC, x.pi1, "o-", ms=2.5, lw=0.8, color=COL.get(ct, "k"), label=ct)
        ax.axvline(sel, ls="--", color="grey", lw=0.6)
        ax.set(xlabel=f"{'expression' if mod == 'RNA' else 'accessibility'} PCs", ylabel="π₁ (held-out chromosomes)",
               title=f"{mod}: phenotype PCs (no other covariates)")
        if k == 1:
            ax.legend(frameon=False, fontsize=5, ncol=2, loc="upper left", bbox_to_anchor=(1.0, 1.0))
        tag(ax, "ab"[k])
    # A2: sex and genotype PCs on top of each set's best A1 phenotype PCs
    a2 = c[c.Stage == "A2"].copy()
    best1 = a1.loc[a1.groupby("Set").pi1.idxmax()].set_index("Set").pi1
    a2["dpi1"] = a2.pi1 - a2.Set.map(best1)
    a2["cfg"] = a2.apply(lambda r: f"{'sex+' if r.Sex else ''}{r.nGenoPC} gPC", axis=1)
    order = [f"{s}{n} gPC" for s in ("", "sex+") for n in (0, 1, 2, 3, 5)]
    order = [o for o in order if o in set(a2.cfg)]
    ax = fig.add_subplot(g[1, 0:2])
    sns.stripplot(a2, x="cfg", y="dpi1", hue="mod", order=order, dodge=True, size=2.5, ax=ax, palette=["#2C7BB6", "#D7191C"],
                  hue_order=["RNA", "ATAC"])
    ax.axhline(0, color="grey", lw=0.6)
    ax.set_xticks(range(len(order)), order, rotation=60, ha="right", fontsize=5.5)
    ax.set(xlabel="", ylabel="Δπ₁ vs best phenotype PCs", title="Sex and genotype PCs")
    ax.legend(frameon=False, fontsize=5)
    tag(ax, "c", x=-0.35)
    # B: window and MAF, at each set's configuration with the most B runs
    # (the MAF 0 / 1 Mb reference point of each configuration comes from stage A1/A2)
    c["cfg"] = list(zip(c.Sex, c.nGenoPC, c.nPhenoPC))
    keep = c[c.Stage == "B"].groupby("Set").cfg.agg(lambda s: s.value_counts().index[0])
    b = c[c.cfg == c.Set.map(keep)].drop_duplicates(["Set", "maf", "window"])
    for k, (var, fix, fv, lab) in enumerate((("window", "maf", 0.0, "cis window (kb)"), ("maf", "window", 1_000_000, "MAF filter"))):
        ax = fig.add_subplot(g[1, 2 + 2 * k:4 + 2 * k])
        x = b[b[fix] == fv]
        for (s, mod), y in x.groupby(["Set", "mod"]):
            y = y.sort_values(var)
            ref = y.pi1.iloc[-1] if var == "window" else y.pi1.iloc[0]
            ax.plot(y[var] / (1000 if var == "window" else 1), y.pi1 - ref, "-", lw=0.7, alpha=0.8,
                    color="#2C7BB6" if mod == "RNA" else "#D7191C")
        ax.axhline(0, color="grey", lw=0.6)
        if var == "window":
            ax.set_xscale("log")
            ax.set_xticks([25, 50, 100, 250, 500, 1000], ["25", "50", "100", "250", "500", "1000"])
            ax.minorticks_off()
        ax.set(xlabel=lab, ylabel=f"Δπ₁ vs {'1 Mb' if var == 'window' else 'no filter'}",
               title="Cis window (MAF unfiltered)" if var == "window" else "MAF filter (1 Mb window)")
        if k == 0:
            ax.legend(handles=[plt.Line2D([], [], color=c, lw=1, label=m) for m, c in (("RNA", "#2C7BB6"), ("ATAC", "#D7191C"))],
                      frameon=False, fontsize=5, loc="lower left")
        tag(ax, "de"[k], x=-0.35)
    save(fig, 6)


# ---------------- S7: mashR sharing ----------------
def s7():
    fig, axs = plt.subplots(1, 3, figsize=(7.5, 3.2), gridspec_kw={"wspace": 0.3, "width_ratios": [1, 1, 0.04]})
    for k, mod in enumerate(("ATAC", "RNA")):
        m = pd.read_csv(f"{BASE}/mash/{mod}_mash_sharing.tsv", sep="\t", index_col=0)
        o = [c for c in CTS if c in m.index]
        m = m.loc[o, o]
        sns.heatmap(m, ax=axs[k], cmap="viridis", vmin=0.5, vmax=1, annot=True, fmt=".2f", annot_kws={"fontsize": 4.5},
                    square=True, cbar=k == 1, cbar_ax=axs[2] if k == 1 else None,
                    cbar_kws={"label": "shared (same sign, within 2-fold)"})
        axs[k].set(xlabel="", ylabel="", title=f"{'caQTL' if mod == 'ATAC' else 'eQTL'} sharing (mashR)")
        axs[k].set_xticklabels([t.get_text().replace("Vascular_endothelial", "Vasc_endo") for t in axs[k].get_xticklabels()], fontsize=5.5)
        axs[k].set_yticklabels([t.get_text().replace("Vascular_endothelial", "Vasc_endo") for t in axs[k].get_yticklabels()], fontsize=5.5)
        tag(axs[k], "ab"[k], x=-0.35)
    save(fig, 7)


# ---------------- S8: brain replication per cell type ----------------
def s8():
    p = pd.read_csv(f"{BASE}/BrainRep/pi1.tsv", sep="\t")
    p = p[p.snp == "best_shared_snp"]
    o = pd.read_csv(f"{BASE}/BrainRep/overlap_enrichment.tsv", sep="\t")
    cts = [c for c in CTS if c in set(p.celltype)]
    sets = {"mash": ("mashR caQTL", "#7B3294"), "mash_specific": ("CT-specific mashR", "#E66101"),
            "perm_p01": ("permutation p < 0.01", "#5E8C31"), "mash_null": ("null peaks", "#BBBBBB")}
    fig, axs = plt.subplots(3, 2, figsize=(7.5, 6.6), gridspec_kw={"hspace": 0.75, "wspace": 0.25})
    for j, bc in enumerate(("neuron", "glia")):
        for i, (metric, lab) in enumerate((("pi1", "π₁"), ("sign_conc_p05", "sign concordance (p < 0.05)"))):
            ax = axs[i, j]
            # CT-specific mashR peaks are too few (1-27 pairs) for pi1; shown only in e-f
            for k, s in enumerate(("mash", "mash_null")):
                v = p[(p.brain_ct == bc) & (p.set == s)].set_index("celltype").reindex(cts)
                x = np.arange(len(cts)) + (k - 0.5) * 0.36
                ax.bar(x, v[metric], 0.36, color=sets[s][1], label=sets[s][0])
                if metric == "pi1":
                    ax.errorbar(x, v.pi1, yerr=[v.pi1 - v.pi1_lo, v.pi1_hi - v.pi1], fmt="none", ecolor="k", lw=0.5)
            if metric != "pi1":
                ax.axhline(0.5, ls="--", lw=0.6, color="grey")
            ax.set_xticks(range(len(cts)), [c.replace("Vascular_endothelial", "Vasc_endo") for c in cts], rotation=45, ha="right")
            ax.set(ylabel=lab, title=f"Adult {bc} caQTL: {lab.split(' ')[0]}", ylim=(0, 1.05))
            if i == j == 0:
                ax.legend(frameon=False, fontsize=5, ncol=2, loc="upper left")
            tag(ax, "abcd"[i * 2 + j], x=-0.12)
        ax = axs[2, j]
        for k, s in enumerate(("mash", "perm_p01")):
            v = o[(o.brain_ct == bc) & (o.set == s)].set_index("celltype").reindex(cts)
            x = np.arange(len(cts)) + (k - 0.5) * 0.3
            ax.errorbar(x, v.adj_OR, yerr=[v.adj_OR - v.adj_OR_lo, v.adj_OR_hi - v.adj_OR], fmt="o", ms=2.5, lw=0.6,
                        color=sets[s][1], label=sets[s][0])
        ax.axhline(1, ls="--", lw=0.6, color="grey")
        ax.set_yscale("log")
        ax.set_xticks(range(len(cts)), [c.replace("Vascular_endothelial", "Vasc_endo") for c in cts], rotation=45, ha="right")
        ax.set_yticks([0.25, 0.5, 1, 2, 4], ["0.25", "0.5", "1", "2", "4"])
        ax.minorticks_off()
        ax.set(ylabel="adjusted OR", ylim=(0.2, 5), title=f"Overlap with adult {bc} caQTL peaks")
        if j == 0:
            ax.legend(frameon=False, fontsize=5, ncol=2, loc="upper left")
        tag(ax, "ef"[j], x=-0.12)
    save(fig, 8)


# ---------------- S9: full S-LDSC (peak, gene, QTL, link sets) ----------------
FAMS = [("peakspec", "CT-specific peaks"), ("genespec", "CT-specific genes"), ("caQTLpeak", "caQTL peaks"),
        ("eGene", "eGenes"), ("linked", "linked peaks"), ("speclinked", "linked CT-specific peaks")]


def s9():
    d = pd.concat([pd.read_csv(f"{BASE}/GWAS/ldsc/sldsc_{b}.tsv", sep="\t") for b in ("spec", "qtl", "link")])
    d = d[d.primary & d.family.isin([f for f, _ in FAMS])].drop_duplicates(["trait_id", "annot"])
    d["ct"] = d.apply(lambda r: r.annot.replace(f"{r.family}_", ""), axis=1)
    fig = plt.figure(figsize=(7.5, 9.2))
    g = fig.add_gridspec(2, 2, height_ratios=[1, 1.45], hspace=0.42, wspace=0.75)
    ax = fig.add_subplot(g[0, :])
    s = d.groupby(["family", "ct"]).agg(n=("p_tau", "size"), k=("p_tau", lambda x: (x < 0.05).sum()),
                                         f=("FDR", lambda x: (x < 0.1).sum())).reset_index()
    xt, xl, x0 = [], [], 0
    for fam, lab in FAMS:
        x = s[s.family == fam].set_index("ct")
        x = x.reindex([c for c in CTS + ["AllCells"] if c in x.index])
        xs = x0 + np.arange(len(x))
        ax.bar(xs, x.k / x.n, color=[COL.get(c, "#999999") for c in x.index], width=0.8)
        for xi, r in zip(xs, x.itertuples()):
            if r.f:
                ax.text(xi, r.k / r.n + 0.003, str(r.f), ha="center", va="bottom", fontsize=5, fontweight="bold")
        ax.text(x0 + (len(x) - 1) / 2, 0.26, lab, ha="center", fontsize=6)
        xt += list(xs)
        xl += [c.replace("Vascular_endothelial", "Vasc_endo") for c in x.index]
        x0 += len(x) + 1.5
    ax.axhline(0.05, ls="--", color="grey", lw=0.6)
    ax.set_xticks(xt, xl, rotation=90, fontsize=5)
    ax.set(ylabel="fraction of traits with τ p < 0.05", ylim=(0, 0.27),
           title="S-LDSC per annotation (number above bar: traits at FDR < 10%)")
    tag(ax, "a", x=-0.06)
    for k, fams in enumerate((("peakspec", "genespec"), ("caQTLpeak", "eGene"))):
        ax = fig.add_subplot(g[1, k])
        x = d[d.family.isin(fams)].copy()
        x["col"] = x.family.map({"peakspec": "peak", "genespec": "gene", "caQTLpeak": "caQTL", "eGene": "eGene"}) + " " + \
            x.ct.str.replace("Vascular_endothelial", "Vasc_endo")
        best = x.groupby("label").p_tau.min().nsmallest(25).index
        m = x[x.label.isin(best)].pivot_table(index="label", columns="col", values="tau_z", aggfunc="max")
        cols = [f"{p} {c.replace('Vascular_endothelial', 'Vasc_endo')}" for p in
                [{"peakspec": "peak", "genespec": "gene", "caQTLpeak": "caQTL", "eGene": "eGene"}[f] for f in fams]
                for c in CTS + ["AllCells"]]
        m = m[[c for c in cols if c in m.columns]].loc[best]
        sig = x[x.label.isin(best)].pivot_table(index="label", columns="col", values="FDR", aggfunc="min").reindex_like(m)
        sns.heatmap(m, ax=ax, cmap="RdBu_r", center=0, vmin=-4, vmax=4, cbar_kws={"label": "τ z", "shrink": 0.5},
                    annot=sig.applymap(lambda v: "*" if v < 0.1 else ""), fmt="", annot_kws={"fontsize": 6})
        ax.set_yticks(np.arange(len(m)) + 0.5, [t if len(t) <= 34 else t[:32] + "…" for t in m.index], fontsize=5)
        ax.set_xticks(np.arange(m.shape[1]) + 0.5, m.columns, fontsize=4.5, rotation=90)
        ax.set(xlabel="", ylabel="", title=("Specific peaks and genes" if k == 0 else "caQTL peaks and eGenes") +
               "\n(25 traits with lowest p; * FDR < 10%)")
        tag(ax, "bc"[k], x=-0.75, y=1.08)
    save(fig, 9)


# ---------------- S10: colocalization audit ----------------
def s10():
    h = pd.read_csv(f"{BASE}/GWAS/results/coloc_audit/hits_audited.tsv", sep="\t")
    txt = open(f"{BASE}/GWAS/results/coloc_audit/R_coloc_check.txt").read()
    rc = pd.DataFrame([dict(py=float(a), R=float(b)) for a, b in re.findall(r"H4 py=([\d.e-]+) R=([\d.e-]+)", txt)])
    fig, axs = plt.subplots(1, 4, figsize=(7.5, 2.2), gridspec_kw={"wspace": 0.6, "width_ratios": [1, 1, 1, 1.3]})
    dc = {"per-CT": "#7B3294", "pooled": "#E66101"}
    ax = axs[0]
    ax.scatter(rc.py, rc.R, s=4, color="k")
    ax.plot([0.5, 1], [0.5, 1], ls="--", color="grey", lw=0.6)
    dd = np.abs(rc.py - rc.R)
    ax.set(xlabel="PP.H4 (this study)", ylabel="PP.H4 (R coloc 5.2.3)",
           title=f"R check: {(dd < 1e-4).sum()} of {len(rc)} identical\nmax |Δ| = {dd.max():.3f}")
    ax = axs[1]
    for k, x in h.groupby("design"):
        ax.scatter(x["PP.H4"], x["PP.H4_p12_5e-6"], s=4, color=dc[k], label=k)
    ax.axhline(0.5, ls="--", color="grey", lw=0.6)
    ax.plot([0.5, 1], [0.5, 1], ls=":", color="grey", lw=0.6)
    ax.set(xlabel="PP.H4 (p12 = 1e-5)", ylabel="PP.H4 (p12 = 5e-6)",
           title=f"Prior sensitivity\n{(h['PP.H4_p12_5e-6'] > 0.5).sum()} of {len(h)} stay > 0.5")
    ax.legend(frameon=False, fontsize=5)
    ax = axs[2]
    bins = np.linspace(0, 4, 17)
    ax.hist([-np.log10(x.pval_beta) for _, x in h.groupby("design")], bins=bins, stacked=True,
            color=[dc[k] for k, _ in h.groupby("design")])
    ax.axvline(-np.log10(0.05), ls="--", color="grey", lw=0.6)
    ax.set(xlabel="−log10 QTL permutation p", ylabel="colocalizations",
           title=f"QTL strength at hits\n{(h.pval_beta < 0.05).sum()} of {len(h)} with p < 0.05")
    ax = axs[3]
    n = h.groupby(["celltype", "modality"]).size().unstack(fill_value=0)
    n = n.reindex([c for c in ["AllCells"] + CTS if c in n.index])
    n.index = [c.replace("Vascular_endothelial", "Vasc_endo").replace("AllCells", "pooled") for c in n.index]
    n[["ATAC", "RNA"]].plot.barh(stacked=True, ax=ax, color=["#D7191C", "#2C7BB6"], width=0.75)
    ax.invert_yaxis()
    ax.set(xlabel="colocalizations (PP.H4 > 0.5)", ylabel="", title="By QTL map")
    ax.legend(["caQTL", "eQTL"], frameon=False, fontsize=5)
    for k, ax in enumerate(axs):
        tag(ax, "abcd"[k], x=-0.35 if k < 3 else -0.45, y=1.25)
    save(fig, 10)


s4()
s6()
s7()
s8()
s9()
s10()
page(f"{PF}/FigureS_peaks_links.pdf", 1, 1, f"{PF}/FigureS_peaks.png")
page(f"{PF}/FigureS_peaks_links.pdf", 2, 2, f"{PF}/FigureS_links.png")
stack([raster(f"{PF}/fig_peaks_ctspec_epigenome.pdf"), raster(f"{PF}/fig_peaks_ctspec_epigenome15.pdf")], 3)
stack([f"{PF}/panels/locus_TDGF1.png", f"{PF}/panels/locus_GATA6.png"], 5, width=6)
for k, x in enumerate(("cascades", "qc", "ldsc")):
    page(f"{PF}/FigureS_tfnet.pdf", k + 1, 11 + k, f"{PF}/FigureS_tfnet_{x}.png")
shutil.rmtree(TMP)
print("supplementary figures:", sorted(f for f in os.listdir(OUT) if f.endswith(".pdf")))
