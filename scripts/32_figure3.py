"""
Figure 3: TF regulatory networks of the HDC (SPIRE, 29_spire.py / 31_tfnet_analysis.py / 31b_caqtl_motif.py).

  a  multiome pairing: edge Jaccard of paired networks with an independent replicate vs with unpaired networks
  b  per-TF cognate motif enrichment among targets (joint network)
  c  caQTL direction predicted by motif disruption of SPIRE-assigned vs other expressed TFs
  d  key TFs per cell type (joint network): expression z (colour) and target accessibility z (circle)
  e  similarity of within-cell-type networks (Spearman of TF scope; Vascular_endothelial excluded)
  f  regulatory cascades: drivers ordered by expression timing; expression z | target activity z (neuronal, neural crest)
  g  drivers re-used by two lineages: target Jaccard vs expression timing
Supplement (FigureS_tfnet.pdf): all five lineage cascades; QC (sparsity, paired vs unpaired, replicate, motif per network); S-LDSC.
Output: Figures/Main/Figure3.{pdf,png}, FigureS_tfnet.pdf, panels/fig3_*.pdf
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

BASE = os.environ["HDC_DIR"] + "/"
RES = os.environ.get("RESULTS", f"{BASE}/TFNet/results/")
F = os.environ.get("FIGDIR", f"{BASE}/PeakAnnot/figures/")
MAINF = os.environ.get("MAINDIR", f"{BASE}/Figures/Main/")
os.makedirs(MAINF, exist_ok=True)
os.makedirs(f"{F}/panels", exist_ok=True)
CTS = ["iPSC", "Ectoderm_1", "Ectoderm_2", "Ectoderm_3", "Neural_crest", "Neuronal", "Retinal", "Mesoderm", "Endoderm",
       "Vascular_endothelial", "Unknown"]
SHORT = {"Vascular_endothelial": "Vasc_endo"}
LIN = {"neur": "Neuronal", "neur_crest": "Neural crest", "ret": "Retinal", "meso": "Mesoderm", "endo": "Endoderm"}
PAIR = "#2C7BB6"
UNP = "#FDAE61"
plt.rcParams.update({"font.size": 7, "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.titlesize": 7, "axes.titleweight": "bold"})


def lab(c):
    return SHORT.get(c, c)


def jac(a, b, sh):
    return sh / (a + b - sh)


def p_signal(ax):
    """edge overlap of the paired network with an independent replicate (new pseudocells + seed) vs with the unpaired
    network (RNA and ATAC from different cells of the same state): the gap is what single-cell pairing adds"""
    pu = pd.read_csv(f"{RES}/paired_vs_unpaired.tsv", sep="\t").set_index("design")
    rep = pd.read_csv(f"{RES}/replicate_ct.tsv", sep="\t").set_index("celltype")
    designs = [("nn_iPSC", "iPSC", "within\niPSC"), ("nn_Neuronal", "Neuronal", "within\nneuronal"),
               ("traj_neur", "traj_neur", "neuronal\ntrajectory"), ("all", None, "11 cell\ntypes")]
    for i, (d, r, name) in enumerate(designs):
        u = pu.loc[d]
        ju = jac(u.edges_paired, u.edges_unpaired, u.shared)
        if r is not None and r in rep.index:
            ax.bar(i - 0.18, rep.loc[r, "edge_jaccard"], 0.36, color="#555555", label="replicate" if i == 0 else None)
        ax.bar(i + 0.18, ju, 0.36, color=UNP, label="unpaired (same state)" if i == 0 else None)
        ax.text(i + 0.18, ju + 0.01, f"{u.edges_unpaired / u.edges_paired:.0%}", ha="center", fontsize=4.8)
    ax.set_xticks(range(len(designs)), [d[2] for d in designs], fontsize=5.5)
    ax.set(ylabel="edge Jaccard with\npaired network", ylim=(0, 1.0),
           title="Pairing matters most\nwithin a cell state")
    if "ct_iPSC" in pu.index:
        c = pu.loc["ct_iPSC"]
        ax.text(0.5, -0.42, f"random pools within iPSC: {int(c.edges_paired)} vs {int(c.edges_unpaired)} edges; "
                "% = unpaired / paired edges", transform=ax.transAxes, fontsize=4.8, ha="center", va="top")
    ax.legend(frameon=False, fontsize=5.2, loc="upper left", ncol=2)


def p_key(ax, n=3):
    k = pd.read_csv(f"{RES}/key_tfs_all.tsv.gz", sep="\t")
    rows = []
    for ct in CTS:
        t = k[(k.celltype == ct) & (k["rank"] <= n)].sort_values("rank").tf.tolist()
        rows += [x for x in t if x not in rows]
    ez = k.pivot(index="tf", columns="celltype", values="expr_z").reindex(index=rows, columns=CTS)
    az = k.pivot(index="tf", columns="celltype", values="activity_z").reindex(index=rows, columns=CTS)
    im = ax.imshow(ez.to_numpy(), cmap="RdBu_r", vmin=-2.5, vmax=2.5, aspect="auto")
    yy, xx = np.meshgrid(np.arange(len(rows)), np.arange(len(CTS)), indexing="ij")
    a = np.clip(az.to_numpy(), 0, None)
    ax.scatter(xx.ravel(), yy.ravel(), s=a.ravel() * 9, facecolor="none", edgecolor="k", lw=0.5)
    ax.set_xticks(range(len(CTS)), [lab(c) for c in CTS], rotation=60, ha="right", fontsize=5.5)
    ax.set_yticks(range(len(rows)), rows, fontsize=4.6)
    ax.set_title("Key TFs per cell type\n(colour: expression z; circle: target accessibility z)")
    cb = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.ax.tick_params(labelsize=5)
    for s_ in (1, 2):
        ax.scatter([], [], s=s_ * 9, facecolor="none", edgecolor="k", lw=0.5, label=f"{s_}")
    ax.legend(frameon=False, fontsize=5, title="target z", title_fontsize=5, loc="upper left", bbox_to_anchor=(1.12, 0.3))


def p_motif_tf(ax):
    m = pd.read_csv(f"{RES}/motif_per_tf.tsv", sep="\t")
    m = m[m.run == "all"]
    k = pd.read_csv(f"{RES}/key_tfs_all.tsv.gz", sep="\t")
    key = set(k[k["rank"] <= 5].tf)
    sig = m.fdr < 0.05
    ax.scatter(np.log10(m.n), np.log2(m.enrichment), s=4, c=np.where(sig, "#C51B7D", "#BBBBBB"), lw=0)
    ax.axhline(0, color="k", lw=0.5)
    lbl = m[sig & m.tf.isin(key)].nsmallest(9, "p").sort_values("n")
    for i, (_, r) in enumerate(lbl.iterrows()):
        ax.annotate(r.tf, (np.log10(r.n), np.log2(r.enrichment)), xytext=(2, 4 if i % 2 else -6),
                    textcoords="offset points", fontsize=4.8)
    ax.set(xlabel="log10 # targets", ylabel="log2 cognate motif enrichment",
           title=f"Motif support per TF\n{int(sig.sum())}/{len(m)} TFs FDR < 5% (pink)")


def p_ctcorr(ax):
    c = pd.read_csv(f"{RES}/scope_corr_ct.tsv", sep="\t", index_col=0)
    order = [x for x in CTS if x in c.index and x != "Vascular_endothelial"]  # 919 cells: neighbourhoods overlap
    c = c.loc[order, order]
    im = ax.imshow(c.to_numpy(), cmap="viridis", vmin=np.nanmin(c.to_numpy()), vmax=1)
    ax.set_xticks(range(len(order)), [lab(x) for x in order], rotation=60, ha="right", fontsize=5.5)
    ax.set_yticks(range(len(order)), [lab(x) for x in order], fontsize=5.5)
    t = "Within-cell-type networks\n(Spearman of TF scope)"
    if os.path.exists(f"{RES}/replicate_ct.tsv"):
        r = pd.read_csv(f"{RES}/replicate_ct.tsv", sep="\t")
        r = r[r.celltype.isin(["iPSC", "Neuronal"])]
        t += "\nreplicate: " + ", ".join(f"{a} {b:.2f}" for a, b in zip(r.celltype, r.scope_spearman))
    ax.set_title(t)
    cb = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
    cb.ax.tick_params(labelsize=5)


def cascade(ax, lin, n=14, title=True):
    t = pd.read_csv(f"{RES}/traj_tfs.tsv.gz", sep="\t")
    t = t[(t.lineage == lin) & t.driver].nsmallest(n, "driver_rank").sort_values("com_expr")
    cv = np.load(f"{RES}/traj_{lin}_curves.npz", allow_pickle=True)
    ti = pd.Series(np.arange(len(cv["tf"])), index=cv["tf"])
    j = ti.reindex(t.tf).to_numpy()
    M = np.hstack([cv["expr_z"][:, j].T, np.full((len(j), 1), np.nan), cv["act_z"][:, j].T])
    ax.imshow(M, cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
    nb = cv["expr_z"].shape[0]
    ax.set_yticks(range(len(j)), t.tf, fontsize=5)
    ax.set_xticks([nb / 2 - 0.5, nb + nb / 2 + 0.5], ["TF expression", "target accessibility"], fontsize=5.5)
    ax.tick_params(axis="x", length=0)
    for x0 in (0, nb + 1):
        ax.annotate("", xy=(x0 + nb - 0.6, len(j) - 0.2), xytext=(x0 - 0.4, len(j) - 0.2),
                    arrowprops=dict(arrowstyle="->", lw=0.6), annotation_clip=False)
    if title:
        ax.set_title(f"{LIN[lin]} lineage\n(iPSC → {LIN[lin].lower()}, pseudotime →)")


def p_shared(ax):
    """drivers re-used by two lineages: target overlap vs the driver's expression timing"""
    o = pd.read_csv(f"{RES}/traj_shared_driver_targets.tsv", sep="\t")
    t = pd.read_csv(f"{RES}/traj_tfs.tsv.gz", sep="\t").set_index(["lineage", "tf"]).com_expr
    nb = 10
    o["timing"] = [(t[(a, f)] + t[(b, f)]) / 2 / (nb - 1) for f, a, b in zip(o.tf, o.lin_a, o.lin_b)]
    y = np.log10(o.jaccard + 1e-3)
    early = o.timing < 0.35
    ax.scatter(o.timing, y, s=7, c=np.where(early, "#E08214", "#542788"), lw=0)
    g = o[o.jaccard > 0].assign(y=y).groupby("tf").agg(x=("timing", "mean"), y=("y", "max"))
    for f, r in g.iterrows():
        ax.text(r.x, r.y + 0.08, f, fontsize=4.6, ha="center")
    z = sorted(o[o.jaccard == 0].tf.unique())
    ax.text(0.99, 0.84, "no shared targets:\n" + "\n".join(", ".join(z[i:i + 3]) for i in range(0, len(z), 3)),
            transform=ax.transAxes, fontsize=4.4, ha="right", va="top", color="#542788")
    ax.set_yticks([-3, -2, -1, 0], ["0", "0.01", "0.1", "1"])
    ax.set(xlabel="driver expression timing\n(pseudotime, 0-1)", ylabel="target Jaccard\nbetween lineages",
           title=f"Re-used drivers rewire targets\n{(~early).sum()} late pairs: median J = {o.jaccard[~early].median():.3f}")


def p_ldsc(ax):
    f = f"{BASE}/GWAS/ldsc/sldsc_tf.tsv"
    if not os.path.exists(f):
        ax.text(0.5, 0.5, "S-LDSC pending", ha="center", transform=ax.transAxes)
        return
    d = pd.read_csv(f, sep="\t")
    d = d[d.primary]
    hit = d[d.FDR < 0.1]
    traits = hit.label.unique().tolist()
    if not traits:
        traits = d.sort_values("p_tau").label.unique()[:12].tolist()
    d = d[d.label.isin(traits)]
    M = d.pivot_table(index="label", columns="annot", values="tau_z")
    cols = [f"tfprog_{c}" for c in CTS if f"tfprog_{c}" in M.columns] + [f"tfdyn_{l}" for l in LIN if f"tfdyn_{l}" in M.columns]
    M = M.reindex(columns=cols)
    im = ax.imshow(M.to_numpy(), cmap="RdBu_r", vmin=-4, vmax=4, aspect="auto")
    fd = d.pivot_table(index="label", columns="annot", values="FDR").reindex(index=M.index, columns=cols)
    yy, xx = np.where(fd.to_numpy() < 0.1)
    ax.scatter(xx, yy, s=6, c="k")
    ax.set_xticks(range(len(cols)), [lab(c.split("_", 1)[1]) if c.startswith("tfprog") else LIN[c[6:]] + " (dyn)"
                                     for c in cols], rotation=60, ha="right", fontsize=5.5)
    ax.set_yticks(range(len(M)), [x[:32] for x in M.index], fontsize=5)
    ax.set_title("TF programs and trait heritability\n(S-LDSC tau z beyond CT-specific peaks)")
    cb = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.ax.tick_params(labelsize=5)


def p_caqtl(ax):
    if not os.path.exists(f"{RES}/caqtl_motif_concordance.tsv"):
        ax.text(0.5, 0.5, "caQTL motif test pending", ha="center", transform=ax.transAxes)
        return
    c = pd.read_csv(f"{RES}/caqtl_motif_concordance.tsv", sep="\t")
    c = c[(c.network == "either") & (c.level == "peak")]
    x = 0
    ticks = []
    for cq, g in c.groupby("caqtl", sort=False):
        for a, col in ((True, "#1B7837"), (False, "#BBBBBB")):
            r = g[g.assigned == a].iloc[0]
            ax.bar(x, r.frac, 0.8, color=col, label=("SPIRE-assigned TF" if a else "other expressed TF") if cq == c.caqtl.iloc[0] else None)
            ax.text(x, r.frac + 0.01, f"n={int(r.n)}", ha="center", fontsize=5)
            x += 1
        ax.text(x - 1.5, g.frac.max() + 0.1, f"p={g.p_fisher.iloc[0]:.2g}", ha="center", fontsize=5.5)
        ticks.append((x - 1.5, cq.replace(", ", "\n")))
        x += 0.6
    ax.axhline(0.5, ls="--", lw=0.5, color="grey")
    ax.set_xticks([t[0] for t in ticks], [t[1] for t in ticks], fontsize=5.5)
    ax.set(ylabel="caQTL direction predicted\nby motif change", ylim=(0.3, 1.12),
           title="Genetic validation of\nTF→peak edges")
    ax.legend(frameon=False, fontsize=5.2, loc="upper left", ncol=1)


PANELS = [p_signal, p_key, p_motif_tf, p_ctcorr, p_shared, p_ldsc, p_caqtl]
if __name__ == "__main__":
    fig = plt.figure(figsize=(7.2, 8.2))
    gs = fig.add_gridspec(3, 6, hspace=0.9, wspace=1.9, height_ratios=[1, 1.75, 1.2])
    axes = {"a": fig.add_subplot(gs[0, 0:2]), "b": fig.add_subplot(gs[0, 2:4]), "c": fig.add_subplot(gs[0, 4:6]),
            "d": fig.add_subplot(gs[1, 0:3]), "e": fig.add_subplot(gs[1, 3:6]),
            "f1": fig.add_subplot(gs[2, 0:2]), "f2": fig.add_subplot(gs[2, 2:4]), "g": fig.add_subplot(gs[2, 4:6])}
    p_signal(axes["a"])
    p_motif_tf(axes["b"])
    p_caqtl(axes["c"])
    p_key(axes["d"])
    p_ctcorr(axes["e"])
    cascade(axes["f1"], "neur")
    cascade(axes["f2"], "neur_crest")
    p_shared(axes["g"])
    for k, ax in axes.items():
        if k == "f2":
            continue
        ax.text(-0.22, 1.16, k[0], transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")
    fig.savefig(f"{MAINF}/Figure3.pdf", bbox_inches="tight")
    fig.savefig(f"{MAINF}/Figure3.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    for fn in PANELS:
        f, ax = plt.subplots(figsize=(3.4, 3))
        fn(ax)
        f.savefig(f"{F}/panels/fig3_{fn.__name__[2:]}.pdf", bbox_inches="tight")
        plt.close(f)

    # supplement
    with PdfPages(f"{F}/FigureS_tfnet.pdf") as pdf:
        fig, axs = plt.subplots(1, 5, figsize=(11, 3.6))
        for ax, lin in zip(axs, LIN):
            if os.path.exists(f"{RES}/traj_{lin}_curves.npz"):
                cascade(ax, lin, n=20)
        fig.suptitle("SPIRE drivers along each lineage (top 20 by # targets, ordered by expression timing)", fontsize=8)
        fig.tight_layout()
        for k, a in enumerate(axs):
            a.text(-0.25, 1.12, "abcde"[k], transform=a.transAxes, fontsize=10, fontweight="bold", va="bottom")
        pdf.savefig(fig, bbox_inches="tight")
        fig.savefig(f"{F}/FigureS_tfnet_cascades.png", dpi=200, bbox_inches="tight")
        plt.close(fig)
        fig, axs = plt.subplots(1, 4, figsize=(12, 3.2))
        mg = pd.read_csv(f"{RES}/motif_global.tsv", sep="\t")
        mg = mg[(mg.sign == "pos") & ~mg.run.str.contains("unpaired|_s2")]
        yy = np.arange(len(mg))
        axs[3].barh(yy, mg.enrichment, color=PAIR, label="observed")
        axs[3].errorbar(mg.null_mean, yy, xerr=2 * mg.null_sd, fmt="|", color="k", ms=4, label="permuted TF labels")
        axs[3].set_yticks(yy, mg.run, fontsize=5)
        axs[3].axvline(1, color="k", lw=0.5)
        axs[3].set(xlabel="cognate motif enrichment (positive edges)", title="Motif validation per network")
        axs[3].legend(frameon=False, fontsize=5)
        s = pd.read_csv(f"{RES}/run_summary.tsv", sep="\t")
        s = s[~s.run.str.contains("unpaired|_s2")]
        axs[0].barh(s.run, s.tfs_per_peak, color=PAIR)
        axs[0].set(xlabel="stable TFs per tested peak", title="Network sparsity")
        axs[0].tick_params(axis="y", labelsize=5)
        if os.path.exists(f"{RES}/paired_vs_unpaired.tsv"):
            pu = pd.read_csv(f"{RES}/paired_vs_unpaired.tsv", sep="\t")
            y = np.arange(len(pu))
            axs[1].barh(y + 0.2, pu.edges_paired, 0.4, color=PAIR, label="paired")
            axs[1].barh(y - 0.2, pu.edges_unpaired, 0.4, color=UNP, label="unpaired")
            axs[1].set_yticks(y, pu.design, fontsize=5.5)
            axs[1].set(xlabel="consensus edges", title="Paired vs unpaired pseudocells")
            axs[1].legend(frameon=False, fontsize=5.5)
        if os.path.exists(f"{RES}/replicate_ct.tsv"):
            r = pd.read_csv(f"{RES}/replicate_ct.tsv", sep="\t")
            axs[2].bar(r.celltype, r.scope_spearman, color="#555555")
            axs[2].set(ylabel="Spearman (scope), independent replicate", ylim=(0, 1), title="Reproducibility")
        fig.tight_layout()
        for k, a in enumerate(axs):
            a.text(-0.3, 1.03, "abcde"[k], transform=a.transAxes, fontsize=10, fontweight="bold", va="bottom")
        pdf.savefig(fig, bbox_inches="tight")
        fig.savefig(f"{F}/FigureS_tfnet_qc.png", dpi=200, bbox_inches="tight")
        plt.close(fig)
        fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.5, 3.8), gridspec_kw={"width_ratios": [2.2, 1], "wspace": 0.35})
        p_ldsc(ax)
        sm = f"{RES}/sldsc_tf_summary.tsv"
        if os.path.exists(sm):
            t = pd.read_csv(sm, sep="\t")
            want = [f"tfprog_{c}" for c in CTS] + [f"tfdyn_{l}" for l in LIN]
            t = t.set_index("set").reindex([w for w in want if w in set(t.set)]).reset_index()
            nm = [lab(x.split("_", 1)[1]) if x.startswith("tfprog") else LIN[x[6:]] + " (dyn)" for x in t.set]
            yy = np.arange(len(t))[::-1]
            ax2.barh(yy, t.n_p05, color=np.where(t.p_excess < 0.05, "#B2182B", "#888888"))
            ax2.axvline(0.05 * t.traits.iloc[0], ls="--", color="k", lw=0.6)
            ax2.set_yticks(yy, nm, fontsize=5.5)
            ax2.set(xlabel=f"traits with τ p < 0.05 (of {t.traits.iloc[0]})", title="Nominal results per set\n(dashed = 5% expected; red = binomial p < 0.05)")
        ax.set_title("S-LDSC τ z of TF programs (top 12 traits by p)\nno result at FDR < 10%")
        for k, a in enumerate((ax, ax2)):
            a.text(-0.2, 1.06, "abcde"[k], transform=a.transAxes, fontsize=10, fontweight="bold", va="bottom")
        pdf.savefig(fig, bbox_inches="tight")
        fig.savefig(f"{F}/FigureS_tfnet_ldsc.png", dpi=200, bbox_inches="tight")
        plt.close(fig)
    print("figures written")
