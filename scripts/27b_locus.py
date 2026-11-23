"""
Locus view of cell-type-specific scMultiMap links: pseudobulk accessibility (CPM, all donors) of every HDC peak per
cell type (shared scale set by the linked peaks, taller peaks clipped), GENCODE v32 genes, and arcs from the gene TSS to its linked peaks in each cell type (solid = peak-promoter
contact in H1 / H1-derived PCHi-C, dashed = no contact / not captured).
Also importable: draw_locus(gene, fig, subplotspec) (used by 27_figures.py for Fig. 1h).
Usage: [COMPACT=1] python 27b_locus.py GENE [GENE ...]  -> PeakAnnot/figures/panels/locus_<GENE>[_compact].pdf (COMPACT: wide/short for Fig. 1)
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402
from matplotlib.patches import Arc  # noqa: E402

BASE = os.environ["HDC_DIR"] + "/"
PA = f"{BASE}/PeakAnnot/"
CTS = ["iPSC", "Ectoderm_1", "Ectoderm_2", "Ectoderm_3", "Neural_crest", "Neuronal", "Retinal", "Mesoderm",
       "Vascular_endothelial", "Endoderm"]
COL = dict(zip(CTS, sns.color_palette("tab10", len(CTS))))
ES = ["H1", "ME", "MSC", "NPC", "TB"]
plt.rcParams.update({"font.size": 6, "pdf.fonttype": 42})

L = pd.read_csv(f"{PA}/results/links_fdr0.1.tsv.gz", sep="\t")
L = L[(L.covar > 0) & L.celltype.isin(CTS)]
pk = pd.read_parquet(f"{PA}/peak_annotation.parquet", columns=["peak", "chr", "start", "end"] + [f"cpm_{c}" for c in CTS])
genes = pd.read_csv(f"{PA}/ref/genes_gencode32.tsv.gz", sep="\t")
genes = genes[genes.gene_type.isin(["protein_coding", "lncRNA"])]
pc = pd.read_csv(f"{PA}/ref/pchic/pchic_po_hg38.tsv.gz", sep="\t")
pc = pc[pc.cell.isin(ES)]

def draw_locus(gname, fig=None, spec=None, compact=False, title_fs=6.5):
    """draw the locus into fig at subplotspec spec (a new figure when fig is None); returns fig"""
    g = genes[genes.gene_name == gname].iloc[0]
    lk = L[L.gene == gname].copy()
    p = lk.peak.str.split("_", expand=True)
    lk["c"] = (p[1].astype(int) + p[2].astype(int)) // 2
    lo, hi = min(lk.c.min(), g.start) - 25_000, max(lk.c.max(), g.end) + 25_000
    reg = pk[(pk.chr == g.chr) & (pk.end > lo) & (pk.start < hi)]
    cts = [c for c in CTS if c in set(lk.celltype)] + [c for c in CTS if c not in set(lk.celltype)]
    hr = [3] + [1] * len(cts) + [1.2]
    if fig is None:
        fig, axes = plt.subplots(len(cts) + 2, 1, figsize=(7.2, 0.17 * len(cts) + 1.5) if compact else (5.5, 0.32 * len(cts) + 2.2),
                                 sharex=True, gridspec_kw={"height_ratios": hr, "hspace": 0.08})
        own = True
    else:
        sg = spec.subgridspec(len(cts) + 2, 1, height_ratios=hr, hspace=0.08)
        axes = [fig.add_subplot(sg[0])]
        axes += [fig.add_subplot(sg[k], sharex=axes[0]) for k in range(1, len(cts) + 2)]
        for a in axes[:-1]:
            a.tick_params(labelbottom=False)
        own = False
    # arcs
    ax = axes[0]
    contact = pc[(pc.gene == gname) & (pc.chr == g.chr)]
    for r in lk.itertuples():
        ok = ((contact.start <= r.c) & (contact.end >= r.c)).any()
        w = abs(r.c - g.tss)
        ax.add_patch(Arc(((r.c + g.tss) / 2, 0), w, w / (hi - lo) * 2.2, theta1=0, theta2=180, color=COL[r.celltype],
                         lw=0.9, ls="-" if ok else "--"))
    ax.set_ylim(0, 1.15)
    ax.axis("off")
    ax.set_title(f"{gname}: {len(lk)} positive links ({', '.join(f'{c} {n}' for c, n in lk.celltype.value_counts().items())});"
                 " solid = PCHi-C contact", fontsize=title_fs)
    # accessibility tracks
    # shared y scale set by the linked peaks (a strong ubiquitous peak would otherwise flatten all tracks); bars clipped
    vmax = 1.3 * reg[reg.peak.isin(lk.peak)][[f"cpm_{c}" for c in cts]].to_numpy().max()
    for ax, ct in zip(axes[1:-1], cts):
        ax.bar((reg.start + reg.end) / 2, reg[f"cpm_{ct}"].clip(upper=vmax), width=(reg.end - reg.start).clip(lower=(hi - lo) / 400),
               color=COL[ct], linewidth=0)
        linked = reg.peak.isin(lk[lk.celltype == ct].peak)
        ax.scatter(((reg.start + reg.end) / 2)[linked], np.full(linked.sum(), vmax * 0.9), marker="v", s=6, color="k")
        ax.set_ylim(0, vmax)
        ax.set_yticks([])
        ax.set_ylabel(ct.replace("Vascular_endothelial", "Vasc. endothelial").replace("_", " "), rotation=0, ha="right",
                      va="center", fontsize=5.5)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
    # genes
    ax = axes[-1]
    gg = genes[(genes.chr == g.chr) & (genes.end > lo) & (genes.start < hi)].sort_values("start")
    for i, r in enumerate(gg.itertuples()):
        y = i % 3
        ax.plot([max(r.start, lo), min(r.end, hi)], [y, y], color="k" if r.gene_name == gname else "grey",
                lw=2 if r.gene_name == gname else 1)
        ax.text(min(max((r.start + r.end) / 2, lo), hi), y + 0.35, r.gene_name, ha="center", fontsize=4.5,
                fontweight="bold" if r.gene_name == gname else "normal")
    ax.set_ylim(-0.6, 3)
    ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.set_xlim(lo, hi)
    ax.set_xlabel(f"{g.chr} (Mb)")
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v / 1e6:.2f}"))
    if own:
        sfx = "_compact" if compact else ""
        fig.savefig(f"{PA}/figures/panels/locus_{gname}{sfx}.pdf", bbox_inches="tight")
        fig.savefig(f"{PA}/figures/panels/locus_{gname}{sfx}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(gname, "done")
    return axes[0]


if __name__ == "__main__":
    for gname in sys.argv[1:]:
        draw_locus(gname, compact=os.environ.get("COMPACT") == "1")
