"""
Do SPIRE TF -> peak edges predict the direction of caQTL effects?

caQTL SNPs: per cell type, every SNV inside the peak (tensorQTL nominal pairs, caQTL_Output_optCov/<CT>) with
  nominal p < 0.05 in a mashR caQTL peak (lfsr < 0.05 in that cell type), or p < 1e-3 in any peak (with 17 donors the
  permutation lead is rarely inside the peak); alleles from GWAS/variant_universe.parquet (tensorQTL slope = per ALT
  allele; REF checked against hg38).
For every HOMER known motif: best log-odds over windows covering the SNP (both strands) for REF and ALT sequence;
motif "hit" if max(REF, ALT) >= HOMER threshold; delta = ALT - REF. Per TF (motifs mapped in 30b): the hit motif with
the largest |delta|.
Peak x TF pairs with a hit and |delta| >= 1 are split into
  SPIRE-assigned   consensus edge TF -> peak in that cell type's network (nn_<CT>) or the joint network (all);
                   predicted sign = sign(delta) x sign(coef)
  not assigned     TF tested in the network (expressed) but no edge; predicted sign = sign(delta) (activator assumption)
Concordance = predicted sign == sign(caQTL slope). Peak-level: one pair per peak (largest |delta|).
Also: the TFs assigned to colocalized caQTL peaks (GWAS/results/coloc_audit/hits_audited.tsv) with motif delta at the
best coloc SNP.
Output TFNet/results/caqtl_motif_pairs.tsv.gz, caqtl_motif_concordance.tsv, coloc_peak_tfs.tsv
"""
import os

import numpy as np
import pandas as pd
from numpy.lib.stride_tricks import sliding_window_view
from scipy.stats import binomtest, fisher_exact

BASE = os.environ["HDC_DIR"] + "/"
RUNS = os.environ.get("RUNS", f"{BASE}/TFNet/runs/")
RES = os.environ.get("RESULTS", f"{BASE}/TFNet/results/")
FA = f"{BASE}/PeakAnnot/ref/hg38.fa"
MOT = os.environ["HOMER_MOTIFS"]
TAG = "noBootstrap_TMM_invnorm"
FL = 30
CTS = ["iPSC", "Ectoderm_1", "Ectoderm_3", "Neural_crest", "Neuronal", "Retinal", "Mesoderm", "Endoderm",
       "Vascular_endothelial", "Unknown"]
log = open(f"{RES}/caqtl_motif_summary.txt", "w")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    log.write(s + "\n")


# ---- fasta access through .fai
fai = pd.read_csv(FA + ".fai", sep="\t", header=None, index_col=0, names=["len", "off", "lb", "lB"])
fh = open(FA, "rb")


def fetch(chrom, start, end):  # 0-based half-open
    ln, off, lb, lB = fai.loc[chrom]
    s = off + (start // lb) * lB + start % lb
    e = off + (end // lb) * lB + end % lb
    fh.seek(s)
    return fh.read(e - s).replace(b"\n", b"").decode().upper()


# ---- motifs
motifs = []
for blk in open(MOT).read().split(">")[1:]:
    lines = blk.strip().split("\n")
    h = lines[0].split("\t")
    P = np.array([[float(v) for v in l.split()] for l in lines[1:]]) + 1e-3
    P /= P.sum(1, keepdims=True)
    motifs.append((h[1], float(h[2]), np.log(P / 0.25)))
mt = pd.read_csv(f"{BASE}/TFNet/motif/motifs.tsv", sep="\t").fillna("")
sym = dict(zip(mt.motif, mt.symbols))
say("motifs", len(motifs))

ENC = {"A": 0, "C": 1, "G": 2, "T": 3}


def onehot(seqs):
    X = np.zeros((len(seqs), 2 * FL + 1, 4), np.float32)
    for i, s in enumerate(seqs):
        for j, c in enumerate(s):
            if c in ENC:
                X[i, j, ENC[c]] = 1
    return X


def scores(X):
    """(n SNPs x motifs) best log-odds over windows covering the centre, both strands."""
    out = np.full((X.shape[0], len(motifs)), -np.inf, np.float32)
    for k, (_, _, M) in enumerate(motifs):
        w = M.shape[0]
        Mr = M[::-1, ::-1]  # reverse complement (A C G T -> T G C A reversed)
        lo, hi = max(FL - w + 1, 0), FL  # window starts covering position FL
        V = sliding_window_view(X[:, lo:hi + w], (w, 4), axis=(1, 2))[:, :, 0]  # n x nwin x w x 4
        f = np.einsum("nkwb,wb->nk", V, M)
        r = np.einsum("nkwb,wb->nk", V, Mr)
        out[:, k] = np.maximum(f.max(1), r.max(1))
    return out


# ---- caQTL SNPs inside their peak (nominal pairs; with 17 donors the permutation lead is rarely in the peak)
import pyarrow.dataset as ds  # noqa: E402

u = pd.read_parquet(f"{BASE}/GWAS/variant_universe.parquet", columns=["id", "chr", "pos38", "ref", "alt"])
u = u.drop_duplicates("id").set_index("id")
mash = pd.read_csv(f"{BASE}/mash/ATAC_mash_strong.tsv.gz", sep="\t", index_col=0)
leads = []
for ct in CTS:
    sig = set(mash.index[mash[f"lfsr_{ct}"] < 0.05]) if f"lfsr_{ct}" in mash else set()
    for c in range(1, 23):
        f = f"{BASE}/caQTL_Output_optCov/{ct}/chr{c}_cis_nominal_{TAG}.cis_qtl_pairs.{c}.parquet"
        t = ds.dataset(f).to_table(columns=["phenotype_id", "variant_id", "start_distance", "end_distance", "pval_nominal",
                                            "slope"],
                                   filter=(ds.field("pval_nominal") < 0.05) & (ds.field("start_distance") >= 0)
                                   & (ds.field("end_distance") <= 0)).to_pandas()
        t["mash"] = t.phenotype_id.isin(sig)
        t = t[t.mash | (t.pval_nominal < 1e-3)]
        leads.append(t.assign(celltype=ct))
    print(ct, sum(len(x) for x in leads if len(x) and x.celltype.iloc[0] == ct), flush=True)
L = pd.concat(leads, ignore_index=True).join(u, on="variant_id").dropna(subset=["pos38"])
L = L[(L.ref.str.len() == 1) & (L.alt.str.len() == 1)]
L["pchr"] = L.phenotype_id.str.split("_").str[0]
L["pos38"] = L.pos38.astype(int)
seqs = [fetch(c, p - 1 - FL, p + FL) for c, p in zip(L.pchr, L.pos38)]
ok = np.array([s[FL] == r for s, r in zip(seqs, L.ref)])
say("in-peak caQTL SNVs (mashR peak & p < 0.05, or p < 1e-3):", len(L), "| mashR:", int(L.mash.sum()),
    "| REF matches hg38:", round(ok.mean(), 4))
L, seqs = L[ok].reset_index(drop=True), [s for s, o in zip(seqs, ok) if o]
alt = [s[:FL] + a + s[FL + 1:] for s, a in zip(seqs, L.alt)]
Sr, Sa = scores(onehot(seqs)), scores(onehot(alt))
thr = np.array([m[1] for m in motifs], np.float32)
hit = np.maximum(Sr, Sa) >= thr
D = np.where(hit, Sa - Sr, 0)
say("SNPs with >= 1 motif hit changed by |delta| >= 1:", int((np.abs(D) >= 1).any(1).sum()))

# per TF: motif with largest |delta|
names = [m[0] for m in motifs]
tf_of = {}
for k, n in enumerate(names):
    for t in sym.get(n, "").split(","):
        if t:
            tf_of.setdefault(t, []).append(k)
rows = []
for t, ks in tf_of.items():
    sub = D[:, ks]
    j = np.abs(sub).argmax(1)
    d = sub[np.arange(len(L)), j]
    m = np.abs(d) >= 1
    if m.any():
        rows.append(pd.DataFrame({"i": np.where(m)[0], "tf": t, "delta": d[m]}))
pairs = pd.concat(rows).join(L[["celltype", "phenotype_id", "variant_id", "slope", "pval_nominal", "mash"]], on="i")


# ---- SPIRE edges
def edges(run):
    p = f"{RUNS}/{run}/edges.tsv.gz"
    if not os.path.exists(p):
        return None, set()
    e = pd.read_csv(p, sep="\t")
    tested = set(np.load(f"{RUNS}/{run}/strata.npz", allow_pickle=True)["tf"])
    return e.groupby(["tf", "peak"]).coef.first(), tested


ea, ta = edges("all")
out = []
for ct, g in pairs.groupby("celltype"):
    ec, tc = edges(f"nn_{ct}")
    for net, e, tested in (("ct", ec, tc), ("all", ea, ta)):
        if e is None:
            continue
        g2 = g[g.tf.isin(tested)].copy()
        key = list(zip(g2.tf, g2.phenotype_id))
        coef = e.reindex(pd.MultiIndex.from_tuples(key)).to_numpy() if key else np.array([])
        g2["coef"] = coef
        g2["network"] = net
        out.append(g2)
pairs = pd.concat(out, ignore_index=True)
pairs["assigned"] = pairs.coef.notna()
pairs["pred"] = np.sign(pairs.delta) * np.where(pairs.assigned, np.sign(pairs.coef.fillna(1)), 1)
pairs["concordant"] = pairs.pred == np.sign(pairs.slope)
pairs.to_csv(f"{RES}/caqtl_motif_pairs.tsv.gz", sep="\t", index=False, float_format="%.4g")

res = []
for net in ("ct", "all", "either"):
    d = pairs if net == "either" else pairs[pairs.network == net]
    if net == "either":
        d = d.sort_values("assigned", ascending=False).drop_duplicates(["celltype", "phenotype_id", "tf"])
    for strong in (False, True):
        dd = d[(d.pval_nominal < 1e-3)] if strong else d[d.mash]
        for lvl in ("pair", "peak"):
            for a in (True, False):
                x = dd[dd.assigned == a]
                if lvl == "peak":
                    x = x.loc[x.delta.abs().groupby([x.celltype, x.phenotype_id]).idxmax()] if len(x) else x
                n, k = len(x), int(x.concordant.sum())
                res.append(dict(network=net, caqtl="all peaks, p<1e-3" if strong else "mashR peaks, p<0.05", level=lvl,
                                assigned=a, n=n, concordant=k, frac=k / n if n else np.nan,
                                p_binom=binomtest(k, n, 0.5, alternative="greater").pvalue if n else np.nan))
res = pd.DataFrame(res)
for (net, cq, lvl), g in res.groupby(["network", "caqtl", "level"]):
    a, b = g[g.assigned].iloc[0], g[~g.assigned].iloc[0]
    if a.n and b.n:
        o, p = fisher_exact([[a.concordant, a.n - a.concordant], [b.concordant, b.n - b.concordant]], alternative="greater")
        res.loc[g.index, "OR_assigned_vs_not"] = o
        res.loc[g.index, "p_fisher"] = p
res.to_csv(f"{RES}/caqtl_motif_concordance.tsv", sep="\t", index=False, float_format="%.4g")
say(res.to_string(index=False))

# ---- colocalized caQTL peaks: assigned TFs and motif change at the best coloc SNP
h = pd.read_csv(f"{BASE}/GWAS/results/coloc_audit/hits_audited.tsv", sep="\t")
h = h[h.modality == "ATAC"].sort_values("PP.H4", ascending=False).drop_duplicates(["celltype", "phenotype_id"])
rows = []
for _, r in h.iterrows():
    tfs = []
    for net, e in (("all", ea), (f"nn_{r.celltype}", edges(f"nn_{r.celltype}")[0])):
        if e is None:
            continue
        try:
            s = e.xs(r.phenotype_id, level=1)
            tfs += [f"{t}{'+' if c > 0 else '-'}[{net}]" for t, c in s.items()]
        except KeyError:
            pass
    chrom = "chr" + str(r.chr)
    pos = int(r.best_pos38)
    vu = u[(u.chr == r.chr) & (u.pos38 == pos)]
    md = ""
    if len(vu) == 1 and len(vu.ref.iloc[0]) == 1 and len(vu.alt.iloc[0]) == 1:
        s = fetch(chrom, pos - 1 - FL, pos + FL)
        a = s[:FL] + vu.alt.iloc[0] + s[FL + 1:]
        sr, sa = scores(onehot([s])), scores(onehot([a]))
        dd = np.where(np.maximum(sr, sa) >= thr, sa - sr, 0)[0]
        top = np.argsort(-np.abs(dd))[:3]
        md = "; ".join(f"{names[k].split('/')[0]} {dd[k]:+.1f}" for k in top if abs(dd[k]) >= 1)
    rows.append(dict(celltype=r.celltype, peak=r.phenotype_id, trait=r.label, PP_H4=r["PP.H4"], best_snp=r.best_snp,
                     snp_in_peak=r.best_pos38 >= int(r.phenotype_id.split("_")[1]) and r.best_pos38 <= int(r.phenotype_id.split("_")[2]),
                     spire_tfs=", ".join(tfs), motif_change_best_snp=md))
cp = pd.DataFrame(rows)
cp.to_csv(f"{RES}/coloc_peak_tfs.tsv", sep="\t", index=False, float_format="%.3g")
say(cp.to_string(index=False))
