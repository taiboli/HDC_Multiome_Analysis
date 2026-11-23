"""
TF regulatory networks from the SPIRE runs (29_spire.py): validation, key TFs per cell type, cross-cell-type
concordance / divergence, dynamics along differentiation, paired vs unpaired, S-LDSC sets.

  1. network summary per run (edges, peaks with a stable TF, TFs per peak, positive fraction)
  2. motif concordance (HOMER scan, 30b): fraction of consensus edges whose peak carries the TF's cognate motif vs the
     TF's motif frequency among the run's tested peaks; null = TF labels permuted across edges (100x; keeps every TF's
     and every peak's number of edges); per TF: binomial test, BH
  3. key TFs per cell type (joint run "all"): expression z across the 11 cell types x target activity
     (signed mean z of the TF's consensus targets' accessibility across cell types); TF expressed above the cell type's
     median TF, >= 30 targets
  4. within-cell-type networks (nn_<CT>, neighbourhood pseudocells): scope (fraction of tested peaks) per TF and cell type; Spearman correlation of
     scope between cell types on TFs tested in both; edge Jaccard on shared TFs x peaks; replicate ceiling (_s2);
     hub TFs (top scope) and cell-type-biased hubs (median-centred log2 scope vs the other cell types)
  5. trajectories (traj_<lin>): per TF expression and target activity over 10 pseudotime bins; dynamic drivers
     (expression range >= 1 log2 unit, >= 50 targets, expression-activity r >= 0.5); timing = centre of mass over bins;
     drivers shared between lineages and their target overlap (Jaccard / random expectation)
  6. paired vs unpaired pseudocells (same strata, RNA and ATAC from different cells)
  7. S-LDSC sets (batch "tf"): tfprog_<CT> = positive targets of the cell type's top 10 key TFs accessible in it
     (detected in >= 2% of its cells); tfdyn_<lin> = positive targets of the lineage's top 10 drivers
Output TFNet/results/
Usage: python 31_tfnet_analysis.py   (RUNS / RESULTS env vars override the input / output folders)
"""
import glob
import os

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import binom, pearsonr, spearmanr

BASE = os.environ["HDC_DIR"] + "/"
RUNS = os.environ.get("RUNS", f"{BASE}/TFNet/runs/")
RES = os.environ.get("RESULTS", f"{BASE}/TFNet/results/")
os.makedirs(RES, exist_ok=True)
CTS = ["iPSC", "Ectoderm_1", "Ectoderm_2", "Ectoderm_3", "Neural_crest", "Neuronal", "Retinal", "Mesoderm", "Endoderm",
       "Vascular_endothelial", "Unknown"]
LINS = ["neur", "neur_crest", "ret", "meso", "endo"]
W = "nn_"  # within-cell-type design (29_spire.py); random pooling (ct_) is the negative control
rng = np.random.default_rng(1)
log = open(f"{RES}/summary.txt", "w")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    log.write(s + "\n")


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty(len(p))
    out[o] = np.minimum(q, 1)
    return out


def pick(M, r, c):
    """M[r, c] for index arrays as a dense 1-d array (empty-safe)"""
    if len(r) == 0:
        return np.zeros(0)
    return np.asarray(M[r, c]).ravel()


def load(run):
    d = f"{RUNS}/{run}/"
    if not os.path.exists(f"{d}/edges.tsv.gz"):
        return None
    st = dict(np.load(f"{d}/strata.npz", allow_pickle=True))
    e = pd.read_csv(f"{d}/edges.tsv.gz", sep="\t")
    sc = pd.read_csv(f"{d}/scope.tsv", sep="\t")
    return {"st": st, "e": e, "sc": sc, "run": run}


runs = sorted(os.path.basename(p.rstrip("/")) for p in glob.glob(f"{RUNS}/*/") if os.path.exists(f"{p}/edges.tsv.gz"))
R = {r: load(r) for r in runs}
say("runs:", ", ".join(runs))

# ---------- 1. summary
rows = []
for r, x in R.items():
    e, st = x["e"], x["st"]
    P, T = len(st["peak"]), len(st["tf"])
    F = sp.load_npz(f"{RUNS}/{r}/freq.npz")
    rows.append(dict(run=r, strata=len(st["strata"]), cells=int(st["n_cells"].sum()), tfs_tested=T, peaks_tested=P,
                     edges_freq_ge_01=F.nnz, edges=len(e), peaks_with_tf=e.peak.nunique(),
                     frac_peaks_with_tf=e.peak.nunique() / P, tfs_with_target=e.tf.nunique(), tfs_per_peak=len(e) / P,
                     tfs_per_targeted_peak=len(e) / max(e.peak.nunique(), 1), frac_positive=(e.coef > 0).mean()))
summ = pd.DataFrame(rows)
summ.to_csv(f"{RES}/run_summary.tsv", sep="\t", index=False, float_format="%.4g")
say(summ.to_string(index=False))

# ---------- 2. motif concordance
fa = pd.read_csv(f"{BASE}/PeakAnnot/sc/features_ATAC.tsv.gz", sep="\t", usecols=["id"]).id
pidx = pd.Series(np.arange(len(fa)), index=fa)
PT = sp.load_npz(f"{BASE}/TFNet/motif/peak_tf_motif.npz").tocsr()
tfm = pd.read_csv(f"{BASE}/TFNet/motif/tf_motif.tsv", sep="\t")
tcol = pd.Series(np.arange(len(tfm)), index=tfm.tf)
glob_rows, tf_rows = [], []
for r, x in R.items():
    e, st = x["e"], x["st"]
    rows_t = pidx.reindex(st["peak"]).to_numpy()
    frac = np.asarray(PT[rows_t].mean(0)).ravel()  # motif frequency among tested peaks, per motif TF
    em = e[e.tf.isin(tcol.index)].copy()
    pr, tc = pidx.reindex(em.peak).to_numpy(), tcol.reindex(em.tf).to_numpy()
    em["hit"] = pick(PT, pr, tc).astype(bool)
    em["exp"] = frac[tc]
    for sign, d in (("all", em), ("pos", em[em.coef > 0]), ("neg", em[em.coef < 0])):
        if len(d) == 0:
            continue
        obs = d.hit.sum() / d.exp.sum()
        pr_d, tc_d = pidx.reindex(d.peak).to_numpy(), tcol.reindex(d.tf).to_numpy()
        null = []
        for _ in range(100):
            tp = rng.permutation(tc_d)
            null.append(pick(PT, pr_d, tp).sum() / frac[tp].sum())
        null = np.array(null)
        glob_rows.append(dict(run=r, sign=sign, edges_with_motif_tf=len(d), obs_rate=d.hit.mean(), exp_rate=d.exp.mean(),
                              enrichment=obs, null_mean=null.mean(), null_sd=null.std(),
                              z=(obs - null.mean()) / null.std() if null.std() > 0 else np.nan))
    g = em[em.coef > 0].groupby("tf").agg(n=("hit", "size"), k=("hit", "sum"))
    g["p_motif"] = frac[tcol.reindex(g.index).to_numpy()]
    g = g[g.n >= 50]
    if len(g):
        g["enrichment"] = (g.k / g.n) / g.p_motif
        g["p"] = binom.sf(g.k - 1, g.n, g.p_motif)
        g["fdr"] = bh(g.p.to_numpy())
        tf_rows.append(g.reset_index().assign(run=r))
mg = pd.DataFrame(glob_rows)
mg.to_csv(f"{RES}/motif_global.tsv", sep="\t", index=False, float_format="%.4g")
mt = pd.concat(tf_rows)
mt.to_csv(f"{RES}/motif_per_tf.tsv", sep="\t", index=False, float_format="%.4g")
say("\nmotif concordance (positive edges):\n" + mg[mg.sign == "pos"].to_string(index=False))
say("per-TF motif enrichment FDR<0.05 (pos edges, >= 50 targets):\n" + mt.groupby("run").apply(
    lambda d: f"{(d.fdr < 0.05).sum()}/{len(d)}").to_string())


# ---------- helpers: z-scored stratum means and target activity
def activity(x, min_targets=1):
    st, e = x["st"], x["e"]
    z = lambda M: (M - M.mean(0)) / np.where(M.std(0) > 0, M.std(0), 1)  # noqa: E731
    ztf, zpk = z(st["mean_tf"].astype(float)), z(st["mean_peak"].astype(float))
    tfs, peaks = st["tf"], st["peak"]
    ti = pd.Series(np.arange(len(tfs)), index=tfs).reindex(e.tf).to_numpy()
    pi = pd.Series(np.arange(len(peaks)), index=peaks).reindex(e.peak).to_numpy()
    n = np.bincount(ti, minlength=len(tfs))
    W = sp.csr_matrix((np.sign(e.coef.to_numpy()) / np.maximum(n[ti], 1), (pi, ti)), shape=(len(peaks), len(tfs)))
    act = np.asarray(zpk @ W) if sp.issparse(zpk) else (W.T @ zpk.T).T
    npos = np.bincount(ti[e.coef.to_numpy() > 0], minlength=len(tfs))
    return ztf, np.asarray(act), n, npos


CANON = {"iPSC": ["POU5F1", "NANOG", "SOX2", "ZSCAN10", "PRDM14", "ZIC3", "DPPA4", "TDGF1"],
         "Ectoderm_1": ["PAX6", "SOX1", "OTX2", "LHX2", "SIX3", "ZIC2", "HES5"],
         "Ectoderm_2": ["PAX6", "SOX1", "OTX2", "LHX2", "RAX", "SIX6"],
         "Ectoderm_3": ["TFAP2A", "TFAP2B", "GATA3", "DLX5", "MSX1", "MSX2", "HAND1"],
         "Neural_crest": ["SOX10", "FOXD3", "TFAP2A", "TFAP2B", "SOX9", "ETS1", "TWIST1"],
         "Neuronal": ["NEUROD1", "NEUROD2", "ONECUT2", "POU3F2", "MYT1L", "INSM1", "ELAVL4", "STMN2", "NHLH1"],
         "Retinal": ["VSX2", "RAX", "SIX6", "OTX2", "LHX2", "MITF", "VAX2", "SIX3"],
         "Mesoderm": ["TWIST1", "PRRX1", "FOXF1", "HAND2", "TBX5", "GATA4", "MEIS2", "FOXC1"],
         "Endoderm": ["SOX17", "FOXA2", "FOXA1", "GATA6", "GATA4", "HNF1B", "HNF4A", "CDX2"],
         "Vascular_endothelial": ["ERG", "FLI1", "SOX18", "SOX17", "SOX7", "ETV2", "TAL1"],
         "Unknown": ["ATF3", "ATF4", "DDIT3", "JUN", "FOS", "HSF1", "EGR1", "KLF6"]}

# ---------- 3. key TFs per cell type (joint network)
key = None
if "all" in R:
    x = R["all"]
    st = x["st"]
    ztf, act, n, npos = activity(x)
    strata, tfs = list(st["strata"]), st["tf"]
    rows = []
    for s, ct in enumerate(strata):
        lvl = st["mean_tf"][s]
        ok = (n >= 30) & (lvl >= np.median(lvl))
        score = np.where(ok & (ztf[s] > 0) & (act[s] > 0), ztf[s] * act[s], 0)
        rows.append(pd.DataFrame({"celltype": ct, "tf": tfs, "expr_logcpm": lvl, "expr_z": ztf[s], "activity_z": act[s],
                                  "n_targets": n, "n_pos": npos, "score": score}))
    key = pd.concat(rows)
    key["rank"] = key.groupby("celltype").score.rank(ascending=False, method="first").astype(int)
    key.loc[key.score == 0, "rank"] = 9999
    # concordance of expression and target activity across cell types, per TF (built in by the joint fit; reported)
    r_ea = np.array([np.corrcoef(ztf[:, j], act[:, j])[0, 1] if n[j] >= 30 else np.nan for j in range(len(tfs))])
    key = key.merge(pd.DataFrame({"tf": tfs, "r_expr_activity": r_ea}), on="tf")
    mtt = mt[mt.run == "all"].set_index("tf")[["enrichment", "fdr"]].rename(columns={"enrichment": "motif_enr", "fdr": "motif_fdr"})
    key = key.join(mtt, on="tf")
    key.sort_values(["celltype", "rank"]).to_csv(f"{RES}/key_tfs_all.tsv.gz", sep="\t", index=False, float_format="%.4g")
    top = key[key["rank"] <= 10].sort_values(["celltype", "rank"])
    say("\nkey TFs (joint network), top 10 per cell type:")
    for ct in CTS:
        t = top[top.celltype == ct]
        say(f"  {ct}: " + ", ".join(f"{a}({b:.1f})" for a, b in zip(t.tf, t.score)))
    can = []
    for ct, g in CANON.items():
        k = key[key.celltype == ct].set_index("tf")
        for t in g:
            can.append(dict(celltype=ct, tf=t, tested=t in k.index, rank=k["rank"].get(t, np.nan),
                            n_targets=k.n_targets.get(t, np.nan), expr_z=k.expr_z.get(t, np.nan),
                            activity_z=k.activity_z.get(t, np.nan)))
    can = pd.DataFrame(can)
    can.to_csv(f"{RES}/canonical_tfs_all.tsv", sep="\t", index=False, float_format="%.3g")
    say("canonical TFs in top 25 / tested:", int((can["rank"] <= 25).sum()), "/", int(can.tested.sum()))
    say(can[can.tested].to_string(index=False))

# ---------- 4. within-cell-type networks
ctr = [c for c in CTS if f"{W}{c}" in R]
if ctr:
    scope = pd.concat([R[f"{W}{c}"]["sc"].assign(celltype=c) for c in ctr])
    S = scope.pivot(index="tf", columns="celltype", values="fraction")[ctr]
    S.to_csv(f"{RES}/scope_ct.tsv.gz", sep="\t", float_format="%.5g")
    lg = np.log10(S + 1e-5)
    cor = pd.DataFrame(np.nan, index=ctr, columns=ctr)
    jac = cor.copy()
    E = {c: R[f"{W}{c}"]["e"] for c in ctr}
    PK = {c: set(R[f"{W}{c}"]["st"]["peak"]) for c in ctr}
    for i, a in enumerate(ctr):
        for b in ctr[i:]:
            m = S[a].notna() & S[b].notna()
            cor.loc[a, b] = cor.loc[b, a] = spearmanr(lg[a][m], lg[b][m])[0]
            tf_s, pk_s = set(S.index[m]), PK[a] & PK[b]
            ea = E[a][E[a].tf.isin(tf_s) & E[a].peak.isin(pk_s)]
            eb = E[b][E[b].tf.isin(tf_s) & E[b].peak.isin(pk_s)]
            ka, kb = set(zip(ea.tf, ea.peak)), set(zip(eb.tf, eb.peak))
            jac.loc[a, b] = jac.loc[b, a] = len(ka & kb) / max(len(ka | kb), 1)
    cor.to_csv(f"{RES}/scope_corr_ct.tsv", sep="\t", float_format="%.4f")
    jac.to_csv(f"{RES}/edge_jaccard_ct.tsv", sep="\t", float_format="%.4f")
    say("\nwithin-CT scope Spearman:\n" + cor.round(2).to_string())
    say("edge Jaccard:\n" + jac.round(3).to_string())
    rep = []
    for c in ["iPSC", "Neuronal", "traj_neur"]:
        if c == "traj_neur" and "traj_neur_s2" in R:
            a, b = R["traj_neur"], R["traj_neur_s2"]
            sa, sb = a["sc"].set_index("tf").fraction, b["sc"].set_index("tf").fraction
            m = sa.index.intersection(sb.index)
            ka, kb = set(zip(a["e"].tf, a["e"].peak)), set(zip(b["e"].tf, b["e"].peak))
            rep.append(dict(celltype=c, scope_spearman=spearmanr(np.log10(sa[m] + 1e-5), np.log10(sb[m] + 1e-5))[0],
                            edge_jaccard=len(ka & kb) / len(ka | kb)))
            continue
        if f"{W}{c}_s2" in R and f"{W}{c}" in R:
            a, b = R[f"{W}{c}"], R[f"{W}{c}_s2"]
            sa, sb = a["sc"].set_index("tf").fraction, b["sc"].set_index("tf").fraction
            m = sa.index.intersection(sb.index)
            ka, kb = set(zip(a["e"].tf, a["e"].peak)), set(zip(b["e"].tf, b["e"].peak))
            rep.append(dict(celltype=c, scope_spearman=spearmanr(np.log10(sa[m] + 1e-5), np.log10(sb[m] + 1e-5))[0],
                            edge_jaccard=len(ka & kb) / len(ka | kb)))
    if rep:
        rep = pd.DataFrame(rep)
        rep.to_csv(f"{RES}/replicate_ct.tsv", sep="\t", index=False, float_format="%.4f")
        say("replicate ceiling:\n" + rep.to_string(index=False))
    # hubs and cell-type-biased hubs
    lr = np.log2(S + 1e-5)
    rel = lr.sub(lr.median(1), axis=0)
    hub = []
    for c in ctr:
        d = pd.DataFrame({"tf": S.index, "scope": S[c], "log2_vs_median_ct": rel[c]}).dropna()
        d["scope_rank"] = d.scope.rank(ascending=False).astype(int)
        hub.append(d.assign(celltype=c))
    hub = pd.concat(hub)
    hub.to_csv(f"{RES}/hubs_ct.tsv.gz", sep="\t", index=False, float_format="%.4g")
    say("within-CT hubs (top 10 scope) and CT-biased hubs (top-200 scope, max log2 vs median CT):")
    for c in ctr:
        d = hub[hub.celltype == c]
        say(f"  {c}: " + ", ".join(d.nsmallest(10, "scope_rank").tf) + " | biased: " +
            ", ".join(d[d.scope_rank <= 200].nlargest(8, "log2_vs_median_ct").tf))

# ---------- 5. trajectories
lin_rows, drv_all = [], {}
for lin in LINS:
    if f"traj_{lin}" not in R:
        continue
    x = R[f"traj_{lin}"]
    st = x["st"]
    ztf, act, n, npos = activity(x)
    mtf = st["mean_tf"]
    rng_e = mtf.max(0) - mtf.min(0)
    idx = np.arange(len(st["strata"]))
    r_ea = np.array([np.corrcoef(ztf[:, j], act[:, j])[0, 1] if n[j] >= 10 else np.nan for j in range(mtf.shape[1])])

    def com(M):
        w = M - M.min(0)
        return (w * idx[:, None]).sum(0) / np.maximum(w.sum(0), 1e-9)

    d = pd.DataFrame({"lineage": lin, "tf": st["tf"], "n_targets": n, "n_pos": npos, "expr_range": rng_e,
                      "expr_max": mtf.max(0), "r_expr_activity": r_ea, "com_expr": com(ztf), "com_activity": com(act),
                      "peak_bin_expr": ztf.argmax(0), "peak_bin_activity": act.argmax(0)})
    d["driver"] = (d.expr_range >= 1) & (d.n_targets >= 50) & (d.r_expr_activity >= 0.5)
    d["driver_rank"] = d.n_targets.where(d.driver, -1).rank(ascending=False, method="first").astype(int)
    lin_rows.append(d)
    drv_all[lin] = d[d.driver]
    np.savez_compressed(f"{RES}/traj_{lin}_curves.npz", tf=st["tf"], expr_z=ztf.astype(np.float32),
                        act_z=act.astype(np.float32), mean_tf=mtf, pt=st["pt_median"], comp=st["comp"],
                        comp_ct=st["comp_ct"], strata=st["strata"])
    dd = d[d.driver]
    say(f"\ntrajectory {lin}: drivers {len(dd)} of {int((d.n_targets >= 50).sum())} TFs with >= 50 targets; "
        f"median com_activity - com_expr {np.median(dd.com_activity - dd.com_expr):.2f} bins; "
        f"activity after expression in {(dd.com_activity > dd.com_expr).mean():.2f}")
    say("  early (com <= 3): " + ", ".join(dd[dd.com_expr <= 3].nlargest(10, "n_targets").tf))
    say("  late (com >= 6):  " + ", ".join(dd[dd.com_expr >= 6].nlargest(15, "n_targets").tf))
if lin_rows:
    traj = pd.concat(lin_rows)
    traj.to_csv(f"{RES}/traj_tfs.tsv.gz", sep="\t", index=False, float_format="%.4g")
    # drivers shared across lineages and target overlap
    ls_ = [l for l in LINS if l in drv_all]
    share = traj[traj.driver].groupby("tf").lineage.apply(lambda s: ",".join(sorted(s))).rename("lineages").reset_index()
    share["n_lineages"] = share.lineages.str.count(",") + 1
    share.to_csv(f"{RES}/traj_driver_sharing.tsv", sep="\t", index=False)
    say("\ndrivers per # lineages:", share.n_lineages.value_counts().sort_index().to_dict())
    say("drivers in >= 3 lineages:", ", ".join(share[share.n_lineages >= 3].tf))
    ov = []
    tgt = {l: R[f"traj_{l}"]["e"].query("coef > 0").groupby("tf").peak.apply(set) for l in ls_}
    tested = {l: set(R[f"traj_{l}"]["st"]["peak"]) for l in ls_}
    for _, s in share[share.n_lineages >= 2].iterrows():
        ll = s.lineages.split(",")
        for i, a in enumerate(ll):
            for b in ll[i + 1:]:
                U = tested[a] & tested[b]
                A, Bs = tgt[a].get(s.tf, set()) & U, tgt[b].get(s.tf, set()) & U
                if not A or not Bs:
                    continue
                inter = len(A & Bs)
                ea = len(A) * len(Bs) / len(U)
                ov.append(dict(tf=s.tf, lin_a=a, lin_b=b, n_a=len(A), n_b=len(Bs), shared=inter,
                               jaccard=inter / len(A | Bs), fold_over_random=inter / ea if ea > 0 else np.nan))
    ov = pd.DataFrame(ov)
    ov.to_csv(f"{RES}/traj_shared_driver_targets.tsv", sep="\t", index=False, float_format="%.4g")
    if len(ov):
        say("shared drivers: median target Jaccard", round(ov.jaccard.median(), 3), "| median fold over random",
            round(ov.fold_over_random.median(), 1))
        say(ov.sort_values("jaccard").head(10).to_string(index=False))
        say(ov.sort_values("jaccard").tail(10).to_string(index=False))
    # lineage-network concordance: scope Spearman on TFs tested in both
    lc = pd.DataFrame(np.nan, index=ls_, columns=ls_)
    for a in ls_:
        for b in ls_:
            sa = R[f"traj_{a}"]["sc"].set_index("tf").fraction
            sb = R[f"traj_{b}"]["sc"].set_index("tf").fraction
            m = sa.index.intersection(sb.index)
            lc.loc[a, b] = spearmanr(np.log10(sa[m] + 1e-5), np.log10(sb[m] + 1e-5))[0]
    lc.to_csv(f"{RES}/scope_corr_traj.tsv", sep="\t", float_format="%.4f")
    say("lineage scope Spearman:\n" + lc.round(2).to_string())

# ---------- 6. paired vs unpaired
pu = []
for r in runs:
    if r.endswith("_unpaired") and r.replace("_unpaired", "") in R:
        a, b = R[r.replace("_unpaired", "")], R[r]
        ka, kb = set(zip(a["e"].tf, a["e"].peak)), set(zip(b["e"].tf, b["e"].peak))
        sa, sb = a["sc"].set_index("tf").fraction, b["sc"].set_index("tf").fraction
        m = sa.index.intersection(sb.index)
        ma = mg[(mg.run == a["run"]) & (mg.sign == "pos")].enrichment
        mb = mg[(mg.run == r) & (mg.sign == "pos")].enrichment
        pu.append(dict(design=a["run"], edges_paired=len(ka), edges_unpaired=len(kb), shared=len(ka & kb),
                       frac_paired_recovered=len(ka & kb) / max(len(ka), 1),
                       peaks_with_tf_paired=a["e"].peak.nunique(), peaks_with_tf_unpaired=b["e"].peak.nunique(),
                       motif_enr_paired=ma.iloc[0] if len(ma) else np.nan, motif_enr_unpaired=mb.iloc[0] if len(mb) else np.nan,
                       scope_spearman=spearmanr(np.log10(sa[m] + 1e-5), np.log10(sb[m] + 1e-5))[0]))
if pu:
    pu = pd.DataFrame(pu)
    pu.to_csv(f"{RES}/paired_vs_unpaired.tsv", sep="\t", index=False, float_format="%.4g")
    say("\npaired vs unpaired:\n" + pu.to_string(index=False))

# ---------- 7. S-LDSC sets (hg19)
f19 = pd.read_csv(f"{BASE}/GWAS/ldsc/features_hg19_ATAC.tsv.gz", sep="\t", dtype={"chr": str}).set_index("id")
sets = []


def add(name, ids):
    d = f19.loc[f19.index.intersection(list(ids))]
    sets.append(d[["chr", "start", "end"]].assign(annot=name))
    say(f"  set {name}: {len(d)} peaks")


if key is not None:
    x = R["all"]
    st, e = x["st"], x["e"]
    det = pd.DataFrame(st["det_peak"].astype(float).T, index=st["peak"], columns=st["strata"])
    for ct in CTS:
        tops = key[(key.celltype == ct) & (key["rank"] <= 10)].tf
        pk = set(e[e.tf.isin(tops) & (e.coef > 0)].peak)
        add(f"tfprog_{ct}", [p for p in pk if det.at[p, ct] >= 0.02])
if lin_rows:
    for lin in LINS:
        if lin in drv_all:
            tops = drv_all[lin].nsmallest(10, "driver_rank").tf
            e = R[f"traj_{lin}"]["e"]
            add(f"tfdyn_{lin}", set(e[e.tf.isin(tops) & (e.coef > 0)].peak))
if sets:
    pd.concat(sets)[["annot", "chr", "start", "end"]].to_csv(f"{RES}/sets_tf.tsv.gz", sep="\t", index=False)
say("done")
