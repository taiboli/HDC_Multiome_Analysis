"""
External replication of HDC caQTLs in adult brain neuron / glia caQTLs (FANS-sorted, data/BraincaQTL;
Synapse syn26708632: nominal stats for all SNP-peak pairs within +-50 kb, hg38; peak coordinates only provided
for FDR 5% caPeaks).

1. Brain peak positions: AMPAD_Data_peaks_<ct>.bed (hg38; the peak IDs of the caQTL files -- every tested peak is
   in it and the FDR 5% peaks of position_for_significant_peaks_FDR0.05_<ct>.RDS, exported to
   BrainRep/sigpeaks_<ct>.bed, have identical coordinates).
2. HDC peak <-> brain peak: interval overlap among brain peaks tested for caQTL (largest overlap if several).
3. HDC discovery sets per cell type: tensorQTL permutation leads with pval_beta < 0.01 ("perm_p01"; genome-wide
   FDR is ~0 with 17 donors) and mashr lfsr < 0.05 top variants ("mash"; 14_mash.R); "mash_specific" = lfsr <
   0.05 here and in <= 20% of the other tested cell types. Background: permutation leads with pval_beta > 0.5
   ("perm_null"), mash features with lfsr > 0.5 in all cell types ("mash_null").
4. Replication SNP: (main) "best_shared_snp" = the SNP with the smallest HDC nominal p (in that cell type) among the
   SNPs also tested for the matched brain peak (brain window is +-50 kb vs 1 Mb in HDC, so most HDC leads are
   absent; choosing on HDC p only keeps the brain p unbiased); (sensitivity) "hdc_lead" = the permutation lead /
   mash top variant itself. Brain p-value of that SNP for the matched brain peak; pi1 = 1 - pi0 (Storey qvalue, R)
   with 100-bootstrap 95% CI; sign concordance among brain p < 0.05 (HDC slope = ALT allele; brain beta = effect
   of its 6th column allele, checked separately for REF/ALT-ordered and swapped rows).
5. Overlap enrichment: HDC peak is a discovery vs tested non-discovery, outcome = overlaps a brain FDR 5% caPeak,
   universe = HDC peaks matched to a tested brain peak (controls for accessibility in brain); Fisher OR and a
   logistic regression adjusting for log mean pseudobulk CPM (AllCells) and peak width.

Outputs (BrainRep/): pairs_{best_shared_snp,hdc_lead}.tsv.gz, pi1.tsv, overlap_enrichment.tsv, brain_peaks.tsv.gz
Usage: python 20_brain_caqtl_replication.py   (needs `module load R/4.2.0; export R_HOME=$(R RHOME)`)
"""
import glob
import os
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, norm
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
from tensorqtl import rfunc

BASE = os.environ["HDC_DIR"] + "/"
BR = os.environ["BRAIN_CAQTL_DIR"] + "/"
OUT = f"{BASE}/BrainRep/"
TAG = "noBootstrap_TMM_invnorm"
BCT = {"neuron": "Neuron_caQTL", "glia": "Glia_caQTL"}
BED = {k: pd.read_csv(f"{BR}/AMPAD_Data_peaks_{k}.bed", sep="\t", header=None, names=["c", "st", "en", "peak"]) for k in BCT}
COLS = ["c", "s", "e", "rsid", "a1", "a2", "peak", "beta", "se", "z"]
PAD = 500
NBOOT = 100
NNULL = 5000
np.random.seed(1)


def safe_pi1(p):
    try:
        return 1 - rfunc.qvalue(p)[1]
    except Exception:  # lambda grid beyond max(p): fall back to lambda = 0.5
        return 1 - rfunc.qvalue(p, lambda_qvalue=0.5)[1]


def logit(y, X, it=50):
    """IRLS logistic regression; returns coefficients and SEs"""
    b = np.zeros(X.shape[1])
    for _ in range(it):
        mu = 1 / (1 + np.exp(-X @ b))
        W = mu * (1 - mu)
        H = X.T @ (X * W[:, None])
        step = np.linalg.solve(H, X.T @ (y - mu))
        b += step
        if np.max(np.abs(step)) < 1e-8:
            break
    return b, np.sqrt(np.diag(np.linalg.inv(H)))


def pid_coords(ids):
    x = pd.Series(ids).str.split("_", expand=True)
    return x[0].str[3:].astype(int).to_numpy(), x[1].astype(int).to_numpy(), x[2].astype(int).to_numpy()


# ---------------- HDC discovery / background sets ----------------
sets = []
for f in sorted(glob.glob(f"{BASE}/caQTL_Output_optCov/*/allchr_cis_{TAG}.cis_qtl.qval.txt.gz")):
    ct = f.split("/")[-2]
    q = pd.read_csv(f, sep="\t", usecols=["phenotype_id", "variant_id", "pval_beta", "slope", "af"])
    sets.append(q[q.pval_beta < 0.01].assign(celltype=ct, set="perm_p01"))
    sets.append(q[q.pval_beta > 0.5].assign(celltype=ct, set="perm_null"))
m = pd.read_csv(f"{BASE}/mash/ATAC_mash_strong.tsv.gz", sep="\t")
mcts = [c[5:] for c in m.columns if c.startswith("lfsr_")]
LF = m[[f"lfsr_{c}" for c in mcts]].to_numpy()
tested = np.isfinite(LF)
sig = (LF < 0.05) & tested
null = np.all(~tested | (LF > 0.5), axis=1)
for i, ct in enumerate(mcts):
    other = (sig.sum(1) - sig[:, i]) / np.maximum(tested.sum(1) - tested[:, i], 1)
    base = m[["phenotype_id", "variant_id"]].assign(slope=m[f"pm_{ct}"], celltype=ct)
    sets.append(base[sig[:, i]].assign(set="mash"))
    sets.append(base[sig[:, i] & (other <= 0.2)].assign(set="mash_specific"))
    sets.append(base[tested[:, i] & null].assign(set="mash_null"))
sets = pd.concat(sets, ignore_index=True)
# cap background sets at NNULL random peaks per cell type
rng = np.random.default_rng(1)
keep = np.ones(len(sets), bool)
for (ct, st), idx in sets.groupby(["celltype", "set"]).groups.items():
    if st.endswith("null") and len(idx) > NNULL:
        keep[np.setdiff1d(idx, rng.choice(idx, NNULL, replace=False))] = False
sets = sets[keep].reset_index(drop=True)
NOM = {ct: f"{BASE}/caQTL_Output_optCov/{ct}/chr{{ch}}_cis_nominal_{TAG}.cis_qtl_pairs.{{ch}}.parquet"
       for ct in sets.celltype.unique()}
sets["chr"], sets["pst"], sets["pen"] = pid_coords(sets.phenotype_id)
print(sets.groupby(["set", "celltype"]).size().unstack(0).to_string(), flush=True)

# all tested HDC peaks (universe for overlap enrichment)
hdc = {}
for f in sorted(glob.glob(f"{BASE}/caQTL_Output_optCov/*/allchr_cis_{TAG}.cis_qtl.qval.txt.gz")):
    ct = f.split("/")[-2]
    hdc[ct] = pd.read_csv(f, sep="\t", usecols=["phenotype_id", "pval_beta"])
allpk = pd.Series(sorted(set().union(*[set(d.phenotype_id) for d in hdc.values()])))
apk = pd.DataFrame({"phenotype_id": allpk})
apk["chr"], apk["pst"], apk["pen"] = pid_coords(apk.phenotype_id)


def match_peaks(hp, bp):
    """hp: HDC peaks (pst pen) on one chr; bp: brain peaks (st en). Brain peak with the largest overlap."""
    bp = bp.sort_values("st").reset_index(drop=True)
    st, en = bp.st.to_numpy(), bp.en.to_numpy()
    res = []
    for a, b in zip(hp.pst.to_numpy(), hp.pen.to_numpy()):
        j = np.searchsorted(st, b)  # brain peaks starting before the HDC end
        k = np.arange(max(j - 20, 0), j)
        ov = np.minimum(en[k], b) - np.maximum(st[k], a)
        res.append(bp.peak.iat[k[np.argmax(ov)]] if len(k) and ov.max() > 0 else None)
    return np.array(res, dtype=object)


def per_chr(ch):
    out, peaks, acc, best = [], [], [], []
    hp = apk[apk.chr == ch]
    s = sets[sets.chr == ch]
    want = s[["phenotype_id", "variant_id"]].drop_duplicates()
    for bct, d in BCT.items():
        if not os.path.exists(f"{BR}/{d}/combined_results_chr{ch}_sorted.gz"):  # neuron chr3/6/12 not downloaded
            print("missing", bct, ch, flush=True)
            continue
        b = pd.read_csv(f"{BR}/{d}/combined_results_chr{ch}_sorted.gz", sep="\t", header=None, names=COLS,
                        usecols=["e", "rsid", "a1", "a2", "peak", "beta", "se", "z"])
        sp = pd.read_csv(f"{OUT}/sigpeaks_{bct}.bed", sep="\t", header=None, names=["c", "st", "en", "peak"])
        sp = sp[sp.c == f"chr{ch}"].set_index("peak")
        g = BED[bct]
        g = g[(g.c == f"chr{ch}") & g.peak.isin(b.peak.unique())].copy()
        g["sig"] = g.peak.isin(sp.index)
        acc.append(pd.DataFrame({"brain_ct": [bct], "chr": [ch], "n_bed": [(BED[bct].c == f"chr{ch}").sum()],
                                 "n_tested": [b.peak.nunique()], "n_tested_in_bed": [len(g)], "n_sig": [len(sp)]}))
        peaks.append(g.assign(brain_ct=bct, chr=ch))
        mp = pd.DataFrame({"phenotype_id": hp.phenotype_id.to_numpy(), "brain_peak": match_peaks(hp, g)}).dropna()
        mp["brain_sig"] = mp.brain_peak.isin(sp.index)
        peaks[-1] = (peaks[-1], mp.assign(brain_ct=bct))
        b = b[b.se > 0].rename(columns={"peak": "brain_peak", "rsid": "variant_id"})
        # (a) the HDC lead / mash top variant itself (only present if within the brain +-50 kb window)
        w = want.merge(mp, on="phenotype_id")
        w = w.merge(b, on=["brain_peak", "variant_id"], how="left")
        out.append(w.assign(brain_ct=bct))
        # (b) HDC-best SNP (min nominal p in that cell type) among SNPs also tested for the matched brain peak
        bb = b.merge(mp[mp.phenotype_id.isin(s.phenotype_id)], on="brain_peak")
        keys = bb[["phenotype_id", "variant_id"]].drop_duplicates()
        for ct, nom in NOM.items():
            ids = s.loc[s.celltype == ct, "phenotype_id"].unique()
            if not len(ids) or not os.path.exists(nom.format(ch=ch)):
                continue
            pf = pq.ParquetFile(nom.format(ch=ch))
            idset = pa.array(ids)
            parts = []
            for rb in pf.iter_batches(batch_size=2_000_000, columns=["phenotype_id", "variant_id", "pval_nominal", "slope"]):
                rb = rb.filter(pc.is_in(rb.column("phenotype_id"), value_set=idset))
                if rb.num_rows:
                    parts.append(rb.to_pandas())
            if not parts:
                continue
            h = pd.concat(parts).merge(keys, on=["phenotype_id", "variant_id"])
            h = h.sort_values("pval_nominal").drop_duplicates("phenotype_id")
            h = h.merge(bb, on=["phenotype_id", "variant_id"], how="left")
            best.append(h.rename(columns={"pval_nominal": "hdc_p", "slope": "hdc_slope"}).assign(celltype=ct, brain_ct=bct))
    return pd.concat(out), peaks, pd.concat(acc), pd.concat(best) if best else pd.DataFrame()


if __name__ == "__main__":
    with Pool(int(os.environ.get("NPROC", 11))) as pool:
        res = pool.map(per_chr, range(1, 23))
    pairs = pd.concat([r[0] for r in res], ignore_index=True)
    match = pd.concat([p[1] for r in res for p in r[1]], ignore_index=True)
    bpk = pd.concat([p[0] for r in res for p in r[1]], ignore_index=True)
    acc = pd.concat([r[2] for r in res], ignore_index=True)
    bpk.to_csv(f"{OUT}/brain_peaks.tsv.gz", sep="\t", index=False)
    accs = acc.groupby("brain_ct")[["n_bed", "n_tested", "n_tested_in_bed", "n_sig"]].sum()
    accs.to_csv(f"{OUT}/brain_peaks_summary.tsv", sep="\t")
    print(accs.to_string(), flush=True)

    # allele alignment: HDC ALT vs brain alleles
    u = pd.read_parquet(f"{BASE}/GWAS/variant_universe.parquet", columns=["id", "ref", "alt"]).drop_duplicates("id", keep=False)
    u = u.rename(columns={"id": "variant_id"})

    def align(x):
        x = x.merge(u, on="variant_id", how="left")
        x["orient"] = np.select([(x.a1 == x.ref) & (x.a2 == x.alt), (x.a1 == x.alt) & (x.a2 == x.ref)], ["ref_alt", "swapped"], "na")
        x["brain_beta_alt"] = np.where(x.orient == "ref_alt", x.beta, np.where(x.orient == "swapped", -x.beta, np.nan))
        x["brain_p"] = 2 * norm.sf(np.abs(x.z))
        return x.drop(columns=["ref", "alt"])

    best = pd.concat([r[3] for r in res], ignore_index=True)
    full = {"best_shared_snp": sets.drop(columns=["variant_id", "slope"]).merge(align(best), on=["celltype", "phenotype_id"], how="left")
                                  .rename(columns={"hdc_slope": "slope"}),
            "hdc_lead": sets.merge(align(pairs), on=["phenotype_id", "variant_id"], how="left")}
    for k, x in full.items():
        x.to_csv(f"{OUT}/pairs_{k}.tsv.gz", sep="\t", index=False, float_format="%.5g")

    # sign concordance by orientation (strong HDC + brain pairs) to fix the brain effect allele
    x = full["best_shared_snp"]
    st = x[(x.hdc_p < 1e-3) & (x.brain_p < 1e-3)].drop_duplicates(["phenotype_id", "variant_id", "brain_ct"])
    for o, d in st.groupby("orient"):
        print("orientation", o, len(d), "concordance (brain beta = 6th col allele):",
              round(float((np.sign(d.slope) == np.sign(d.brain_beta_alt)).mean()), 3), flush=True)

    rows = []
    for k, x in full.items():
        x = x[x.brain_p.notna()]
        for (ct, st_, bct), d in x.groupby(["celltype", "set", "brain_ct"]):
            d = d.drop_duplicates("phenotype_id")
            p = d.brain_p.to_numpy()
            r = dict(snp=k, celltype=ct, set=st_, brain_ct=bct, n_pairs=len(p))
            if len(p) >= 30:
                bs = [safe_pi1(np.random.choice(p, len(p))) for _ in range(NBOOT)]
                dd = d[(d.brain_p < 0.05) & d.brain_beta_alt.notna()]
                r.update(pi1=safe_pi1(p), pi1_lo=np.nanpercentile(bs, 2.5), pi1_hi=np.nanpercentile(bs, 97.5),
                         frac_p05=(p < 0.05).mean(), n_p05=len(dd),
                         sign_conc_p05=float((np.sign(dd.slope) == np.sign(dd.brain_beta_alt)).mean()),
                         frac_brain_sig_peak=float(d.brain_sig.mean()))
            rows.append(r)
    pi = pd.DataFrame(rows)
    pi.to_csv(f"{OUT}/pi1.tsv", sep="\t", index=False, float_format="%.4g")
    print(pi.round(3).to_string(), flush=True)

    # overlap enrichment
    cnt = pd.read_csv(f"{BASE}/PseudoBulk/AllCells_ATAC_SumCount.tsv.gz", sep="\t", index_col=0)
    lcpm = np.log1p((cnt / cnt.sum() * 1e6).mean(1)).rename("lcpm")
    rows = []
    for bct in BCT:
        mm = match[match.brain_ct == bct].set_index("phenotype_id").brain_sig
        for ct, q in hdc.items():
            for name, disc in [("perm_p01", q.set_index("phenotype_id").pval_beta < 0.01),
                               ("mash", None), ("mash_specific", None)]:
                if disc is None:
                    if ct not in mcts:
                        continue
                    ids = set(sets[(sets.celltype == ct) & (sets.set == name)].phenotype_id)
                    disc = pd.Series(q.phenotype_id.isin(ids).to_numpy(), index=q.phenotype_id)
                d = pd.DataFrame({"disc": disc.astype(int)}).join(mm, how="inner").join(lcpm, how="left").dropna()
                _, ps, pe = pid_coords(d.index)
                d["lw"] = np.log(pe - ps)
                d["brain_sig"] = d.brain_sig.astype(int)
                t = pd.crosstab(d.disc, d.brain_sig).reindex(index=[0, 1], columns=[0, 1], fill_value=0)
                orr, pf = fisher_exact(t.to_numpy(), alternative="greater")
                try:
                    X = np.column_stack([np.ones(len(d)), d[["disc", "lcpm", "lw"]].to_numpy(float)])
                    cf, se = logit(d.brain_sig.to_numpy(float), X)
                    lor, lse = cf[1], se[1]
                    lp = 2 * norm.sf(abs(lor / lse))
                except Exception:
                    lor = lse = lp = np.nan
                rows.append(dict(brain_ct=bct, celltype=ct, set=name, n_universe=len(d), n_disc=int(d.disc.sum()),
                                 frac_overlap_disc=d.brain_sig[d.disc == 1].mean(), frac_overlap_bg=d.brain_sig[d.disc == 0].mean(),
                                 fisher_OR=orr, fisher_p=pf, adj_OR=np.exp(lor), adj_OR_lo=np.exp(lor - 1.96 * lse),
                                 adj_OR_hi=np.exp(lor + 1.96 * lse), adj_p=lp))
    ov = pd.DataFrame(rows)
    ov.to_csv(f"{OUT}/overlap_enrichment.tsv", sep="\t", index=False, float_format="%.4g")
    print(ov.round(4).to_string())
