"""
caQTL covariate optimization scored by external replication in adult brain neuron / glia caQTLs
(data/BraincaQTL; see 20_brain_caqtl_replication.py) instead of pi1 alone.

Permutation mapping of every configuration is too slow on CPU, so each configuration is scored on the HDC peaks
that overlap a brain-tested peak, using the SNPs tested in both (brain window +-50 kb), genome-wide:
  HDC nominal p per shared pair: covariate-residualized correlation, t test with dof = n - 2 - n_cov (= tensorQTL)
  disc     peak's best shared-SNP p x n shared SNPs (Bonferroni within peak) < 0.01
  null     the same > 0.5
  rep      brain p < 0.05 and same direction at the HDC-best shared SNP
  r0       replication rate among null peaks (baseline from brain power / LD)
  excess   n_rep(disc) - n_disc * r0          (objective: neuron + glia excess)
HDC slope = ALT dosage (variant universe REF/ALT); brain beta = effect of the brain files' 6th-column allele.
To avoid overfitting the choice to the replication data, configurations are ranked on odd chromosomes and the
improvement is reported on even chromosomes (held out).

Grid: PhenoPC {0..10} x Sex {0,1} x GenoPC {0,1,2,3,5} x log nCell {0,1} (skip residual dof < 3).
The chosen configuration is then mapped with tensorQTL (03_map_cis.py) and checked with 20.

Output: BrainRep/covopt/summary_<CT>.tsv, BrainRep/covopt/best.tsv
Usage: python 21_covopt_brainrep.py [CellType,...]   (default all 10 cell types; ~100 GB RAM)
"""
import glob
import os
import sys
import time
from itertools import product

import numpy as np
import pandas as pd
from scipy.stats import norm, t as tdist
from tensorqtl import genotypeio, read_phenotype_bed

GENO = os.environ["GENOTYPE_DIR"] + "/Genotype/plink/chr{}"
BASE = os.environ["HDC_DIR"] + "/"
BR = os.environ["BRAIN_CAQTL_DIR"] + "/"
TAG = "noBootstrap_TMM_invnorm"
OUT = f"{BASE}/BrainRep/covopt/"
BCT = {"neuron": "Neuron_caQTL", "glia": "Glia_caQTL"}
COLS = ["c", "s", "e", "rsid", "a1", "a2", "peak", "beta", "se", "z"]
os.makedirs(OUT, exist_ok=True)

cts = sys.argv[1].split(",") if len(sys.argv) > 1 else \
    sorted(f.split("/")[-2] for f in glob.glob(f"{BASE}/caQTL/*/tensorqtl.atacExpr_{TAG}.bed.gz") if "AllCells" not in f)
U = pd.read_parquet(f"{BASE}/GWAS/variant_universe.parquet", columns=["id", "ref", "alt"]).drop_duplicates("id", keep=False).set_index("id")
bed = {k: pd.read_csv(f"{BR}/AMPAD_Data_peaks_{k}.bed", sep="\t", header=None, names=["c", "st", "en", "peak"]) for k in BCT}


def overlap_match(hp, bp):
    bp = bp.sort_values("st").reset_index(drop=True)
    st, en = bp.st.to_numpy(), bp.en.to_numpy()
    res = []
    for a, b in zip(hp.pst.to_numpy(), hp.pen.to_numpy()):
        j = np.searchsorted(st, b)
        k = np.arange(max(j - 20, 0), j)
        ov = np.minimum(en[k], b) - np.maximum(st[k], a)
        res.append(bp.peak.iat[k[np.argmax(ov)]] if len(k) and ov.max() > 0 else None)
    return np.array(res, dtype=object)


# ---------------- genotypes and brain stats, loaded once ----------------
t0 = time.time()
geno, brain = {}, {}
for chrom in map(str, range(1, 23)):
    gdf, _ = genotypeio.load_genotypes(GENO.format(chrom))
    gdf = gdf[~gdf.index.duplicated(keep=False)]  # duplicated rsIDs (multi-allelic) are ambiguous
    gdf = gdf[gdf.index.isin(U.index)]
    for bct, d in BCT.items():
        b = pd.read_csv(f"{BR}/{d}/combined_results_chr{chrom}_sorted.gz", sep="\t", header=None, names=COLS,
                        usecols=["rsid", "a1", "a2", "peak", "beta", "se", "z"])
        b = b[(b.se > 0) & b.rsid.isin(gdf.index)]
        al = U.reindex(b.rsid)
        ref, alt = al.ref.to_numpy(), al.alt.to_numpy()
        a1, a2 = b.a1.to_numpy(), b.a2.to_numpy()
        fwd, swp = (a1 == ref) & (a2 == alt), (a1 == alt) & (a2 == ref)
        b = b.assign(brain_beta_alt=np.where(fwd, b.beta, -b.beta), brain_p=2 * norm.sf(np.abs(b.z)))[fwd | swp]
        brain[(chrom, bct)] = b[["rsid", "peak", "brain_beta_alt", "brain_p"]].rename(columns={"rsid": "variant_id", "peak": "brain_peak"})
    geno[chrom] = gdf[gdf.index.isin(set(brain[(chrom, "neuron")].variant_id) | set(brain[(chrom, "glia")].variant_id))]
    print(f"chr{chrom}: {len(geno[chrom])} shared-candidate variants; {time.time() - t0:.0f}s", flush=True)

best_rows = []
for ct in cts:
    t1 = time.time()
    P_df, pos = read_phenotype_bed(f"{BASE}/caQTL/{ct}/tensorqtl.atacExpr_{TAG}.bed.gz")
    n = P_df.shape[1]
    base_cov = pd.read_csv(f"{BASE}/caQTL/{ct}/tensorqtl.atacCov_{TAG}.txt", sep="\t", index_col=0)[P_df.columns]
    nc = pd.read_csv(f"{BASE}/PseudoBulk/{ct}_NCell.tsv", sep="\t", index_col=0).iloc[:, 0]
    lognc = np.log(nc.reindex(P_df.columns).astype(float)).to_numpy()
    w, v = np.linalg.eigh(np.cov(P_df.values, rowvar=False))
    V = v[:, np.argsort(w)[::-1]]
    sh = []
    for chrom in map(str, range(1, 23)):
        sel = pos["chr"] == chrom
        if not sel.any():
            continue
        hp = pd.DataFrame({"phenotype_id": P_df.index[sel]})
        x = hp.phenotype_id.str.split("_", expand=True)
        hp["pst"], hp["pen"] = x[1].astype(int), x[2].astype(int)
        for bct in BCT:
            b = brain[(chrom, bct)]
            g = bed[bct][(bed[bct].c == f"chr{chrom}") & bed[bct].peak.isin(b.brain_peak.unique())]
            mp = pd.DataFrame({"phenotype_id": hp.phenotype_id, "brain_peak": overlap_match(hp, g)}).dropna()
            sh.append(mp.merge(b, on="brain_peak").assign(brain_ct=bct, chr=int(chrom)))
    sh = pd.concat(sh, ignore_index=True)
    vid = pd.Index(sh.variant_id.unique())
    G = pd.concat([geno[c][P_df.columns] for c in geno]).reindex(vid).to_numpy(float)
    G = np.where(np.isnan(G), np.nanmean(G, 1, keepdims=True), G)
    keepv = G.std(1) > 0  # monomorphic in these donors
    sh = sh[keepv[vid.get_indexer(sh.variant_id)]].reset_index(drop=True)
    pid = pd.Index(sh.phenotype_id.unique())
    P = P_df.loc[pid].to_numpy(float)
    ip, iv = pid.get_indexer(sh.phenotype_id), vid.get_indexer(sh.variant_id)
    sh["nshared"] = sh.groupby(["phenotype_id", "brain_ct"]).variant_id.transform("size")
    grp = sh.groupby(["phenotype_id", "brain_ct"]).ngroup().to_numpy()
    print(f"{ct}: {len(sh)} shared pairs, {sh.groupby('brain_ct').phenotype_id.nunique().to_dict()} peaks; {time.time() - t1:.0f}s", flush=True)

    rows = []
    for k, sex, gpc, ncell in product(range(11), (0, 1), (0, 1, 2, 3, 5), (0, 1)):
        C = [V[:, :k]] + ([base_cov.loc[["Sex"]].to_numpy(float).T] if sex else []) + \
            ([base_cov.loc[[f"GenoPC_{i + 1}" for i in range(gpc)]].to_numpy(float).T] if gpc else []) + \
            ([lognc[:, None]] if ncell else [])
        X = np.column_stack([np.ones(n)] + C)
        ncov = X.shape[1] - 1
        dof = n - 2 - ncov
        if dof < 3:
            continue
        Q, _ = np.linalg.qr(X)
        rg = G - (G @ Q) @ Q.T
        rp = P - (P @ Q) @ Q.T
        rg /= np.linalg.norm(rg, axis=1, keepdims=True)
        rp /= np.linalg.norm(rp, axis=1, keepdims=True)
        r = np.clip(np.einsum("ij,ij->i", rg[iv], rp[ip]), -0.999999, 0.999999)
        p = 2 * tdist.sf(np.abs(r * np.sqrt(dof / (1 - r ** 2))), dof)
        o = np.lexsort((p, grp))
        first = np.r_[True, grp[o][1:] != grp[o][:-1]]
        best = sh.iloc[o[first]].assign(hdc_r=r[o[first]], hdc_p=p[o[first]])
        best["padj"] = np.minimum(best.hdc_p * best.nshared, 1)
        best["rep"] = (best.brain_p < 0.05) & (np.sign(best.hdc_r) == np.sign(best.brain_beta_alt))
        out = dict(celltype=ct, nPhenoPC=k, Sex=sex, nGenoPC=gpc, logNCell=ncell, nCov=ncov, dof=dof)
        for half, hd in (("odd", best[best.chr % 2 == 1]), ("even", best[best.chr % 2 == 0])):
            for bct, d in hd.groupby("brain_ct"):
                disc, nul = d[d.padj < 0.01], d[d.padj > 0.5]
                r0 = nul.rep.mean()
                out.update({f"{half}_{bct}_nDisc": len(disc), f"{half}_{bct}_rep": int(disc.rep.sum()),
                            f"{half}_{bct}_r0": r0, f"{half}_{bct}_excess": disc.rep.sum() - len(disc) * r0,
                            f"{half}_{bct}_repRate": disc.rep.mean()})
            out[f"{half}_obj"] = out[f"{half}_neuron_excess"] + out[f"{half}_glia_excess"]
        rows.append(out)
    d = pd.DataFrame(rows)
    d.to_csv(f"{OUT}/summary_{ct}.tsv", sep="\t", index=False, float_format="%.4g")
    sel = d.sort_values("odd_obj", ascending=False).iloc[0]
    cur = d[(d.nPhenoPC == 3) & (d.Sex == 0) & (d.nGenoPC == 0) & (d.logNCell == 0)].iloc[0]  # production optCov
    best_rows.append(dict(celltype=ct, **{f"sel_{c}": sel[c] for c in ["nPhenoPC", "Sex", "nGenoPC", "logNCell"]},
                          sel_odd_obj=sel.odd_obj, cur_odd_obj=cur.odd_obj, sel_even_obj=sel.even_obj, cur_even_obj=cur.even_obj,
                          sel_even_neuron_rep=sel.even_neuron_rep, cur_even_neuron_rep=cur.even_neuron_rep,
                          sel_even_glia_rep=sel.even_glia_rep, cur_even_glia_rep=cur.even_glia_rep,
                          rank_cur_odd=int((d.odd_obj > cur.odd_obj).sum()) + 1, n_configs=len(d)))
    print(ct, f"{len(d)} configs, {time.time() - t1:.0f}s; selected", sel[["nPhenoPC", "Sex", "nGenoPC", "logNCell"]].to_dict(),
          f"odd {sel.odd_obj:.1f} vs current {cur.odd_obj:.1f}; held-out even {sel.even_obj:.1f} vs current {cur.even_obj:.1f}", flush=True)
    pd.DataFrame(best_rows).to_csv(f"{OUT}/best.tsv", sep="\t", index=False, float_format="%.4g")
print(pd.DataFrame(best_rows).to_string())
