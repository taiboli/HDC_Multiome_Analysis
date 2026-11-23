"""
Aggregate (locus-free) test of trait relevance of each cell type x modality: are lead variants of the
strongest HDC QTLs enriched for GWAS association compared with lead variants of null features?
This pools information over thousands of features, so it is informative even where few QTLs pass FDR.

For each cell type x modality (lead variants from 07_coloc_smr.py lead tables; one feature per lead variant,
the one with smallest pval_beta) and each GWAS trait:
  strong = pval_beta < 0.01, null = pval_beta > 0.5
  (1) chi2 enrichment: mean GWAS chi2 at strong leads / mean chi2 at null leads re-weighted to the strong
      leads' composition of lead-variant MAF (in the 17 donors) x distance-to-feature strata (direct standardization)
  (2) Mantel-Haenszel OR for GWAS p < 1e-3 at the lead variant, same strata (reported only when >= 5 events in
      each group; rare-event ORs at stricter thresholds were unstable)
  Standard errors: delete-one-block jackknife over 100 contiguous genomic blocks (LD / clustering of GWAS loci).
Pooled rows "ALL_<MOD>" combine all cell types (union of lead variants, best pval_beta per variant).

Output: GWAS/results/enrichment_qtl_leads.tsv (+ BH FDR across all tests)
"""
import glob
import os

import numpy as np
import pandas as pd
from scipy.stats import norm

BASE = os.environ["HDC_DIR"] + "/"
R = f"{BASE}/GWAS/results/"
os.makedirs(R, exist_ok=True)
man = pd.read_csv(f"{BASE}/GWAS/gwas_manifest.tsv", sep="\t")[["trait_id", "source", "ancestry", "label"]]
Z_SIG = norm.isf(1e-3 / 2)
MIN_EVENTS = 5
NBLOCK = 100


def mh(a, b, c, d):
    n = a + b + c + d
    with np.errstate(divide="ignore", invalid="ignore"):
        num, den = np.nansum(a * d / n), np.nansum(b * c / n)
    return num / den if den > 0 else np.nan


def table(df, z):
    """per-stratum 2x2: a strong&sig, b strong&not, c null&sig, d null&not"""
    sig = (z * z > Z_SIG ** 2)
    s = df.strong.to_numpy()
    g = df.stratum.to_numpy()
    k = g.max() + 1
    a = np.bincount(g, (s & sig), k)
    b = np.bincount(g, (s & ~sig), k)
    c = np.bincount(g, (~s & sig), k)
    d = np.bincount(g, (~s & ~sig), k)
    return a, b, c, d


def chi2_ratio(df, z):
    """mean chi2 strong / stratum-standardized mean chi2 null"""
    s = df.strong.to_numpy()
    g = df.stratum.to_numpy()
    k = g.max() + 1
    c2 = z * z
    ns, nn = np.bincount(g, s, k), np.bincount(g, ~s, k)
    ss, sn = np.bincount(g, c2 * s, k), np.bincount(g, c2 * ~s, k)
    ok = (ns > 0) & (nn > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        return ss[ok].sum() / ns[ok].sum() / (np.sum(ns[ok] * sn[ok] / nn[ok]) / ns[ok].sum())


def jk_se(vals):
    v = np.asarray(vals, float)
    v = v[np.isfinite(v)]
    if len(v) < NBLOCK * 0.9:
        return np.nan
    return np.sqrt((len(v) - 1) / len(v) * np.sum((v - v.mean()) ** 2))


def analyze(lead, label, mod):
    lead = lead.sort_values("pval_beta").drop_duplicates("variant_id")
    lead = lead[(lead.pval_beta < 0.01) | (lead.pval_beta > 0.5)].copy()
    lead["strong"] = lead.pval_beta < 0.01
    maf = np.minimum(lead.af, 1 - lead.af)
    mb = np.digitize(maf, [0.1, 0.2, 0.3])
    db = np.digitize(lead.start_distance.abs(), [10_000, 50_000, 200_000])
    lead["stratum"] = mb * 4 + db
    lead = lead.sort_values(["chr", "pos38"])
    lead["block"] = np.arange(len(lead)) * NBLOCK // len(lead)
    rows = []
    traits = [c for c in lead.columns if "__" in c]
    for t in traits:
        d = lead[np.isfinite(lead[t])]
        z = d[t].to_numpy(float)
        if len(d) < 200 or d.strong.sum() < 20:
            continue
        a, b, c, dd = table(d, z)
        rr = chi2_ratio(d, z)
        orr = mh(a, b, c, dd) if min(a.sum(), c.sum()) >= MIN_EVENTS else np.nan
        jr, jo = [], []
        blk = d.block.to_numpy()
        for k in range(NBLOCK):
            m = blk != k
            jr.append(np.log(chi2_ratio(d[m], z[m])))
            if np.isfinite(orr):
                aa, bb, cc, ddd = table(d[m], z[m])
                with np.errstate(divide="ignore"):
                    jo.append(np.log(mh(aa, bb, cc, ddd)))
        se_r = jk_se(jr)
        se_o = jk_se(jo) if np.isfinite(orr) else np.nan
        lr = np.log(rr)
        lo = np.log(orr) if np.isfinite(orr) and orr > 0 else np.nan
        s = d.strong.to_numpy()
        rows.append(dict(celltype=label, modality=mod, trait_id=t, n_strong=int(s.sum()), n_null=int((~s).sum()),
                         mean_chi2_strong=np.mean(z[s] ** 2), mean_chi2_null=np.mean(z[~s] ** 2),
                         chi2_ratio=rr, log_chi2_ratio=lr, se_jk=se_r,
                         p=norm.sf(lr / se_r) if se_r > 0 else np.nan,  # one-sided: enrichment
                         sig_strong=int(a.sum()), sig_null=int(c.sum()), OR_MH_p1e3=orr, se_jk_logOR=se_o,
                         p_OR=2 * norm.sf(abs(lo / se_o)) if np.isfinite(lo) and se_o > 0 else np.nan))
    return rows


out = []
for mod in ("RNA", "ATAC"):
    pooled = []
    for ct_dir in sorted(glob.glob(f"{BASE}/GWAS/coloc/{mod}/*/")):
        ct = ct_dir.rstrip("/").split("/")[-1]
        lead_files = sorted(glob.glob(f"{ct_dir}/chr*.lead.parquet"))
        if not lead_files:
            continue
        lead = pd.concat([pd.read_parquet(f).assign(chr=int(os.path.basename(f).split(".")[0][3:]))
                          for f in lead_files], ignore_index=True)
        out += analyze(lead, ct, mod)
        pooled.append(lead)
        print(mod, ct, flush=True)
    out += analyze(pd.concat(pooled, ignore_index=True), f"ALL_{mod}", mod)

res = pd.DataFrame(out).merge(man, on="trait_id", how="left")
p = res.p.fillna(1).to_numpy()
o = np.argsort(p)
q = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
res["FDR"] = np.nan
res.loc[res.index[o], "FDR"] = np.minimum(q, 1)
res.sort_values("p").to_csv(f"{R}/enrichment_qtl_leads.tsv", sep="\t", index=False, float_format="%.4g")
print(res.sort_values("p").head(30).to_string())
