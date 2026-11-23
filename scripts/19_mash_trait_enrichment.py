"""
Trait enrichment of mashr QTLs (14_mash.R): does borrowing strength across the 10 cell types give QTLs whose
variants carry more GWAS signal than null features?

Per modality, each feature has one top cross-cell-type variant (13_mash_extract.py); its GWAS z for every trait
is taken from the 07_coloc_smr.py lead tables (the variant is the permutation lead of at least one cell type).
Sets per cell type:
  sig      lfsr < 0.05 in the cell type
  specific sig in the cell type and significant in <= 20% of the other cell types in which it was tested
  shared   sig in >= 80% of the tested cell types
  null     lfsr > 0.5 in every tested cell type (common background for all sets)
Statistic: mean GWAS chi2 in the set / mean chi2 of null features standardized to the set's composition of
top-variant MAF x |distance to feature| strata; SE by delete-one-block jackknife (100 contiguous blocks);
one-sided p. Pooled "anyCT" row: sig in any cell type. BH FDR over all tests.

Output: GWAS/results/mash_trait_enrichment.tsv
"""
import glob
import os

import numpy as np
import pandas as pd
from scipy.stats import norm

BASE = os.environ["HDC_DIR"] + "/"
R = f"{BASE}/GWAS/results/"
NBLOCK = 100
man = pd.read_csv(f"{BASE}/GWAS/gwas_manifest.tsv", sep="\t")[["trait_id", "source", "ancestry", "label"]]


def ratio(z2, s, g, k):
    ns, nn = np.bincount(g, s, k), np.bincount(g, ~s, k)
    ss, sn = np.bincount(g, z2 * s, k), np.bincount(g, z2 * ~s, k)
    ok = (ns > 0) & (nn > 0)
    if ns[ok].sum() == 0:
        return np.nan
    return (ss[ok].sum() / ns[ok].sum()) / (np.sum(ns[ok] * sn[ok] / nn[ok]) / ns[ok].sum())


def test(d, s, traits):
    """d: features (sorted by position) with stratum, block, trait z columns; s: boolean set membership"""
    k = d.stratum.max() + 1
    g, blk = d.stratum.to_numpy(), d.block.to_numpy()
    out = []
    for t in traits:
        z = d[t].to_numpy(float)
        ok = np.isfinite(z)
        if s[ok].sum() < 20:
            continue
        z2, ss, gg, bb = z[ok] ** 2, s[ok], g[ok], blk[ok]
        r = ratio(z2, ss, gg, k)
        jk = np.array([np.log(ratio(z2[bb != b], ss[bb != b], gg[bb != b], k)) for b in range(NBLOCK)])
        jk = jk[np.isfinite(jk)]
        se = np.sqrt((len(jk) - 1) / len(jk) * np.sum((jk - jk.mean()) ** 2)) if len(jk) > 0.9 * NBLOCK else np.nan
        lr = np.log(r)
        out.append(dict(trait_id=t, n_set=int(ss.sum()), n_null=int((~ss).sum()), chi2_ratio=r, log_ratio=lr,
                        se_jk=se, p=norm.sf(lr / se) if se > 0 else np.nan,
                        mean_chi2_set=z2[ss].mean(), mean_chi2_null=z2[~ss].mean()))
    return out


rows = []
for mod in ("RNA", "ATAC"):
    m = pd.read_csv(f"{BASE}/mash/{mod}_mash_strong.tsv.gz", sep="\t")
    cts = [c[5:] for c in m.columns if c.startswith("lfsr_")]
    # GWAS z / af / distance of the top variant from the lead tables of any cell type
    lead = []
    for ct in cts:
        for f in glob.glob(f"{BASE}/GWAS/coloc/{mod}/{ct}/chr*.lead.parquet"):
            lead.append(pd.read_parquet(f).assign(chr=int(os.path.basename(f).split(".")[0][3:])))
    lead = pd.concat(lead, ignore_index=True)
    traits = [c for c in lead.columns if "__" in c]
    lead = lead.drop_duplicates(["phenotype_id", "variant_id"])
    d = m.merge(lead[["phenotype_id", "variant_id", "chr", "pos38", "af", "start_distance"] + traits],
                on=["phenotype_id", "variant_id"], how="inner")
    print(mod, len(m), "features;", len(d), "with GWAS z at the top variant", flush=True)
    maf = np.minimum(d.af, 1 - d.af)
    d["stratum"] = np.digitize(maf, [0.1, 0.2, 0.3]) * 4 + np.digitize(d.start_distance.abs(), [10_000, 50_000, 200_000])
    d = d.sort_values(["chr", "pos38"]).reset_index(drop=True)
    d["block"] = np.arange(len(d)) * NBLOCK // len(d)
    LF = d[[f"lfsr_{c}" for c in cts]].to_numpy()
    tested = np.isfinite(LF)
    sig = (LF < 0.05) & tested
    null = np.all(~tested | (LF > 0.5), axis=1)
    frac_sig = sig.sum(1) / np.maximum(tested.sum(1), 1)
    sets = {"anyCT": sig.any(1)}
    for i, ct in enumerate(cts):
        other_sig = (sig.sum(1) - sig[:, i]) / np.maximum(tested.sum(1) - tested[:, i], 1)
        sets[f"{ct}"] = sig[:, i]
        sets[f"{ct}:specific"] = sig[:, i] & (other_sig <= 0.2)
    sets["shared"] = frac_sig >= 0.8
    for name, s in sets.items():
        keep = s | null
        dd = d[keep].reset_index(drop=True)
        for r in test(dd, s[keep], traits):
            rows.append(dict(modality=mod, set=name, **r))
        print(mod, name, int(s.sum()), flush=True)

res = pd.DataFrame(rows).merge(man, on="trait_id", how="left")
p = res.p.to_numpy(float)
ok = np.isfinite(p)
o = np.argsort(p[ok])
q = np.minimum.accumulate((p[ok][o] * ok.sum() / np.arange(1, ok.sum() + 1))[::-1])[::-1]
fdr = np.full(len(p), np.nan)
tmp = np.empty(ok.sum())
tmp[o] = np.minimum(q, 1)
fdr[ok] = tmp
res["FDR"] = fdr
res = res.sort_values("p")
res.to_csv(f"{R}/mash_trait_enrichment.tsv", sep="\t", index=False, float_format="%.4g")
print(res.groupby(["modality", "set"]).agg(n_traits=("p", "size"), med_ratio=("chi2_ratio", "median"),
                                           frac_gt1=("chi2_ratio", lambda x: (x > 1).mean()),
                                           nom05=("p", lambda x: (x < 0.05).sum()),
                                           fdr10=("FDR", lambda x: (x < 0.1).sum())).round(3).to_string())
top = res.head(30).copy()
top["label"] = top.label.str[:35]
print(top[["modality", "set", "source", "label", "n_set", "chi2_ratio", "p", "FDR"]].to_string())
