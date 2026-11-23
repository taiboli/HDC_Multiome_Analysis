"""
LDSC .sumstats.gz for every harmonized GWAS (06b), restricted to HapMap3 regression SNPs (the union of the
1000G Phase 3 EUR and EAS hm3_noMHC weight SNPs), strand-ambiguous SNPs removed.
Columns: SNP A1 A2 Z N with A1 = HDC ALT allele (the allele the harmonized beta refers to), Z = beta / se,
N = manifest N (total sample size; observed-scale h2 for binary traits, which does not affect enrichment).
Missing manifest N: median of the source file's N column if it has one (BBJ quantitative traits, CAD, CeAn,
Epilepsy, PrCa); UKB_meta_hq -> N of the UKB_EUR GWAS of the same phenotype; otherwise (BBJ Asthma, Glaucoma,
CPD) a placeholder N = 100,000 flagged in sumstats_qc.tsv (N_source) -- h2 is then not interpretable, but
S-LDSC enrichment and tau z-scores do not depend on N.

Output: GWAS/ldsc/sumstats/<trait_id>.sumstats.gz
"""
import glob
import os

import numpy as np
import pandas as pd

BASE = os.environ["HDC_DIR"] + "/"
REF = os.environ["LDSC_REF_DIR"] + "/"
OUT = f"{BASE}/GWAS/ldsc/sumstats/"
os.makedirs(OUT, exist_ok=True)

snps = set()
for pat in ["1000G_Phase3_weights_hm3_no_MHC/*/weights.hm3_noMHC.*.l2.ldscore.gz",
            "1000G_Phase3_EAS_weights_hm3_no_MHC/*/weights.EAS.hm3_noMHC.*.l2.ldscore.gz"]:
    for f in glob.glob(REF + pat):
        snps |= set(pd.read_csv(f, sep="\t", usecols=["SNP"]).SNP)
print("hm3 regression SNPs:", len(snps), flush=True)

u = pd.read_parquet(f"{BASE}/GWAS/variant_universe.parquet", columns=["vidx", "id", "ref", "alt"])
u = u[u.id.isin(snps) & u.ref.str.len().eq(1) & u.alt.str.len().eq(1)]
amb = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}
u = u[[(r, a) not in amb for r, a in zip(u.ref, u.alt)]].set_index("vidx")
print("HDC variants usable:", len(u), flush=True)

man = pd.read_csv(f"{BASE}/GWAS/gwas_manifest.tsv", sep="\t")
man["N_source"] = np.where(man.N.notna(), "manifest", "")
e = man[man.trait_id.str.startswith("UKB_EUR__")]
eurN = dict(zip(e.trait_id.str.replace("UKB_EUR__", "UKB_meta_hq__", regex=False), e.N))
for i, m in man[man.N.isna()].iterrows():
    n, src = np.nan, ""
    if m.trait_id.startswith("UKB_meta_hq__") and pd.notna(eurN.get(m.trait_id)) and m.trait_id in eurN:
        n, src = eurN[m.trait_id], "UKB_EUR_same_pheno"
    elif not glob.glob(m.file.replace("*", "1")) == []:
        f0 = glob.glob(m.file.replace("*", "1"))[0]
        hdr = pd.read_csv(f0, sep=r"\s+", nrows=0, compression="gzip").columns
        if "N" in hdr:
            n, src = pd.read_csv(f0, sep=r"\s+", usecols=["N"], compression="gzip").N.median(), "file_median"
    if pd.isna(n):
        n, src = 100_000, "placeholder"
    man.loc[i, ["N", "N_source"]] = [n, src]
rows = []
for _, m in man.iterrows():
    f = f"{BASE}/GWAS/harmonized/{m.trait_id}.parquet"
    if not os.path.exists(f) or os.path.exists(f"{OUT}/{m.trait_id}.sumstats.gz"):
        continue
    g = pd.read_parquet(f, columns=["vidx", "beta", "se"])
    g = g[g.vidx.isin(u.index) & (g.se > 0) & np.isfinite(g.beta)]
    d = u.loc[g.vidx].reset_index(drop=True)
    s = pd.DataFrame({"SNP": d.id, "A1": d.alt, "A2": d.ref, "Z": (g.beta / g.se).to_numpy(), "N": int(m.N)})
    s = s.drop_duplicates("SNP", keep=False)
    s.to_csv(f"{OUT}/{m.trait_id}.sumstats.gz", sep="\t", index=False, float_format="%.5g")
    rows.append(dict(trait_id=m.trait_id, n_snp=len(s), N=int(m.N), N_source=m.N_source,
                     mean_chi2=float(np.mean(s.Z ** 2))))
    print(m.trait_id, len(s), flush=True)
qc = f"{OUT}/sumstats_qc.tsv"
old = pd.read_csv(qc, sep="\t") if os.path.exists(qc) else pd.DataFrame()
new = pd.DataFrame(rows)
if len(old) and "N_source" not in old:
    old["N_source"] = "manifest"
pd.concat([old, new]).drop_duplicates("trait_id", keep="last").to_csv(qc, sep="\t", index=False)
print("done")
