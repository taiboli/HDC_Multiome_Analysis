"""
Harmonize GWAS summary statistics (FinnGen R10, Pan-UKB, BiobankJapan, PAGE; files collected for the EB
analyses under EBQTL/GWAS/) onto the HDC genotype variant universe (06a_variant_universe.R).

  python 06b_harmonize_gwas.py manifest      -> GWAS/gwas_manifest.tsv (one row per trait x ancestry analysis)
  python 06b_harmonize_gwas.py run <row>     -> GWAS/harmonized/<trait_id>.parquet  (chr, vidx, beta, se)
                                                GWAS/harmonized/qc/<trait_id>.tsv     (match rate, lambda, sdY)

Matching: chr + position (pos38 for FinnGen/PAGE, pos19 for Pan-UKB/BBJ) + alleles (either orientation);
beta is re-signed to the HDC ALT allele (= tensorQTL effect allele). Unmatched / ambiguous-allele rows dropped.
sdY (needed by coloc for quantitative traits) is estimated as in coloc:::sdY.est from varbeta, MAF and N.
"""
import glob
import os
import sys

import numpy as np
import pandas as pd

G = os.environ["GWAS_DIR"] + "/"
OUT = os.environ["HDC_DIR"] + "/GWAS/"
MANIFEST = f"{OUT}/gwas_manifest.tsv"


def build_manifest():
    rows = []
    # FinnGen R10 (hg38, Finnish); drop redundant alternative T1D / T2D definitions
    fg = pd.read_csv(f"{G}/FinnGen/FinnGen_summary_stats_R10_manifest.tsv", sep="\t").set_index("phenocode")
    for f in sorted(glob.glob(f"{G}/FinnGen/GWAS_*_FinnGen_*.gz")):
        code = os.path.basename(f).split("_FinnGen_")[1][:-3]
        if code in ("T1D_EARLY", "T1D_STRICT2", "T1D_WIDE", "T2D_WIDE"):
            continue
        m = fg.loc[code]
        cc = m.num_controls > 0
        rows.append(dict(trait_id=f"FinnGen__{code}", source="FinnGen", ancestry="FIN", label=m.phenotype,
                         file=f, build="hg38", format="finngen", pop="", type="cc" if cc else "quant",
                         N=int(m.num_cases + m.num_controls), ncase=int(m.num_cases) if cc else 0))
    # Pan-UKB (hg19): multi-ancestry high-quality meta, AFR (ancestry-matched to YRI donors), EUR
    man = pd.read_csv(f"{G}/Pan-UK Biobank phenotype manifest - phenotype_manifest.csv", low_memory=False)
    for f in sorted(glob.glob(f"{G}/PAGE/UKB/*.tsv.bgz")):
        b = os.path.basename(f)[:-8]
        ttype, code = b.split("-")[0], b.split("-")[1]
        hit = man[(man.trait_type == ttype) & (man.phenocode.astype(str) == code)]
        m = hit.iloc[0] if len(hit) else None
        label = code if m is None else (m.description if pd.isna(m.coding_description)
                                        else f"{m.description}: {m.coding_description}")
        cc = ttype in ("categorical", "icd10", "phecode")
        for pop in ("meta_hq", "AFR", "EUR"):
            if m is None:
                n, nc = np.nan, 0
            elif pop == "meta_hq":
                nc = m.n_cases_hq_cohort_both_sexes
                n = nc + (np.nansum([m[f"n_controls_{p}"] for p in ("AFR", "AMR", "CSA", "EAS", "EUR", "MID")]) if cc else 0)
            else:
                nc = m[f"n_cases_{pop}"]
                n = nc + (m[f"n_controls_{pop}"] if cc else 0)
            rows.append(dict(trait_id=f"UKB_{pop}__{code}", source=f"UKB_{pop}",
                             ancestry={"meta_hq": "multi", "AFR": "AFR", "EUR": "EUR"}[pop], label=label, file=f,
                             build="hg19", format="panukb", pop=pop, type="cc" if cc else "quant",
                             N=n, ncase=nc if cc else 0))
    # BiobankJapan (hg19, East Asian)
    bbj = [("Height", "2019_BBJ_Height_autosomes_BOLT.txt.gz", "bbj_height", "quant", 191787),
           ("BMI", "All_2017_BMI_BBJ_autosome.txt.gz", "bbj_std", "quant", 158284),
           ("CPD", "All_2018_CPD_BBJ_autosome_Pcorrected.txt.gz", "bbj_cpd", "quant", np.nan),
           ("T2D", "BBJ_BetaBased1.MAF_001.AtLeast2studies.AllChr_v2.txt.gz", "bbj_std", "cc", 191764),
           ("CAD", "CAD.auto.rsq07.mac10.txt.gz", "bbj_saige", "cc", np.nan),
           ("CeAn", "CeAn.auto.rsq07.mac10.txt.gz", "bbj_saige", "cc", np.nan),
           ("Epilepsy", "Epilepsy.auto.rsq07.mac10.txt.gz", "bbj_saige", "cc", np.nan),
           ("PrCa", "PrCa.auto.rsq07.mac10.txt.gz", "bbj_saige", "cc", np.nan),
           ("Asthma", "Asthma.chr*.auto_sex_stratified.txt.gz", "bbj_sexstrat", "cc", np.nan),
           ("Glaucoma", "Glaucoma.chr*.auto_sex_stratified.txt.gz", "bbj_sexstrat", "cc", np.nan)]
    for t in ["CRP", "DBP", "eGFR", "HbA1c", "HDL-C", "LDL-C", "MCHC", "Plt", "SBP", "TC", "TG", "WBC"]:
        bbj.append((t, f"BBJ.{t}.autosome.txt.gz", "bbj_std", "quant", np.nan))  # N column in file
    for t, f, fmt, ty, n in bbj:
        rows.append(dict(trait_id=f"BBJ__{t}", source="BBJ", ancestry="EAS", label=t, file=f"{G}/BBJ/{f}",
                         build="hg19", format=fmt, pop="", type=ty, N=n, ncase=0))
    # PAGE (Wojcik 2019; hg38 harmonized GWAS Catalog files; multi-ethnic, ~35% African American)
    pm = pd.read_csv(f"{G}/PAGE/PAGE_UKB_EFO_mapping_type.csv")
    for _, m in pm.iterrows():
        f = glob.glob(f"{G}/PAGE/31217584-{m.Accession}-*.h.tsv.gz")[0]
        rows.append(dict(trait_id=f"PAGE__{m.Accession}", source="PAGE", ancestry="multi(AFR-rich)",
                         label=m.Trait, file=f, build="hg38", format="page", pop="",
                         type="cc" if m.Case > 0 else "quant", N=m.Sample_Size, ncase=m.Case))
    df = pd.DataFrame(rows)
    os.makedirs(f"{OUT}/harmonized/qc", exist_ok=True)
    df.to_csv(MANIFEST, sep="\t", index=False)
    print(df.groupby("source").size())


class PopulationMissing(Exception):
    pass


def read_gwas(m):
    """Return DataFrame chr, pos, a_other, a_eff, beta, se, af (effect allele freq, may be NaN), N (may be NaN)."""
    fmt = m.format
    if fmt == "finngen":
        ren = {"#chrom": "chr", "pos": "pos", "ref": "a_other", "alt": "a_eff", "beta": "beta", "sebeta": "se",
               "af_alt": "af"}
        d = pd.read_csv(m.file, sep="\t", usecols=list(ren)).rename(columns=ren)
    elif fmt == "panukb":
        p = m["pop"]
        hdr = pd.read_csv(m.file, sep="\t", nrows=0, compression="gzip").columns
        if p == "meta_hq" and "beta_meta_hq" not in hdr:  # no high-quality meta for this trait -> full meta
            p = "meta"
        if f"beta_{p}" not in hdr:
            raise PopulationMissing(p)
        af = f"af_{p}" if f"af_{p}" in hdr else f"af_controls_{p}"
        ren = {"chr": "chr", "pos": "pos", "ref": "a_other", "alt": "a_eff", f"beta_{p}": "beta", f"se_{p}": "se",
               af: "af"}
        cols = list(ren) + ([f"low_confidence_{p}"] if not p.startswith("meta") else [])
        d = pd.read_csv(m.file, sep="\t", usecols=cols, compression="gzip", dtype={"chr": str}, low_memory=False)
        if not p.startswith("meta"):
            d = d[d[f"low_confidence_{p}"].astype(str).str.lower() == "false"].drop(columns=f"low_confidence_{p}")
        d = d.rename(columns=ren)
    elif fmt == "page":
        ren = {"hm_chrom": "chr", "hm_pos": "pos", "hm_other_allele": "a_other", "hm_effect_allele": "a_eff",
               "hm_beta": "beta", "standard_error": "se", "hm_effect_allele_frequency": "af"}
        d = pd.read_csv(m.file, sep="\t", usecols=list(ren), dtype={"hm_chrom": str}, low_memory=False).rename(columns=ren)
    elif fmt == "bbj_std":  # SNP CHR POS REF ALT Frq ... BETA SE P [N]; effect = ALT
        d = pd.read_csv(m.file, sep=r"\s+", usecols=lambda c: c in ("CHR", "POS", "REF", "ALT", "Frq", "BETA", "SE", "N"))
        d = d.rename(columns={"CHR": "chr", "POS": "pos", "REF": "a_other", "ALT": "a_eff", "BETA": "beta",
                              "SE": "se", "Frq": "af"})
    elif fmt == "bbj_height":
        ren = {"CHR": "chr", "POS": "pos", "REF": "a_other", "ALT": "a_eff", "ALT_freq": "af", "BETA": "beta", "SE": "se"}
        d = pd.read_csv(m.file, sep=" ", usecols=list(ren)).rename(columns=ren)
    elif fmt == "bbj_cpd":  # A1 A2 A1Frq: effect allele taken as A1 (sign used only for MR direction)
        d = pd.read_csv(m.file, sep="\t", usecols=["CHR", "POS", "A1", "A2", "A1Frq", "BETA", "SE"])
        d = d.rename(columns={"CHR": "chr", "POS": "pos", "A1": "a_eff", "A2": "a_other", "A1Frq": "af",
                              "BETA": "beta", "SE": "se"})
    elif fmt == "bbj_saige":  # Allele2 = ALT = effect
        d = pd.read_csv(m.file, sep=" ", usecols=["CHR", "POS", "Allele1", "Allele2", "AF_Allele2", "N", "BETA", "SE"])
        d = d.rename(columns={"CHR": "chr", "POS": "pos", "Allele1": "a_other", "Allele2": "a_eff",
                              "AF_Allele2": "af", "BETA": "beta", "SE": "se"})
    elif fmt == "bbj_sexstrat":  # sex-combined fixed-effect meta of male / female GWAS
        ren = {"CHR": "chr", "POS": "pos", "Allele1": "a_other", "Allele2": "a_eff", "AF.Controls.x": "af",
               "metabeta": "beta", "metase": "se"}
        d = pd.concat([pd.read_csv(f, sep="\t", usecols=list(ren)) for f in sorted(glob.glob(m.file))]).rename(columns=ren)
    else:
        raise ValueError(fmt)
    d["chr"] = pd.to_numeric(d["chr"].astype(str).str.replace("chr", ""), errors="coerce")
    d = d.dropna(subset=["chr", "pos", "beta", "se"])
    d = d[(d.se > 0) & d.chr.between(1, 22)]
    d["chr"] = d.chr.astype(np.int64)
    d["pos"] = d.pos.astype(np.int64)
    return d


def sdy_est(varbeta, maf, n):
    ok = np.isfinite(varbeta) & np.isfinite(maf) & (maf > 0.01) & np.isfinite(n)
    oneover, nvx = 1 / varbeta[ok], 2 * n[ok] * maf[ok] * (1 - maf[ok])
    cf = np.sum(oneover * nvx) / np.sum(oneover ** 2)  # lm(nvx ~ oneover - 1)
    return np.sqrt(cf) if cf > 0 else np.nan


def harmonize(row):
    m = pd.read_csv(MANIFEST, sep="\t").iloc[row]
    out = f"{OUT}/harmonized/{m.trait_id}.parquet"
    try:
        d = read_gwas(m)
    except PopulationMissing as e:  # e.g. no AFR analysis for this Pan-UKB trait: record and skip
        pd.DataFrame([dict(trait_id=m.trait_id, n_input=0, n_matched=0, note=f"no {e} analysis")]).to_csv(
            f"{OUT}/harmonized/qc/{m.trait_id}.tsv", sep="\t", index=False)
        print(m.trait_id, "skipped: no", e)
        return
    n_in = len(d)
    pos_col = "pos38" if m.build == "hg38" else "pos19"
    u = pd.read_parquet(f"{OUT}/variant_universe.parquet", columns=["vidx", "chr", pos_col, "ref", "alt"])
    u = u.dropna(subset=[pos_col])
    u["key"] = u.chr.astype(np.int64) * 1_000_000_000 + u[pos_col].astype(np.int64)
    d["key"] = d.chr * 1_000_000_000 + d.pos
    x = d.merge(u[["key", "vidx", "ref", "alt"]], on="key")
    a_eff, a_oth = x.a_eff.astype(str).str.upper(), x.a_other.astype(str).str.upper()
    same = (a_eff == x.alt) & (a_oth == x.ref)
    flip = (a_eff == x.ref) & (a_oth == x.alt)
    x = x[same | flip].copy()
    sign = np.where(same[same | flip], 1.0, -1.0)
    x["beta"] = x.beta * sign
    x["af"] = np.where(sign > 0, x.af, 1 - x.af)
    x["z2"] = (x.beta / x.se) ** 2
    x = x.sort_values("z2", ascending=False).drop_duplicates("vidx").sort_values("vidx")
    x["chr"] = x.key // 1_000_000_000
    res = pd.DataFrame({"chr": x.chr.astype(np.int8), "vidx": x.vidx.astype(np.int32),
                        "beta": x.beta.astype(np.float32), "se": x.se.astype(np.float32)})
    res.to_parquet(out + ".tmp", index=False, row_group_size=2_000_000)
    os.replace(out + ".tmp", out)

    from scipy.stats import chi2
    lam = np.median(x.z2) / chi2.ppf(0.5, 1)
    sdy = np.nan
    if m.type == "quant":
        n = x["N"].to_numpy(float) if "N" in x and x["N"].notna().any() else np.full(len(x), m.N, float)
        maf = np.minimum(x.af, 1 - x.af).to_numpy(float)
        sdy = sdy_est((x.se ** 2).to_numpy(float), maf, n)
    qc = dict(trait_id=m.trait_id, n_input=n_in, n_matched=len(res), frac_universe=len(res) / len(u),
              n_flipped=int((sign < 0).sum()), lambda_gc=lam, n_gws=int((x.z2 > 29.716).sum()), sdY=sdy)
    pd.DataFrame([qc]).to_csv(f"{OUT}/harmonized/qc/{m.trait_id}.tsv", sep="\t", index=False)
    print(qc, flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "manifest":
        build_manifest()
    else:
        harmonize(int(sys.argv[2]))
