#!/usr/bin/env python3
# =============================================================================
# 06_diagnose_residual.py   (generalized: any deleted copy P / M1 / M2)
#   For each deletion experiment present in the manifest, ask why any
#   deleted-copy-allele signal survives in that deletion's RNA-seq. Ground
#   truth is the deletion line's WGS: if copy D is truly gone, D's private
#   allele must have ~0 reads in the dD genomic DNA.
#
#   Auto-detects every deleted_copy in $OUT/rna/manifest.tsv (P|M1|M2) that has
#   RNA samples, and writes one diagnosis per copy.
#
#   Per D-class site (resolved_copy==D):
#     rna_del_frac     D-allele fraction in pooled D-deletion RNA  (~0 expected)
#     wgs_delline_frac D-allele fraction in dD WGS                 (~0 expected)
#     wgs_retainA/B    D-allele fraction in the two retaining lines (~0.5)
#     wgs_T21_frac     D-allele fraction in T21 WGS                (~1/3)
#   cause for rna_del_frac>RNA_THR: MISMAP_WGS (dD WGS not clean) | RNA_ONLY
#
#   Output: $OUT/fig/residual_diagnosis_<D>.tsv
# =============================================================================
import os, sys
import pysam
import pandas as pd
import numpy as np

OUT     = os.environ.get("OUT", os.path.expanduser("~/Desktop/DS/ASE_out"))
RNA_THR = float(os.environ.get("RES_RNA_THR", "0.02"))
WGS_THR = float(os.environ.get("RES_WGS_THR", "0.05"))
TRIO    = f"{OUT}/vcf/trio.chr21.vcf.gz"
T21QC   = f"{OUT}/vcf/T21.qc.vcf.gz"
DISOMY  = {"P": "dP", "M1": "dM1", "M2": "dM2"}     # copy -> deletion line

# ---- manifest ---------------------------------------------------------------
man = pd.read_csv(f"{OUT}/rna/manifest.tsv", sep="\t", comment="#", header=None,
                  names=["sid", "r1", "r2", "condition", "deleted_copy"])
man["deleted_copy"] = man["deleted_copy"].astype(str)

# ---- AD loaders (one pass each) ---------------------------------------------
def load_ad(path):
    d = {}
    if not os.path.exists(path):
        sys.stderr.write(f"[06] WARN: {path} missing\n"); return d
    for r in pysam.VariantFile(path):
        if r.alts is None or len(r.alts) != 1:
            continue
        rec = {}
        for s in r.samples:
            ad = r.samples[s].get("AD")
            if ad is not None and len(ad) >= 2:
                rec[s] = (int(ad[0]), int(ad[1]))
        d[(r.pos, r.ref, r.alts[0])] = rec
    return d

trio_ad = load_ad(TRIO)
t21_ad  = load_ad(T21QC)

def frac(adt, idx):
    if adt is None:
        return np.nan
    tot = adt[0] + adt[1]
    return adt[idx] / tot if tot > 0 else np.nan

# ---- sites ------------------------------------------------------------------
sites = pd.read_csv(f"{OUT}/snp/sites_with_gene.tsv", sep="\t")

def pooled_rna(sids):
    """pos -> [refSum, altSum] over the given RNA samples."""
    acc = {}
    for sid in sids:
        p = f"{OUT}/ase/{sid}.ase.table"
        if not os.path.exists(p):
            sys.stderr.write(f"[06] WARN: {p} missing\n"); continue
        t = pd.read_csv(p, sep="\t")
        for _, row in t.iterrows():
            a = acc.setdefault(int(row["position"]), [0, 0])
            a[0] += int(row["refCount"]); a[1] += int(row["altCount"])
    return acc

# ---- per deleted copy -------------------------------------------------------
targets = [c for c in ("P", "M1", "M2")
           if (man["deleted_copy"] == c).any()]
if not targets:
    sys.exit("[06] no deletion conditions found in manifest (deleted_copy in P/M1/M2).")

for D in targets:
    del_line = DISOMY[D]
    retain   = [v for k, v in DISOMY.items() if k != D]   # two retaining lines
    del_sids = man.loc[man["deleted_copy"] == D, "sid"].tolist()
    rna = pooled_rna(del_sids)

    sub = sites[sites["resolved_copy"] == D].copy()
    sub["idx"] = sub["odd_is_alt"].astype(int)            # D-allele index in AD
    rows = []
    for _, r in sub.iterrows():
        pos = int(r["pos"]); idx = int(r["idx"]); k = (pos, r["ref"], r["alt"])
        rc = rna.get(pos)
        if rc is None:
            rfrac = reads = tot = np.nan
        else:
            tot = rc[0] + rc[1]; reads = rc[idx]
            rfrac = reads / tot if tot > 0 else np.nan
        tad = trio_ad.get(k, {}); qad = t21_ad.get(k, {})
        rows.append({
            "copy": D, "gene": r["gene"], "pos": pos, "ref": r["ref"], "alt": r["alt"],
            "rna_del_frac": rfrac, "rna_reads": reads, "rna_total": tot,
            "wgs_delline_frac": frac(tad.get(del_line), idx),
            "wgs_delline_reads": (tad.get(del_line) or (np.nan, np.nan))[idx],
            "wgs_retainA_frac": frac(tad.get(retain[0]), idx),
            "wgs_retainB_frac": frac(tad.get(retain[1]), idx),
            "wgs_T21_frac": frac(next(iter(qad.values()), None), idx),
        })
    df = pd.DataFrame(rows)

    def classify(row):
        if pd.isna(row["rna_del_frac"]) or row["rna_del_frac"] <= RNA_THR:
            return "ok/low_floor"
        if not pd.isna(row["wgs_delline_frac"]) and row["wgs_delline_frac"] > WGS_THR:
            return "MISMAP_WGS"
        return "RNA_ONLY"
    df["cause"] = df.apply(classify, axis=1)
    df = df.sort_values("rna_del_frac", ascending=False)
    out = f"{OUT}/fig/residual_diagnosis_{D}.tsv"
    df.to_csv(out, sep="\t", index=False)

    flagged = df[df["rna_del_frac"] > RNA_THR]
    print(f"\n[06] === deleted copy {D}  (line {del_line}; RNA n={len(del_sids)}) ===")
    print(f"[06] {D}-informative sites examined: {len(df)}")
    print(f"[06] sites with del-RNA {D}-fraction > {RNA_THR}: {len(flagged)}")
    for c, k in flagged["cause"].value_counts().items():
        print(f"        {c}: {k}")
    genes = sorted(flagged["gene"].dropna().astype(str)
                   [flagged["gene"].astype(str) != "NA"].unique().tolist())
    print(f"[06] genes affected: {genes}")
    print(f"[06] median del-RNA {D}-fraction (all sites): "
          f"{df['rna_del_frac'].median():.4f}")
    print(f"[06] -> {out}")
