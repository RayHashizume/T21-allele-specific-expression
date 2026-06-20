#!/usr/bin/env python3
# =============================================================================
# 08_statistics.py   (publication statistics; PRIMARY unit = biological replicate)
#   n=3 replicates per group are treated as the unit of replication (NOT pooled
#   to one number and NOT pseudoreplicated at the read level). For every
#   biological replicate we compute three chr21-wide summaries, then compare
#   baseline vs each condition as an UNPAIRED two-sample problem (baseline and
#   the deletion line are independent samples):
#
#     PRIMARY (replicate-level, n per group):
#       frac_C   chr21-wide allele fraction of copy C (sites pooled within the
#                replicate, then P:M1:M2 normalized to sum 1)
#       total    chr21 expression dosage = sum(chr21 gene counts)/non-chr21
#                mapped reads x 1e6
#       abs_C    absolute expression of copy C = total x frac_C
#     Tests per class/metric: Welch t-test (primary) + Mann-Whitney U (n=3 is
#     low-powered for MWU: min two-sided p ~0.1, so report both).
#
#     SECONDARY (gene as unit): paired Wilcoxon signed-rank across genes of the
#     read-pooled per-gene fraction (baseline vs condition), per class.
#
#   Inputs: manifest, sites_with_gene.tsv, ase/*.ase.table, expr/gene_counts.txt,
#           expr/libsizes.tsv  (run 04 and 07 first).
#   Outputs ($OUT/stats): replicate_metrics.tsv, stats_summary.tsv,
#           fig_stats_fraction_<cond>.pdf, fig_stats_absolute_<cond>.pdf
# =============================================================================
import os, sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT    = os.environ.get("OUT", os.path.expanduser("~/Desktop/DS/ASE_out"))
CHR21  = os.environ.get("CHR21", "chr21")
MINTOT = int(os.environ.get("ASE_MIN_TOTAL", "20"))     # gene-level (secondary)
COPIES = ["P", "M1", "M2"]
PAL    = {"P": "#3B6FB6", "M1": "#E08214", "M2": "#7B3294"}
ST     = f"{OUT}/stats"; os.makedirs(ST, exist_ok=True)

man = pd.read_csv(f"{OUT}/rna/manifest.tsv", sep="\t", comment="#", header=None,
                  names=["sid", "r1", "r2", "condition", "deleted_copy"])
man["deleted_copy"] = man["deleted_copy"].astype(str)
man["group"] = np.where(man["deleted_copy"] == "none", "baseline", man["condition"])
sid2grp = dict(zip(man["sid"], man["group"]))
comp = [g for g in man["group"].unique() if g != "baseline"]

# ---- sites + per-sample ASE -------------------------------------------------
sites = pd.read_csv(f"{OUT}/snp/sites_with_gene.tsv", sep="\t")
sites_g = sites[sites["gene"].notna() & (sites["gene"] != "NA")].copy()
sites_g["gene"] = sites_g["gene"].str.split(",")
sites_g = sites_g.explode("gene", ignore_index=True)

def load_ase(sid):
    t = pd.read_csv(f"{OUT}/ase/{sid}.ase.table", sep="\t").rename(
        columns={"position": "pos"})
    # use ALL informative sites (gene or not) for the chr21-wide replicate frac
    m = sites[["pos", "resolved_copy", "odd_is_alt"]].merge(
        t[["pos", "refCount", "altCount", "totalCount"]], on="pos")
    m["odd"] = np.where(m["odd_is_alt"] == 1, m["altCount"], m["refCount"])
    m["tot"] = m["totalCount"]; m["sid"] = sid
    return m[m["tot"] > 0]

ase = pd.concat([load_ase(s) for s in man["sid"]], ignore_index=True)

# ---- chr21 total dosage per replicate (from featureCounts + libsizes) -------
fc = pd.read_csv(f"{OUT}/expr/gene_counts.txt", sep="\t", comment="#")
meta = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
ccols = [c for c in fc.columns if c not in meta]
ren = {c: os.path.basename(c).replace(".bam", "") for c in ccols}
fc = fc.rename(columns=ren)
chr21_total = fc[list(ren.values())].sum()                 # total chr21 counts/sample
libs = pd.read_csv(f"{OUT}/expr/libsizes.tsv", sep="\t").set_index("sid")["nonchr21_mapped"]
dosage = {s: chr21_total[s] / libs[s] * 1e6 for s in man["sid"]}

# ---- PRIMARY: per-replicate metrics -----------------------------------------
rows = []
for s in man["sid"]:
    sub = ase[ase["sid"] == s]
    fr = {}
    for c in COPIES:
        d = sub[sub["resolved_copy"] == c]
        fr[c] = d["odd"].sum() / d["tot"].sum() if d["tot"].sum() > 0 else np.nan
    ssum = sum(v for v in fr.values() if not np.isnan(v))
    rec = {"sid": s, "group": sid2grp[s], "total_chr21": dosage[s]}
    for c in COPIES:
        rec[f"frac_{c}"] = fr[c] / ssum if ssum > 0 else np.nan       # normalized
        rec[f"abs_{c}"] = dosage[s] * (fr[c] / ssum) if ssum > 0 else np.nan
    rows.append(rec)
rep = pd.DataFrame(rows)
rep.to_csv(f"{ST}/replicate_metrics.tsv", sep="\t", index=False)

# ---- gene-level pooled fractions (SECONDARY) --------------------------------
genefr = (sites_g[["pos", "resolved_copy", "odd_is_alt", "gene"]]
          .merge(pd.concat([load_ase(s).assign(group=sid2grp[s])
                            for s in man["sid"]]), on=["pos", "resolved_copy"]))
gp = (genefr.groupby(["group", "gene", "resolved_copy"], as_index=False)
            .agg(odd=("odd", "sum"), tot=("tot", "sum")))
gp = gp[gp["tot"] >= MINTOT]
gp["frac"] = gp["odd"] / gp["tot"]

# ---- tests ------------------------------------------------------------------
def two_sample(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    out = {"baseline_mean": np.mean(a), "baseline_sd": np.std(a, ddof=1),
           "cond_mean": np.mean(b), "cond_sd": np.std(b, ddof=1),
           "n_base": len(a), "n_cond": len(b)}
    try:
        out["welch_p"] = stats.ttest_ind(a, b, equal_var=False).pvalue
    except Exception:
        out["welch_p"] = np.nan
    try:
        out["mwu_p"] = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
    except Exception:
        out["mwu_p"] = np.nan
    return out

results = []
for grp in comp:
    b = rep[rep["group"] == "baseline"]; d = rep[rep["group"] == grp]
    # total dosage
    r = two_sample(b["total_chr21"], d["total_chr21"])
    r.update({"condition": grp, "metric": "total_chr21_dosage", "copy": "-",
              "ratio_cond_base": r["cond_mean"] / r["baseline_mean"]})
    results.append(r)
    for c in COPIES:
        for metric, col in (("fraction", f"frac_{c}"), ("absolute", f"abs_{c}")):
            r = two_sample(b[col], d[col])
            r.update({"condition": grp, "metric": metric, "copy": c,
                      "ratio_cond_base": (r["cond_mean"] / r["baseline_mean"]
                                          if r["baseline_mean"] else np.nan)})
            results.append(r)
        # SECONDARY gene-level paired Wilcoxon (genes as unit)
        wb = gp[(gp["group"] == "baseline") & (gp["resolved_copy"] == c)] \
            .set_index("gene")["frac"]
        wd = gp[(gp["group"] == grp) & (gp["resolved_copy"] == c)] \
            .set_index("gene")["frac"]
        common = wb.index.intersection(wd.index)
        wp = (stats.wilcoxon(wb[common], wd[common]).pvalue
              if len(common) >= 6 else np.nan)
        results.append({"condition": grp, "metric": "gene_wilcoxon", "copy": c,
                        "baseline_mean": wb[common].median(),
                        "cond_mean": wd[common].median(),
                        "n_base": len(common), "n_cond": len(common),
                        "welch_p": np.nan, "mwu_p": np.nan, "wilcoxon_p": wp,
                        "ratio_cond_base": np.nan})
res = pd.DataFrame(results)
res.to_csv(f"{ST}/stats_summary.tsv", sep="\t", index=False)

# ---- figures: per-replicate points + mean/SD, p annotated -------------------
def metric_fig(grp, metric, ylabel, fname):
    col = {"fraction": "frac_{}", "absolute": "abs_{}"}[metric]
    fig, axes = plt.subplots(1, 3, figsize=(9, 4), sharey=(metric == "fraction"))
    for ax, c in zip(axes, COPIES):
        for j, g in enumerate(["baseline", grp]):
            vals = rep[rep["group"] == g][col.format(c)].dropna().values
            ax.scatter(np.full(len(vals), j + 1) + np.random.normal(0, 0.05, len(vals)),
                       vals, s=28, color=PAL[c], zorder=3)
            ax.errorbar(j + 1, np.mean(vals), yerr=np.std(vals, ddof=1),
                        fmt="_", color="black", capsize=6, ms=20, zorder=2)
        rr = res[(res["condition"] == grp) & (res["metric"] == metric) &
                 (res["copy"] == c)]
        p = rr["welch_p"].iat[0] if len(rr) else np.nan
        ax.set_title(f"{c}-allele  (Welch p={p:.2g})")
        ax.set_xticks([1, 2]); ax.set_xticklabels(["baseline", grp])
        if metric == "fraction":
            for y in (1/3, 1/2, 0):
                ax.axhline(y, ls="--", lw=0.6, color="grey")
    axes[0].set_ylabel(ylabel)
    fig.suptitle(f"{grp}: per-replicate {metric} (n=3 biological replicates)",
                 fontsize=11)
    fig.tight_layout(); fig.savefig(f"{ST}/{fname}"); plt.close(fig)

for grp in comp:
    metric_fig(grp, "fraction", "chr21-wide allele fraction",
               f"fig_stats_fraction_{grp}.pdf")
    metric_fig(grp, "absolute", "absolute expr (per 1e6 non-chr21 reads)",
               f"fig_stats_absolute_{grp}.pdf")

# ---- console ----------------------------------------------------------------
pd.set_option("display.width", 160)
for grp in comp:
    print(f"\n[08] === {grp} vs baseline  (n=3 biological replicates) ===")
    sub = res[res["condition"] == grp]
    for _, r in sub.iterrows():
        tag = f"{r['metric']}/{r['copy']}"
        if r["metric"] == "gene_wilcoxon":
            print(f"   {tag:<22} gene-paired Wilcoxon p={r['wilcoxon_p']:.2g} "
                  f"(median {r['baseline_mean']:.3f} -> {r['cond_mean']:.3f}, "
                  f"n={r['n_base']} genes)")
        else:
            print(f"   {tag:<22} {r['baseline_mean']:.3f}+/-{r['baseline_sd']:.3f}"
                  f" -> {r['cond_mean']:.3f}+/-{r['cond_sd']:.3f}  "
                  f"ratio={r['ratio_cond_base']:.2f}  "
                  f"Welch p={r['welch_p']:.2g}  MWU p={r['mwu_p']:.2g}")
print(f"\n[08] -> {ST}/replicate_metrics.tsv ; stats_summary.tsv ; "
      f"fig_stats_(fraction|absolute)_<cond>.pdf")
