#!/usr/bin/env python3
# =============================================================================
# 05_estimate_and_plot.py  (condition-driven; supports KNOWN and UNKNOWN edits)
#   Reads $OUT/rna/manifest.tsv (sid,R1,R2,condition,deleted_copy) and, for each
#   non-baseline condition, quantifies per-copy (P/M1/M2) expression and
#   compares to baseline. Works whether or not the targeted copy is known:
#     deleted_copy = none      -> baseline (trisomy, ~1/3 each)
#     deleted_copy = P|M1|M2   -> known target (highlighted as "expected")
#     deleted_copy = unknown   -> blind: the affected copy is INFERRED from data
#
#   Outputs ($OUT/fig):
#     dosage_summary.tsv          per condition: chr21-wide x_P,x_M1,x_M2,
#                                 normalized, ratio-to-baseline, auto-CALL
#     class_shift_<cond>.tsv/.pdf baseline vs condition, per class
#     gene_allele_fractions.tsv   per gene x condition (full genes)
#     fig_stacked_full.pdf        per-gene 3-allele stacks across conditions
#
#   RNA-seq measures expressed allele DOSAGE; it cannot distinguish physical
#   deletion from silencing/inactivation (both remove the RNA). Use DNA /
#   epigenetic assays for the mechanism.
# =============================================================================
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT    = os.environ.get("OUT", os.path.expanduser("~/Desktop/DS/ASE_out"))
MINTOT = int(os.environ.get("ASE_MIN_TOTAL", "20"))
# fractional drop (relative to baseline) below which a copy is "called" reduced
CALL_DROP = float(os.environ.get("CALL_DROP", "0.40"))
COPIES = ["P", "M1", "M2"]
PAL    = {"P": "#3B6FB6", "M1": "#E08214", "M2": "#7B3294"}
os.makedirs(f"{OUT}/fig", exist_ok=True)

# ---- manifest + groups ------------------------------------------------------
man = pd.read_csv(f"{OUT}/rna/manifest.tsv", sep="\t", comment="#", header=None,
                  names=["sid", "r1", "r2", "condition", "deleted_copy"])
man["deleted_copy"] = man["deleted_copy"].astype(str)
# group: baseline if deleted_copy==none, else the condition label (arbitrary)
man["group"] = np.where(man["deleted_copy"] == "none", "baseline", man["condition"])
if "baseline" not in set(man["group"]):
    sys.exit("[05] no baseline samples (deleted_copy=none) in manifest.")
# expected target per group (P/M1/M2 if known, else 'unknown'/'none')
expected = (man.groupby("group")["deleted_copy"]
            .agg(lambda s: s.mode().iat[0]).to_dict())
comp_groups = [g for g in man["group"].unique() if g != "baseline"]

# ---- sites + ASE ------------------------------------------------------------
sites = pd.read_csv(f"{OUT}/snp/sites_with_gene.tsv", sep="\t")
sites = sites[sites["gene"].notna() & (sites["gene"] != "NA")].copy()
sites["gene"] = sites["gene"].str.split(",")
sites = sites.explode("gene", ignore_index=True)
sites = sites[["pos", "resolved_copy", "odd_is_alt", "gene"]]

def load_ase(sid, group):
    t = pd.read_csv(f"{OUT}/ase/{sid}.ase.table", sep="\t").rename(
        columns={"position": "pos"})
    m = sites.merge(t[["pos", "refCount", "altCount", "totalCount"]], on="pos")
    m["oddCount"] = np.where(m["odd_is_alt"] == 1, m["altCount"], m["refCount"])
    m["total"] = m["totalCount"]; m["group"] = group
    return m[m["total"] > 0]

ase = pd.concat([load_ase(r.sid, r.group) for r in man.itertuples()],
                ignore_index=True)

# per (group, gene, class) -- for per-gene / stacked figures
pooled = (ase.groupby(["group", "gene", "resolved_copy"], as_index=False)
             .agg(odd=("oddCount", "sum"), tot=("total", "sum")))
pooled["frac"] = np.where(pooled["tot"] >= MINTOT,
                          pooled["odd"] / pooled["tot"], np.nan)
gene_pos = sites.groupby("gene")["pos"].min().rename("pos")

# =============================================================================
#  chr21-wide DOSAGE SUMMARY (read-weighted over all expressed informative
#  sites of each class) -- the headline read-out for blind/unknown samples
# =============================================================================
wide = (ase.groupby(["group", "resolved_copy"], as_index=False)
           .agg(odd=("oddCount", "sum"), tot=("total", "sum")))
wide["frac"] = wide["odd"] / wide["tot"]
dose = wide.pivot(index="group", columns="resolved_copy", values="frac")
for c in COPIES:
    if c not in dose.columns:
        dose[c] = np.nan
dose = dose[COPIES]
# normalize the three measured fractions to sum 1 -> relative dosage P:M1:M2
norm = dose.div(dose.sum(axis=1), axis=0)
base = norm.loc["baseline"]

rows = []
for g in dose.index:
    r = {"condition": g, "expected": expected.get(g, "none")}
    for c in COPIES:
        r[f"x{c}"] = dose.loc[g, c]
        r[f"x{c}_norm"] = norm.loc[g, c]
        r[f"x{c}_vs_base"] = norm.loc[g, c] / base[c] if base[c] > 0 else np.nan
    if g == "baseline":
        r["call"] = "baseline"
    else:
        # a copy is "reduced" if its normalized dosage fell by >CALL_DROP vs base
        drops = {c: 1 - r[f"x{c}_vs_base"] for c in COPIES
                 if not np.isnan(r[f"x{c}_vs_base"])}
        reduced = [c for c, d in drops.items() if d > CALL_DROP]
        if not reduced:
            r["call"] = "no copy clearly reduced"
        else:
            parts = [f"{c} reduced to {r[f'x{c}_vs_base']*100:.0f}% of baseline"
                     for c in sorted(reduced, key=lambda c: r[f'x{c}_vs_base'])]
            r["call"] = "; ".join(parts)
    rows.append(r)
dosage = pd.DataFrame(rows)
dosage.to_csv(f"{OUT}/fig/dosage_summary.tsv", sep="\t", index=False)

print("[05] chr21-wide allele dosage (normalized P:M1:M2) and call:")
for _, r in dosage.iterrows():
    print(f"      {r['condition']:<16} "
          f"P={r['xP_norm']:.3f} M1={r['xM1_norm']:.3f} M2={r['xM2_norm']:.3f}"
          f"   expected={r['expected']:<8} -> {r['call']}")

# =============================================================================
#  per-condition class-shift figures (baseline vs condition)
# =============================================================================
def class_shift_figure(grp):
    exp = expected.get(grp, "none")
    w = (pooled[pooled["group"].isin(["baseline", grp])]
         .pivot_table(index=["gene", "resolved_copy"], columns="group",
                      values="frac").reset_index())
    for col in ("baseline", grp):
        if col not in w.columns:
            w[col] = np.nan
    w.to_csv(f"{OUT}/fig/class_shift_{grp}.tsv", sep="\t", index=False)
    # which class actually dropped most (data-driven)
    inferred = min(COPIES, key=lambda c: dose.loc[grp, c]
                   if not np.isnan(dose.loc[grp, c]) else 9)
    fig, axes = plt.subplots(1, 3, figsize=(9, 4), sharey=True)
    for ax, cls in zip(axes, COPIES):
        s = w[w["resolved_copy"] == cls]
        data = [s["baseline"].dropna().values, s[grp].dropna().values]
        for j, vals in enumerate(data):
            ax.scatter(np.random.normal(j + 1, 0.06, len(vals)), vals,
                       s=10, alpha=0.5, color=PAL[cls], zorder=3)
        if any(len(d) for d in data):
            ax.boxplot(data, positions=[1, 2], widths=0.5, showfliers=False)
        for y in (1/3, 1/2, 0):
            ax.axhline(y, ls="--", lw=0.6, color="grey")
        tag = ""
        if cls == exp and exp in COPIES:
            tag = "  <- expected"
        elif exp not in COPIES and cls == inferred:
            tag = "  <- inferred reduced"
        ax.set_title(f"{cls}-allele (n={s['baseline'].notna().sum()}){tag}")
        ax.set_xticks([1, 2]); ax.set_xticklabels(["baseline", grp])
        ax.set_ylim(-0.05, 1.05)
    axes[0].set_ylabel("expression fraction of the copy")
    ttl = (f"{grp}: expected target {exp}" if exp in COPIES
           else f"{grp}: blind -> inferred reduced copy = {inferred}")
    fig.suptitle(ttl, fontsize=11)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig/fig_class_shift_{grp}.pdf")
    plt.close(fig)

for g in comp_groups:
    class_shift_figure(g)

# =============================================================================
#  per-gene 3-allele stacked bars across baseline + conditions
# =============================================================================
def fractions(group):
    sub = pooled[(pooled["group"] == group) & (pooled["frac"].notna())]
    return {g: dict(zip(d["resolved_copy"], d["frac"]))
            for g, d in sub.groupby("gene")}

def complete_fracs(d, deleted):
    """Complete a per-gene {copy:frac} dict to all three copies.
    The deleted/affected copy must be DIRECTLY measured (never residual-filled,
    so its ~0 is real); at most ONE retained copy may be residual-completed."""
    missing = [c for c in COPIES if c not in d]
    if deleted in COPIES and deleted in missing:
        return None                      # deleted copy lacks a direct site -> drop
    if len(missing) == 0:
        v = dict(d)
    elif len(missing) == 1:
        v = dict(d); v[missing[0]] = max(0.0, 1 - sum(d.values()))
    else:
        return None                      # >1 copy missing -> not resolvable
    s = sum(v.values())
    return {c: v[c] / s for c in COPIES} if s > 0 else None

rows = []
for group in ["baseline"] + comp_groups:
    deleted = expected.get(group, "none")
    # blind (unknown) conditions: the affected copy is not known a priori, so
    # require all three direct to avoid completing the unknown affected copy.
    strict = deleted not in COPIES and deleted != "none"
    for gene, d in fractions(group).items():
        if strict:
            if len(d) < 3:
                continue
            s = sum(d.values()); v = {c: d[c] / s for c in COPIES}
        else:
            v = complete_fracs(d, deleted if deleted in COPIES else None)
            if v is None:
                continue
        rows.append({"gene": gene, "group": group,
                     "xP": v["P"], "xM1": v["M1"], "xM2": v["M2"]})
res = pd.DataFrame(rows)
if not res.empty:
    res = res.merge(gene_pos, on="gene", how="left")
res.to_csv(f"{OUT}/fig/gene_allele_fractions.tsv", sep="\t", index=False)

groups_all = ["baseline"] + comp_groups
keep = [g for g, sub in res.groupby("gene")
        if set(groups_all).issubset(set(sub["group"]))] if not res.empty else []
plot = res[res["gene"].isin(keep)]
if not plot.empty:
    order = plot.drop_duplicates("gene").sort_values("pos")["gene"].tolist()
    n = len(groups_all)
    fig, axes = plt.subplots(n, 1, figsize=(max(8, len(order) * 0.22), 2.2 * n),
                             sharex=True)
    if n == 1:
        axes = [axes]
    for ax, grp in zip(axes, groups_all):
        sub = plot[plot["group"] == grp].set_index("gene").reindex(order)
        x = np.arange(len(order)); bottom = np.zeros(len(order))
        for c in COPIES:
            ax.bar(x, sub[f"x{c}"].values, bottom=bottom, width=0.85,
                   color=PAL[c], label=c)
            bottom += np.nan_to_num(sub[f"x{c}"].values)
        for y in (1/3, 2/3):
            ax.axhline(y, ls="--", lw=0.5, color="grey")
        ax.set_ylim(0, 1); ax.set_ylabel(f"{grp}\nfraction")
    axes[0].legend(title="chr21 copy", ncol=3, fontsize=8, loc="upper right")
    axes[-1].set_xticks(np.arange(len(order)))
    axes[-1].set_xticklabels(order, rotation=90, fontsize=5)
    fig.suptitle("Per-gene chr21 allele fractions", fontsize=11)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig/fig_stacked_full.pdf")
    plt.close(fig)

# =============================================================================
#  Bansal-style 2-way figure: affected copy vs the other two (maximises genes)
#  Needs only the affected copy's own site -> all genes with that site qualify.
# =============================================================================
for grp in comp_groups:
    D = expected.get(grp, "none")
    if D not in COPIES:                  # blind: infer the most-reduced copy
        D = min(COPIES, key=lambda c: dose.loc[grp, c]
                if (grp in dose.index and not np.isnan(dose.loc[grp, c])) else 9)
    pb = pooled[(pooled["group"] == "baseline") & (pooled["resolved_copy"] == D)
                & (pooled["frac"].notna())].set_index("gene")["frac"]
    pdc = pooled[(pooled["group"] == grp) & (pooled["resolved_copy"] == D)
                 & (pooled["frac"].notna())].set_index("gene")["frac"]
    common = pb.index.intersection(pdc.index)
    if len(common) == 0:
        continue
    order = gene_pos.reindex(common).sort_values().index.tolist()
    fig, axes = plt.subplots(2, 1, figsize=(max(8, len(order) * 0.18), 5),
                             sharex=True, sharey=True)
    for ax, (lab, ser) in zip(axes, [("baseline", pb), (grp, pdc)]):
        x = np.arange(len(order)); vals = ser.reindex(order).values
        ax.bar(x, vals, width=0.85, color=PAL[D], label=f"{D} (affected)")
        ax.bar(x, 1 - vals, bottom=vals, width=0.85, color="#BBBBBB",
               label="other two")
        ax.set_ylim(0, 1); ax.set_ylabel(f"{lab}\nfraction")
    axes[0].legend(ncol=2, fontsize=8, loc="upper right")
    axes[-1].set_xticks(np.arange(len(order)))
    axes[-1].set_xticklabels(order, rotation=90, fontsize=5)
    fig.suptitle(f"{grp}: {D}-allele fraction vs the other two "
                 f"(all genes with a {D} site, n={len(order)})", fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig/fig_affected_vs_rest_{grp}.pdf"); plt.close(fig)
    print(f"[05] 2-way figure {grp}: affected={D}, n={len(order)} genes "
          f"-> fig_affected_vs_rest_{grp}.pdf")

print(f"[05] conditions: {comp_groups}")
print(f"[05] genes plotted as full stacks (all groups): {len(keep)}")
print(f"[05] -> {OUT}/fig/dosage_summary.tsv ; fig_class_shift_<cond>.pdf ; "
      f"fig_stacked_full.pdf ; fig_affected_vs_rest_<cond>.pdf")
