#!/usr/bin/env python3
# =============================================================================
# 07_absolute_expression.py
#   ABSOLUTE per-allele expression (not just fractions). For each gene:
#       absolute_C = (normalized gene expression) x (directly-measured fraction of copy C)
#   so that e.g. M2 deletion shows P and M1 ~unchanged while M2 -> 0.
#
#   Gene expression is counted with featureCounts over chr21 exons and
#   normalized by NON-chr21 mapped reads (robust to the chr21 dosage change
#   itself). Allele fractions use ONLY directly-measured class sites -- no
#   residual-fill -- so a copy without its own site is simply not plotted for
#   that gene (it never invents a segment).
#
#   Deps: featureCounts (subread) + samtools.  Install: mamba install -y -c bioconda subread
#   Outputs ($OUT/expr): gene_counts.txt, libsizes.tsv, absolute_expression.tsv,
#                        fig_absolute_<cond>.pdf, fig_absolute_stacked.pdf
# =============================================================================
import os, sys, gzip, subprocess, shlex
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT     = os.environ.get("OUT", os.path.expanduser("~/Desktop/DS/ASE_out"))
REF     = os.environ["REF"]
CHR21   = os.environ.get("CHR21", "chr21")
GFF     = os.environ["ANNOT_GFF"]
NCID    = os.environ.get("CHR21_NCID", "NC_060945.1")  # GFF contig name for chr21
THREADS = os.environ.get("THREADS", "8")
MINTOT  = int(os.environ.get("ASE_MIN_TOTAL", "20"))
COPIES  = ["P", "M1", "M2"]
PAL     = {"P": "#3B6FB6", "M1": "#E08214", "M2": "#7B3294"}
EXPR    = f"{OUT}/expr"; os.makedirs(EXPR, exist_ok=True)
SAF     = f"{EXPR}/chr21.saf"
COUNTS  = f"{EXPR}/gene_counts.txt"
LIBS    = f"{EXPR}/libsizes.tsv"

man = pd.read_csv(f"{OUT}/rna/manifest.tsv", sep="\t", comment="#", header=None,
                  names=["sid", "r1", "r2", "condition", "deleted_copy"])
man["deleted_copy"] = man["deleted_copy"].astype(str)
man["group"] = np.where(man["deleted_copy"] == "none", "baseline", man["condition"])
bams = {r.sid: f"{OUT}/rna/{r.sid}.bam" for r in man.itertuples()}

# ---- 1. SAF (chr21 exons; rename GFF contig -> CHR21) ------------------------
def attr(s, key):
    for f in s.split(";"):
        if f.startswith(key + "="):
            return f[len(key) + 1:]
    return None

if not os.path.exists(SAF):
    op = gzip.open if GFF.endswith(".gz") else open
    n = 0
    with op(GFF, "rt") as fh, open(SAF, "w") as out:
        out.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] != NCID or f[2] != "exon":
                continue
            g = attr(f[8], "gene") or attr(f[8], "gene_id") or attr(f[8], "Parent")
            if not g:
                continue
            out.write(f"{g}\t{CHR21}\t{f[3]}\t{f[4]}\t{f[6]}\n"); n += 1
    sys.stderr.write(f"[07] wrote SAF with {n} chr21 exon rows\n")

# ---- 2. featureCounts (paired, fragment counts) -----------------------------
if not os.path.exists(COUNTS):
    bam_list = " ".join(shlex.quote(bams[s]) for s in man["sid"])
    cmd = (f"featureCounts -F SAF -a {shlex.quote(SAF)} -o {shlex.quote(COUNTS)} "
           f"-T {THREADS} -p --countReadPairs {bam_list}")
    sys.stderr.write(f"[07] {cmd}\n")
    subprocess.run(cmd, shell=True, check=True)

# ---- 3. library sizes = NON-chr21 mapped reads (samtools idxstats) ----------
if not os.path.exists(LIBS):
    rows = []
    for s in man["sid"]:
        idx = subprocess.run(["samtools", "idxstats", bams[s]],
                             capture_output=True, text=True, check=True).stdout
        nonchr21 = sum(int(l.split("\t")[2]) for l in idx.strip().splitlines()
                       if l.split("\t")[0] != CHR21)
        rows.append({"sid": s, "nonchr21_mapped": nonchr21})
    pd.DataFrame(rows).to_csv(LIBS, sep="\t", index=False)
libs = pd.read_csv(LIBS, sep="\t").set_index("sid")["nonchr21_mapped"]

# ---- 4. normalized gene expression per sample -------------------------------
fc = pd.read_csv(COUNTS, sep="\t", comment="#")
meta = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
count_cols = [c for c in fc.columns if c not in meta]
ren = {c: os.path.basename(c).replace(".bam", "") for c in count_cols}
fc = fc.rename(columns=ren)
expr = fc.set_index("Geneid")[list(ren.values())]
# normalize: counts per million non-chr21 mapped reads
E = expr.divide([libs[s] for s in expr.columns], axis=1) * 1e6   # gene x sample
# mean per group
sid2grp = dict(zip(man["sid"], man["group"]))
Eg = E.T.groupby([sid2grp[s] for s in E.columns]).mean().T            # gene x group

# ---- 5. directly-measured allele fractions per (group, gene, class) ---------
sites = pd.read_csv(f"{OUT}/snp/sites_with_gene.tsv", sep="\t")
sites = sites[sites["gene"].notna() & (sites["gene"] != "NA")].copy()
sites["gene"] = sites["gene"].str.split(",")
sites = sites.explode("gene", ignore_index=True)
sites = sites[["pos", "resolved_copy", "odd_is_alt", "gene"]]

def load_ase(sid):
    t = pd.read_csv(f"{OUT}/ase/{sid}.ase.table", sep="\t").rename(
        columns={"position": "pos"})
    m = sites.merge(t[["pos", "refCount", "altCount", "totalCount"]], on="pos")
    m["odd"] = np.where(m["odd_is_alt"] == 1, m["altCount"], m["refCount"])
    m["tot"] = m["totalCount"]; m["group"] = sid2grp[sid]
    return m[m["tot"] > 0]

ase = pd.concat([load_ase(s) for s in man["sid"]], ignore_index=True)
frac = (ase.groupby(["group", "gene", "resolved_copy"], as_index=False)
           .agg(odd=("odd", "sum"), tot=("tot", "sum")))
frac = frac[frac["tot"] >= MINTOT]
frac["frac"] = frac["odd"] / frac["tot"]

# ---- 6. absolute per-allele expression = E_gene(group) x measured fraction --
recs = []
for r in frac.itertuples():
    if r.gene in Eg.index and r.group in Eg.columns:
        recs.append({"gene": r.gene, "group": r.group, "copy": r.resolved_copy,
                     "frac": r.frac, "gene_expr": Eg.loc[r.gene, r.group],
                     "abs_expr": Eg.loc[r.gene, r.group] * r.frac})
absdf = pd.DataFrame(recs)
absdf.to_csv(f"{OUT}/expr/absolute_expression.tsv", sep="\t", index=False)

# ---- 7. per-class figure: baseline vs each condition (absolute) -------------
comp = [g for g in man["group"].unique() if g != "baseline"]
for grp in comp:
    fig, axes = plt.subplots(1, 3, figsize=(9, 4), sharey=True)
    for ax, cls in zip(axes, COPIES):
        b = absdf[(absdf["group"] == "baseline") & (absdf["copy"] == cls)] \
            .set_index("gene")["abs_expr"]
        d = absdf[(absdf["group"] == grp) & (absdf["copy"] == cls)] \
            .set_index("gene")["abs_expr"]
        common = b.index.intersection(d.index)
        for g in common:                       # paired lines per gene
            ax.plot([1, 2], [b[g], d[g]], color=PAL[cls], alpha=0.3, lw=0.6,
                    zorder=1)
        ax.scatter(np.ones(len(common)), b[common], s=12, color=PAL[cls], zorder=3)
        ax.scatter(np.full(len(common), 2), d[common], s=12, color=PAL[cls],
                   zorder=3)
        ax.set_xticks([1, 2]); ax.set_xticklabels(["baseline", grp])
        ax.set_title(f"{cls}-allele (n={len(common)} genes)")
    axes[0].set_ylabel("absolute allele expression\n(per 1e6 non-chr21 reads)")
    fig.suptitle(f"{grp}: absolute per-allele expression "
                 f"(deleted copy -> 0; others ~unchanged)", fontsize=11)
    fig.tight_layout(); fig.savefig(f"{OUT}/expr/fig_absolute_{grp}.pdf")
    plt.close(fig)

# ---- 8. per-gene absolute stacked bars (deleted copy direct; <=1 retained residual)
expected = (man.groupby("group")["deleted_copy"]
            .agg(lambda s: s.mode().iat[0]).to_dict())

def complete_fracs(d, deleted):
    """deleted copy must be directly measured; <=1 retained copy may be residual."""
    missing = [c for c in COPIES if c not in d]
    if deleted in COPIES and deleted in missing:
        return None
    if len(missing) == 0:
        v = dict(d)
    elif len(missing) == 1:
        v = dict(d); v[missing[0]] = max(0.0, 1 - sum(d.values()))
    else:
        return None
    s = sum(v.values())
    return {c: v[c] / s for c in COPIES} if s > 0 else None

fr_by = {(g, gene): dict(zip(s["resolved_copy"], s["frac"]))
         for (g, gene), s in frac.groupby(["group", "gene"])}
groups_all = ["baseline"] + comp
absrows = []
for (grp, gene), d in fr_by.items():
    if gene not in Eg.index or grp not in Eg.columns:
        continue
    deleted = expected.get(grp, "none")
    strict = deleted not in COPIES and deleted != "none"
    if strict:
        if len(d) < 3:
            continue
        ssum = sum(d.values()); vv = {c: d[c] / ssum for c in COPIES}
    else:
        vv = complete_fracs(d, deleted if deleted in COPIES else None)
        if vv is None:
            continue
    E = Eg.loc[gene, grp]
    absrows.append({"gene": gene, "group": grp,
                    "P": E * vv["P"], "M1": E * vv["M1"], "M2": E * vv["M2"]})
pv = pd.DataFrame(absrows)
keep = ([g for g, s in pv.groupby("gene")
         if set(groups_all).issubset(set(s["group"]))] if not pv.empty else [])
pf = pv[pv["gene"].isin(keep)]
if not pf.empty:
    order = (pf.drop_duplicates("gene")
             .merge(sites.groupby("gene")["pos"].min().rename("pos"), on="gene")
             .sort_values("pos")["gene"].tolist())
    fig, axes = plt.subplots(len(groups_all), 1,
                             figsize=(max(8, len(order) * 0.22), 2.2 * len(groups_all)),
                             sharex=True, sharey=True)
    if len(groups_all) == 1:
        axes = [axes]
    for ax, grp in zip(axes, groups_all):
        sub = pf[pf["group"] == grp].set_index("gene").reindex(order)
        x = np.arange(len(order)); bottom = np.zeros(len(order))
        for c in COPIES:
            ax.bar(x, sub[c].values, bottom=bottom, width=0.85, color=PAL[c],
                   label=c)
            bottom += np.nan_to_num(sub[c].values)
        ax.set_ylabel(f"{grp}\nabs. expr")
    axes[0].legend(title="chr21 copy", ncol=3, fontsize=8, loc="upper right")
    axes[-1].set_xticks(np.arange(len(order)))
    axes[-1].set_xticklabels(order, rotation=90, fontsize=5)
    fig.suptitle("Per-gene ABSOLUTE allele expression "
                 "(bar height = total; deleted copy removed)", fontsize=11)
    fig.tight_layout(); fig.savefig(f"{OUT}/expr/fig_absolute_stacked.pdf")
    plt.close(fig)

# ---- console summary: PAIRED comparison (genes measured in both groups) -----
# total chr21 expression per group (unique gene-level, summed) -> dosage ratio
gtot = (absdf.drop_duplicates(["gene", "group"])
              .groupby("group")["gene_expr"].sum())
print("[07] PAIRED per-allele absolute expression (common genes only):")
for grp in comp:
    base_tot = gtot.get("baseline", np.nan); grp_tot = gtot.get(grp, np.nan)
    print(f"  -- {grp} vs baseline --   total chr21 dosage ratio "
          f"({grp}/baseline) = {grp_tot/base_tot:.2f}  "
          f"(1.00=full compensation, 0.67=no compensation)")
    for c in COPIES:
        b = absdf[(absdf["group"] == "baseline") & (absdf["copy"] == c)] \
            .set_index("gene")["abs_expr"]
        d = absdf[(absdf["group"] == grp) & (absdf["copy"] == c)] \
            .set_index("gene")["abs_expr"]
        common = b.index.intersection(d.index)
        if len(common) == 0:
            print(f"     {c}: (no common genes)"); continue
        ratio = (d[common] / b[common].replace(0, np.nan)).median()
        print(f"     {c}: baseline med={b[common].median():.2f}  "
              f"{grp} med={d[common].median():.2f}  "
              f"median paired ratio={ratio:.2f}  (n={len(common)} genes)")
print(f"[07] -> {OUT}/expr/absolute_expression.tsv ; fig_absolute_<cond>.pdf ; "
      f"fig_absolute_stacked.pdf")
