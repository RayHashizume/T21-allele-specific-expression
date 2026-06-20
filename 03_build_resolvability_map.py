#!/usr/bin/env python3
# =============================================================================
# 03_build_resolvability_map.py   (v2: fast overlap + gene-level biotype)
#   Assign each informative SNP (from 02) to overlapping transcribed features
#   (exons -- which include UTRs -- and ncRNA exons) and summarise, per gene,
#   which of P / M1 / M2 are individually resolvable.
#
#   Resolvability rule (biallelic reality, fractions sum to 1):
#     - sites of >=2 distinct classes in a gene -> ALL THREE estimable ("full")
#     - only ONE class present -> that copy vs the merged rest ("partial:<copy>")
#     - no informative site -> "none"
#
#   Handles RefSeq GFF (seqid NC_060945.1) and chr-named GFF/GTF alike; the
#   chr21 contig is recognised via aliases and CHROM_ALIAS.
#
#   Outputs:
#     $OUT/snp/sites_with_gene.tsv      informative sites + gene/biotype
#     $OUT/snp/gene_resolvability.tsv   per-gene counts + status
# =============================================================================
import os, sys, gzip, bisect
from collections import defaultdict

OUT     = os.environ.get("OUT", os.path.expanduser("~/Desktop/DS/ASE_out"))
CHR21   = os.environ.get("CHR21", "chr21")
GFF     = os.environ.get("ANNOT_GFF", "")
ALIAS   = os.environ.get("CHROM_ALIAS", "")
SITES   = f"{OUT}/snp/informative_sites.tsv"

CHR21_ALIASES = {CHR21, "chr21", "21", "NC_060945.1"}
if ALIAS and os.path.exists(ALIAS):
    with open(ALIAS) as fh:
        for ln in fh:
            names = set(ln.rstrip("\n").split("\t"))
            if names & CHR21_ALIASES:
                CHR21_ALIASES |= names

def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)

def parse_attr(field, fmt):
    d = {}
    if fmt == "gff3":
        for kv in field.rstrip(";").split(";"):
            if "=" in kv:
                k, v = kv.split("=", 1); d[k.strip()] = v.strip()
    else:
        for kv in field.strip().rstrip(";").split(";"):
            kv = kv.strip()
            if " " in kv:
                k, v = kv.split(" ", 1); d[k.strip()] = v.strip().strip('"')
    return d

def gene_of(d):
    return (d.get("gene_name") or d.get("gene") or d.get("gene_id")
            or d.get("Name") or d.get("Parent") or "NA")

def biotype_of(d):
    return (d.get("gene_biotype") or d.get("gene_type") or d.get("biotype")
            or d.get("gbkey") or "NA")

def load_annotation():
    if not GFF or not os.path.exists(GFF):
        sys.exit(f"[03] annotation not found: {GFF}  (set ANNOT_GFF in 00_config.sh)")
    fmt = "gff3" if (".gff" in GFF) else "gtf"
    exons = []                       # (start, end, gene)
    gene_biotype = {}                # gene -> biotype (from 'gene' feature lines)
    with opener(GFF) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] not in CHR21_ALIASES:
                continue
            ft = f[2]
            if ft == "gene":
                d = parse_attr(f[8], fmt)
                gene_biotype[gene_of(d)] = biotype_of(d)
            elif ft == "exon":
                d = parse_attr(f[8], fmt)
                exons.append((int(f[3]), int(f[4]), gene_of(d)))
    if not exons:
        sys.exit("[03] no chr21 exons parsed -- check contig naming / CHROM_ALIAS.")
    exons.sort()
    starts = [e[0] for e in exons]
    ends   = [e[1] for e in exons]
    genes  = [e[2] for e in exons]
    max_len = max(b - a for a, b in zip(starts, ends))
    sys.stderr.write(f"[03] chr21 exons={len(exons)}  genes_with_biotype="
                     f"{len(gene_biotype)}  max_exon_len={max_len}\n")
    return starts, ends, genes, gene_biotype, max_len

def main():
    starts, ends, genes, gbt, max_len = load_annotation()
    n_lines = sum(1 for _ in open(SITES)) - 1

    per_gene = defaultdict(lambda: {"P": 0, "M1": 0, "M2": 0})
    out_sites = open(f"{OUT}/snp/sites_with_gene.tsv", "w")
    done = assigned = 0
    with open(SITES) as fh:
        header = fh.readline().rstrip("\n")
        out_sites.write(header + "\tgene\tbiotype\n")
        for ln in fh:
            done += 1
            if done % 5000 == 0:
                sys.stderr.write(f"\r[03] assigning sites {done}/{n_lines}")
                sys.stderr.flush()
            f = ln.rstrip("\n").split("\t")
            pos = int(f[1]); copy = f[4]
            lo = bisect.bisect_left(starts, pos - max_len)
            hi = bisect.bisect_right(starts, pos)
            hits = [genes[j] for j in range(lo, hi)
                    if starts[j] <= pos <= ends[j]]
            if not hits:
                out_sites.write(ln.rstrip("\n") + "\tNA\tNA\n")
                continue
            ghit = sorted(set(hits))
            bhit = sorted({gbt.get(g, "NA") for g in ghit})
            out_sites.write(ln.rstrip("\n") + "\t" + ",".join(ghit) +
                            "\t" + ",".join(bhit) + "\n")
            assigned += 1
            for g in ghit:
                per_gene[g][copy] += 1
    out_sites.close()
    sys.stderr.write(f"\r[03] sites assigned to >=1 gene: {assigned}/{n_lines}\n")

    full = part = 0
    with open(f"{OUT}/snp/gene_resolvability.tsv", "w") as out:
        out.write("gene\tbiotype\tn_P\tn_M1\tn_M2\tn_classes\tstatus\n")
        for g, d in sorted(per_gene.items()):
            present = [c for c in ("P", "M1", "M2") if d[c] > 0]
            n = len(present)
            if n >= 2:
                status = "full"; full += 1
            elif n == 1:
                status = f"partial:{present[0]}"; part += 1
            else:
                status = "none"
            out.write(f"{g}\t{gbt.get(g,'NA')}\t{d['P']}\t{d['M1']}\t{d['M2']}"
                      f"\t{n}\t{status}\n")
    sys.stderr.write(
        f"[03] genes fully resolvable (P,M1,M2 separable): {full};  "
        f"partial (one copy vs rest): {part}\n"
        f"[03] -> {OUT}/snp/sites_with_gene.tsv\n"
        f"[03] -> {OUT}/snp/gene_resolvability.tsv\n")

if __name__ == "__main__":
    main()
