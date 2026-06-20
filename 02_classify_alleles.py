#!/usr/bin/env python3
# =============================================================================
# 02_classify_alleles.py
#   CORE STEP.  For every biallelic chr21 SNV, decide which single copy
#   (P, M1, or M2) carries the "private" (minority) allele, using the genotypes
#   of the three disomy lines.
#
#   Biology recap (this line is MI-nondisjunction -> M1 and M2 are distinct
#   homologs, heterozygous chr21-wide, NOT identical sister chromatids):
#       dP  observes the diploid genotype of {M1, M2}
#       dM1 observes the diploid genotype of {P,  M2}
#       dM2 observes the diploid genotype of {P,  M1}
#
#   At a biallelic site the three copies hold only two bases, so the only
#   informative pattern is a 2-1 split: exactly ONE copy is the odd one out.
#   We recover the per-copy alleles by enumerating the 8 possible
#   (aP, aM1, aM2) in {0,1}^3, predicting each line's UNORDERED genotype, and
#   keeping the assignment(s) consistent with what was observed.  A site is
#   accepted only if all consistent assignments agree on the SAME odd copy and
#   the SAME odd allele; otherwise it is dropped as ambiguous.
#
#   This deletion-mapping is deliberately the same logic as the experiment:
#   the copy's private allele is the one that DISAPPEARS in the line where that
#   copy was deleted and is RETAINED in the other two lines.  So at the RNA
#   level the M2-deleted cells should show the M2-private allele drop to ~0 at
#   M2-informative sites.
#
#   Outputs:
#     $OUT/snp/informative_sites.tsv   one row per accepted site:
#         chrom pos ref alt resolved_copy odd_allele odd_is_alt
#         merged_copies aP aM1 aM2 dp_dP dp_dM1 dp_dM2
#       where odd_is_alt=1 means the odd copy carries ALT (so in RNA counts the
#       odd-copy expression = altCount), 0 means it carries REF.
#     $OUT/snp/hetsites.vcf.gz (+.tbi)  the same positions, GT=0/1, one dummy
#       sample, for STAR-WASP and GATK ASEReadCounter.
#
#   Requires: pysam  (pip install pysam --break-system-packages)
# =============================================================================
import os, sys, itertools, gzip
import pysam

OUT      = os.environ.get("OUT", os.path.expanduser("~/Desktop/DS/ASE_out"))
CHR21    = os.environ.get("CHR21", "chr21")
MIN_DP   = int(os.environ.get("MIN_DP_WGS", "12"))
MIN_GQ   = int(os.environ.get("MIN_GQ", "20"))
CENSAT   = os.environ.get("CENSAT_BED", "")
# max fraction of a copy's private allele tolerated in the deletion line's WGS;
# above this the locus mismaps and the site is dropped (see WGS allele-depth gate)
WGS_AD_THR = float(os.environ.get("WGS_AD_THR", "0.05"))

VCF_IN   = f"{OUT}/vcf/trio.chr21.vcf.gz"
TSV_OUT  = f"{OUT}/snp/informative_sites.tsv"
VCF_OUT  = f"{OUT}/snp/hetsites.vcf.gz"

COPIES = ("P", "M1", "M2")

# ---- optional satellite/rDNA mask (drop sites inside) -----------------------
def load_mask(bed, chrom):
    ivs = []
    if bed and os.path.exists(bed):
        op = gzip.open if bed.endswith(".gz") else open
        with op(bed, "rt") as fh:
            for ln in fh:
                if ln.startswith(("#", "track", "browser")):
                    continue
                f = ln.split("\t")
                if len(f) >= 3 and f[0] == chrom:
                    ivs.append((int(f[1]), int(f[2])))
        ivs.sort()
    return ivs

def in_mask(ivs, pos):  # pos is 1-based; bed is 0-based half-open
    p = pos - 1
    lo, hi = 0, len(ivs)
    while lo < hi:
        mid = (lo + hi) // 2
        if ivs[mid][1] <= p:
            lo = mid + 1
        else:
            hi = mid
    return lo < len(ivs) and ivs[lo][0] <= p < ivs[lo][1]

# ---- per-copy resolution from the disomy trio genotypes ---------------------
def content_set(gt):
    """gt is a tuple of allele indices, e.g. (0,1). Returns frozenset of alleles."""
    a = [x for x in gt if x is not None]
    if len(a) != 2:           # require a confident diploid call
        return None
    return frozenset(a)

def resolve(c_dP, c_dM1, c_dM2):
    """Return (resolved_copy, odd_allele, (aP,aM1,aM2)) or None if uninformative
    /ambiguous.  c_* are frozensets of {0,1}. odd_allele is 0(REF) or 1(ALT)."""
    cands = []
    for aP, aM1, aM2 in itertools.product((0, 1), repeat=3):
        if (frozenset((aM1, aM2)) == c_dP and
            frozenset((aP,  aM2)) == c_dM1 and
            frozenset((aP,  aM1)) == c_dM2):
            cands.append((aP, aM1, aM2))
    if not cands:
        return None
    decided = set()
    for combo in cands:
        # which allele is the minority (appears exactly once)?
        ones = sum(combo)
        if ones == 1:
            odd_allele = 1
            odd_idx = combo.index(1)
        elif ones == 2:
            odd_allele = 0
            odd_idx = combo.index(0)
        else:
            return None                     # monomorphic -> not informative
        decided.add((COPIES[odd_idx], odd_allele, combo))
    # require unanimity on (copy, allele); allow combos to differ only in the
    # identity of the two MERGED copies (which a biallelic site cannot separate)
    keys = {(d[0], d[1]) for d in decided}
    if len(keys) != 1:
        return None
    copy, allele = next(iter(keys))
    combo = sorted(decided)[0][2]
    return copy, allele, combo

# ---- main -------------------------------------------------------------------
def main():
    os.makedirs(f"{OUT}/snp", exist_ok=True)
    mask = load_mask(CENSAT, CHR21)
    vcf = pysam.VariantFile(VCF_IN)
    smp = list(vcf.header.samples)
    for s in ("dP", "dM1", "dM2"):
        if s not in smp:
            sys.exit(f"[02] sample '{s}' not in VCF header {smp}; fix 01 reheader.")

    # Pre-scan: positions carrying >1 variant record (multiallelic sites split
    # by `bcftools norm -m -any`, i.e. effectively tri-allelic among samples).
    # Our P/M1/M2 model is biallelic and GATK ASEReadCounter requires exactly
    # one variant per position, so drop these positions entirely.
    from collections import Counter
    posc = Counter()
    v0 = pysam.VariantFile(VCF_IN)
    for r in v0.fetch(CHR21):
        posc[r.pos] += 1
    v0.close()
    dup_pos = {p for p, c in posc.items() if c > 1}

    # writer for the het-sites VCF (single dummy sample, GT 0/1)
    hdr = pysam.VariantHeader()
    for c in vcf.header.contigs.values():
        hdr.contigs.add(c.name, length=c.length)
    hdr.add_meta('FORMAT', items=[('ID', 'GT'), ('Number', '1'),
                                  ('Type', 'String'), ('Description', 'Genotype')])
    hdr.add_sample('ASE')
    vout = pysam.VariantFile(VCF_OUT, "wz", header=hdr)

    n_in = n_masked = n_ambig = n_ok = n_dup = n_adfail = 0
    counts = {"P": 0, "M1": 0, "M2": 0}
    with open(TSV_OUT, "w") as out:
        out.write("chrom\tpos\tref\talt\tresolved_copy\todd_allele\todd_is_alt"
                  "\tmerged_copies\taP\taM1\taM2\tdp_dP\tdp_dM1\tdp_dM2\n")
        for rec in vcf.fetch(CHR21):
            n_in += 1
            if rec.pos in dup_pos:
                n_dup += 1
                continue                                  # multi-record position
            if len(rec.ref) != 1 or rec.alts is None or len(rec.alts) != 1 \
               or len(rec.alts[0]) != 1:
                continue                                  # SNV only
            if mask and in_mask(mask, rec.pos):
                n_masked += 1
                continue
            g = {s: rec.samples[s] for s in ("dP", "dM1", "dM2")}
            # quality / depth gates
            ok = True
            dps = {}
            for s in ("dP", "dM1", "dM2"):
                dp = g[s].get("DP")
                gq = g[s].get("GQ")
                dps[s] = dp if dp is not None else 0
                if dp is None or dp < MIN_DP:
                    ok = False
                if gq is not None and gq < MIN_GQ:
                    ok = False
            if not ok:
                continue
            cs = {s: content_set(g[s]["GT"]) for s in ("dP", "dM1", "dM2")}
            if any(v is None for v in cs.values()):
                continue
            res = resolve(cs["dP"], cs["dM1"], cs["dM2"])
            if res is None:
                n_ambig += 1
                continue
            copy, odd_allele, combo = res

            # WGS allele-depth gate (ground truth from genomic DNA):
            # the line that DELETED this copy must carry ~0 reads of the copy's
            # private allele. If it does carry the allele, the locus mismaps
            # (paralog/repeat/CN), so the site is unreliable -> drop. This is
            # independent of the RNA result (non-circular) and also cleans the
            # P and M1 assignments, not only M2.
            del_line = "d" + copy            # P->dP, M1->dM1, M2->dM2
            ad = g[del_line].get("AD")
            if ad is None or len(ad) < 2:
                n_adfail += 1
                continue
            priv = ad[odd_allele]            # private-allele reads in deletion line
            tot_ad = ad[0] + ad[1]
            priv_frac = (priv / tot_ad) if tot_ad > 0 else 1.0
            if priv_frac > WGS_AD_THR:
                n_adfail += 1
                continue

            odd_is_alt = 1 if odd_allele == 1 else 0
            merged = [c for c in COPIES if c != copy]
            counts[copy] += 1
            n_ok += 1
            out.write("\t".join(map(str, [
                rec.chrom, rec.pos, rec.ref, rec.alts[0], copy, odd_allele,
                odd_is_alt, ",".join(merged),
                combo[0], combo[1], combo[2],
                dps["dP"], dps["dM1"], dps["dM2"]])) + "\n")
            # write het site
            nr = vout.new_record(contig=rec.chrom, start=rec.start,
                                 alleles=(rec.ref, rec.alts[0]))
            nr.samples['ASE']['GT'] = (0, 1)
            vout.write(nr)

    vout.close()
    pysam.tabix_index(VCF_OUT, preset="vcf", force=True)
    sys.stderr.write(
        f"[02] sites in={n_in} multi_record_dropped={n_dup} masked={n_masked} "
        f"ambiguous_dropped={n_ambig} wgs_ad_gate_dropped={n_adfail} "
        f"accepted={n_ok}\n[02] informative per copy: "
        f"P={counts['P']}  M1={counts['M1']}  M2={counts['M2']}\n"
        f"[02] -> {TSV_OUT}\n[02] -> {VCF_OUT}\n")

if __name__ == "__main__":
    main()
