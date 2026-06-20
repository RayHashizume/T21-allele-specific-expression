#!/usr/bin/env bash
# =============================================================================
# 01_call_disomy_variants.sh
#   Joint (multi-sample) SNV calling on chr21 for the three DISOMY lines.
#
#   WHY only the disomy lines: dP/dM1/dM2 are bona-fide DIPLOID for chr21
#   (single-copy deletion, confirmed by G-band/FISH/STR), so an ordinary diploid
#   caller is exactly correct.  The parental T21 line carries 3 chr21 copies, so
#   its chr21 allele fractions sit at ~1/3 or ~2/3 and a diploid caller would
#   mis-genotype it.  We therefore use T21 ONLY as an orthogonal QC (every
#   informative site should be heterozygous in T21 with the minority allele
#   near 1/3).  All allele->copy assignment is done from the disomy trio in 02.
#
#   Output: $OUT/vcf/trio.chr21.vcf.gz  with FORMAT/GT,AD,DP for dP,dM1,dM2.
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$HERE/00_config.sh"

# Sample order is fixed and recorded; 02 reads sample names from the VCF header.
$BCFTOOLS mpileup -f "$REF" -r "$CHR21" \
    -a FORMAT/AD,FORMAT/DP,FORMAT/SP \
    -q 20 -Q 20 --threads "$THREADS" \
    "$BAM_dP" "$BAM_dM1" "$BAM_dM2" \
  | $BCFTOOLS call -m -v --threads "$THREADS" -Oz -o "$OUT/vcf/trio.raw.vcf.gz"

# Force readable sample names (mpileup uses file paths / RG by default).
printf '%s\n%s\n%s\n' dP dM1 dM2 > "$OUT/vcf/samples.txt"
$BCFTOOLS reheader -s "$OUT/vcf/samples.txt" "$OUT/vcf/trio.raw.vcf.gz" \
  | $BCFTOOLS norm -f "$REF" -m -any -Ou \
  | $BCFTOOLS view -v snps -m2 -M2 \
        -i "QUAL>=$MIN_QUAL" -Oz -o "$OUT/vcf/trio.chr21.vcf.gz"
$BCFTOOLS index -t "$OUT/vcf/trio.chr21.vcf.gz"

echo "[01] biallelic SNVs on $CHR21:" >&2
$BCFTOOLS view -H "$OUT/vcf/trio.chr21.vcf.gz" | wc -l >&2

# --- optional: T21 QC pileup at the same sites (allele fractions) ------------
# Triploid AF check is done in 02 if this file is present.
$BCFTOOLS mpileup -f "$REF" -r "$CHR21" -a FORMAT/AD,FORMAT/DP \
    -q 20 -Q 20 --threads "$THREADS" \
    -T "$OUT/vcf/trio.chr21.vcf.gz" "$BAM_T21" \
  | $BCFTOOLS call -m -Oz -o "$OUT/vcf/T21.qc.vcf.gz" || \
  echo "[01] T21 QC pileup skipped (non-fatal)" >&2
$BCFTOOLS index -t "$OUT/vcf/T21.qc.vcf.gz" 2>/dev/null || true

echo "[01] done -> $OUT/vcf/trio.chr21.vcf.gz" >&2
