#!/usr/bin/env bash
# =============================================================================
# 04_align_count.sh   (HISAT2; manifest-driven, supports any deletion condition)
#   Aligns every RNA-seq sample listed in $OUT/rna/manifest.tsv to CHM13v2 with
#   HISAT2, then counts ref/alt reads at the informative sites with GATK
#   ASEReadCounter. A default manifest (baseline + delM2) is created on first
#   run; ADD ROWS for future delP / delM1 experiments and re-run -- existing
#   BAMs are reused, only new samples are aligned.
#
#   manifest.tsv columns (tab-separated, '#'=comment, no header row):
#       sid   R1   R2   condition   deleted_copy(none|P|M1|M2)
#
#   Reference-bias note + index-without-annotation rationale: see prior header.
#   Output: $OUT/ase/<sid>.ase.table
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$HERE/00_config.sh"

HIDX="$OUT/rna/hisat2_idx/chm13v2"
HETVCF="$OUT/snp/hetsites.vcf.gz"
DICT="${REF%.fa}.dict"
MANIFEST="$OUT/rna/manifest.tsv"

# ---- default manifest on first run (edit/add rows for new experiments) ------
if [ ! -e "$MANIFEST" ]; then
  cat > "$MANIFEST" <<EOF
# sid	R1	R2	condition	deleted_copy(none|P|M1|M2)
# Add rows for future deletion experiments (e.g. delP / delM1) and re-run.
KH2B2_1	$RNA_DIR/KH2B2_1_1.fq.gz	$RNA_DIR/KH2B2_1_2.fq.gz	baseline	none
KH2B2_2	$RNA_DIR/KH2B2_2_1.fq.gz	$RNA_DIR/KH2B2_2_2.fq.gz	baseline	none
KH2B2_3	$RNA_DIR/KH2B2_3_1.fq.gz	$RNA_DIR/KH2B2_3_2.fq.gz	baseline	none
deltaM2_8	$RNA_DIR/deltaM2_8_1.fq.gz	$RNA_DIR/deltaM2_8_2.fq.gz	delM2	M2
deltaM2_14	$RNA_DIR/deltaM2_14_1.fq.gz	$RNA_DIR/deltaM2_14_2.fq.gz	delM2	M2
deltaM2_21	$RNA_DIR/deltaM2_21_1.fq.gz	$RNA_DIR/deltaM2_21_2.fq.gz	delM2	M2
EOF
  echo "[04] wrote default manifest -> $MANIFEST  (edit to add experiments)" >&2
fi

# ---- GATK sequence dictionary -----------------------------------------------
[ -e "$DICT" ] || [ -e "$REF.dict" ] || $GATK CreateSequenceDictionary -R "$REF"

# ---- build HISAT2 index once ------------------------------------------------
if [ ! -e "${HIDX}.1.ht2" ]; then
  mkdir -p "$(dirname "$HIDX")"
  echo "[04] building HISAT2 index (one-time)..." >&2
  hisat2-build -p "$THREADS" "$REF" "$HIDX"
fi

[ -e "${HETVCF}.tbi" ] || $BCFTOOLS index -t "$HETVCF"

# ---- align + count every manifest sample ------------------------------------
grep -v '^#' "$MANIFEST" | while IFS=$'\t' read -r SID R1 R2 COND DEL; do
  [ -z "${SID:-}" ] && continue
  echo "[04] processing $SID (condition=$COND deleted=$DEL)" >&2
  PFX="$OUT/rna/$SID."
  if [ ! -e "${PFX}bam" ]; then
    echo "[04]   aligning $SID" >&2
    hisat2 -p "$THREADS" -x "$HIDX" -1 "$R1" -2 "$R2" \
           --no-unal --rg-id "$SID" --rg "SM:$SID" \
      | $SAMTOOLS sort -@ "$THREADS" -o "${PFX}bam" -
    $SAMTOOLS index "${PFX}bam"
  else
    echo "[04]   ${PFX}bam exists -- reusing" >&2
    [ -e "${PFX}bam.bai" ] || $SAMTOOLS index "${PFX}bam"
  fi
  $GATK ASEReadCounter -R "$REF" -I "${PFX}bam" -V "$HETVCF" \
        -L "$CHR21" --min-mapping-quality 60 --min-base-quality 20 \
        -O "$OUT/ase/$SID.ase.table"
done

echo "[04] done -> $OUT/ase/*.ase.table" >&2
