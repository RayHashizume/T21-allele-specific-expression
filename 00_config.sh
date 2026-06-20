#!/usr/bin/env bash
# =============================================================================
# 00_config.sh  --  central configuration, sourced by the other scripts.
#   Edit the paths/parameters here only.  `source 00_config.sh` at the top of
#   each step.
# =============================================================================
set -a

# ---- reference --------------------------------------------------------------
REF=/Users/ray/Ref_files/chm13v2.0.fa            # must have .fai and .dict
# chr21 contig name AS IT APPEARS IN YOUR BAMs / REF.  Confirm with:
#   samtools idxstats /Users/ray/Desktop/DS/HiNova.BAM/HiNova.CHM13/dP.HiNova.CHM13.bam | cut -f1 | grep -i 21
# Common possibilities: "chr21" (UCSC/marbl analysis set) | "21" | "NC_060945.1" (RefSeq)
CHR21=chr21

# ---- WGS BAMs (CHM13v2-aligned, indexed) ------------------------------------
BAM_T21=/Users/ray/Desktop/DS/HiNova.BAM/HiNova.CHM13/T21.HiNova.CHM13.bam   # QC only
BAM_dP=/Users/ray/Desktop/DS/HiNova.BAM/HiNova.CHM13/dP.HiNova.CHM13.bam
BAM_dM1=/Users/ray/Desktop/DS/HiNova.BAM/HiNova.CHM13/dM1.HiNova.CHM13.bam
BAM_dM2=/Users/ray/Desktop/DS/HiNova.BAM/HiNova.CHM13/dM2.HiNova.CHM13.bam
BAM_ONT=/Users/ray/Desktop/DS/ONT/Long-T2T/KH2B2_ONT_methyl_CHM13v2.bam      # optional phasing QC

# ---- RNA-seq FASTQ (paired) -------------------------------------------------
RNA_DIR=/Users/ray/Desktop/DS/RNA-seq/2023.3.16_RNA-seq/rawdata
# sample_id<TAB>fastq_R1<TAB>fastq_R2<TAB>group   is built in 04; edit there if names differ.

# ---- annotation -------------------------------------------------------------
# Provide ONE of these. Script 03/04 auto-detect format and remap contig names.
#  (a) RefSeq GTF (for STAR) + GFF (for resolvability), contigs = NC_0609xx.1
#      https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/
#  (b) CAT/Liftoff GENCODEv35 GFF3, contigs = chr1..chrY (matches a chr-named ref)
#      https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_GENCODEv35_CAT_Liftoff.gff3
ANNOT_GTF=/Users/ray/Ref_files/chm13v2.0.gtf       # for STAR --sjdbGTFfile
ANNOT_GFF=/Users/ray/Ref_files/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz
# chromAlias to convert RefSeq accession <-> chrN if your annotation and ref disagree.
# (NCBI: getChromInfoFromNCBI; or UCSC GCA_009914755.4.chromAlias.txt). Leave empty if not needed.
CHROM_ALIAS=

# ---- output -----------------------------------------------------------------
OUT=/Users/ray/Desktop/DS/ASE_out
mkdir -p "$OUT" "$OUT/vcf" "$OUT/snp" "$OUT/rna" "$OUT/ase" "$OUT/fig"

# ---- thresholds -------------------------------------------------------------
MIN_DP_WGS=12          # min per-line depth at a site to trust the disomy genotype
MIN_GQ=20              # min genotype quality
MIN_QUAL=30            # min site QUAL
THREADS=8

# ---- masks (strongly recommended for chr21 acrocentric p-arm / segdups) -----
# chm13v2.0_censat_v2.0.bed (satellite/rDNA) from marbl/CHM13. Sites inside are dropped.
CENSAT_BED=/Users/ray/Ref_files/chm13v2.0_censat_v2.0.bed   # leave nonexistent to skip

# ---- tools (override if not on PATH) ----------------------------------------
SAMTOOLS=samtools
BCFTOOLS=bcftools
BEDTOOLS=bedtools
STAR=STAR
GATK=gatk
PYTHON=python3
RSCRIPT=Rscript

echo "[config] REF=$REF  CHR21=$CHR21  OUT=$OUT" >&2
