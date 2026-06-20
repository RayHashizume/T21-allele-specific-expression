# T21_ASE_pipeline
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20772770.svg)](https://doi.org/10.5281/zenodo.20772770)
Allele-specific quantification of the three chromosome 21 homologs in trisomy 21
(T21) induced pluripotent stem cells from RNA-seq, for validating
allele-specific chromosome 21 elimination (or transcriptional inactivation) at
the level of expressed allelic dosage.

Given RNA-seq from a parental trisomy line and from an edited line, the pipeline
reports, for each chr21 gene and chromosome-wide, the fraction and the absolute
expression contributed by each of the three homologs (P, M1, M2), together with
replicate-level statistics.

---

## Overview

The trisomy 21 line carries three chr21 homologs: **P** (paternal) and **M1**
and **M2** (the two maternal homologs). Because the trisomy arose through a
meiosis I non-disjunction event, M1 and M2 are non-identical homologs and are
therefore distinguishable along the whole chromosome.

Three isogenic **induced-disomy** lines, each retaining a defined pair of
homologs after elimination of a single chr21 copy — **dP** (M1+M2), **dM1**
(P+M2) and **dM2** (P+M1) — are used to phase chr21 variants. The single-
chromosome composition of each disomy line was confirmed by G-banding, FISH and
STR analysis.

**Phasing principle.** The *private allele* of a copy C is the allele that is
absent from the disomy line lacking C and present in the other two; at a
biallelic site this resolves copy C against the other two.

**Mechanism note.** RNA-seq quantifies expressed allelic dosage and does **not**
distinguish physical chromosome elimination from transcriptional silencing of an
intact homolog (both abolish the corresponding allele's transcripts). The
pipeline treats deletion and inactivation experiments identically; the mechanism
must be established by orthogonal DNA/cytogenetic assays.

---

## Requirements

A conda (mambaforge) environment. All Python steps must be run with the conda
`python` interpreter (which provides `pysam`), not a system Python.

| Tool | Version used |
|------|--------------|
| bcftools | 1.20 (HTSlib 1.23.1) |
| SAMtools | 1.18 |
| HISAT2 | 2.2.2 |
| GATK | 4.4.0.0 |
| Subread / featureCounts | 2.1.1 |
| Python | 3.10 (pysam 0.24.0, pandas, NumPy, SciPy, matplotlib) |

```bash
mamba install -y -c bioconda   bcftools samtools hisat2 gatk4 subread
mamba install -y -c conda-forge pysam pandas numpy scipy matplotlib
```

HISAT2 (rather than STAR) is used so that the genome index fits in 32 GB RAM.

---

## Inputs

- **Reference:** T2T-CHM13v2.0 FASTA (`chr21`-named, with `.fai`).
- **Annotation:** NCBI RefSeq GFF for T2T-CHM13v2.0 (assembly GCF_009914755.1;
  chr21 contig `NC_060945.1`).
- **WGS BAMs** (CHM13v2-aligned, indexed): the three disomy lines dP, dM1, dM2
  (required); the parental trisomy line (quality control only).
- **RNA-seq FASTQ** (paired-end): the parental trisomy line and the edited
  line(s), n = 3 biological replicates per condition.

Paths are set in `00_config.sh`.

---

## Configuration (`00_config.sh`)

| Variable | Default | Meaning |
|----------|---------|---------|
| `REF` | — | reference FASTA |
| `CHR21` | `chr21` | target contig name |
| `ANNOT_GFF` | — | RefSeq GFF |
| `OUT` | — | output root |
| `THREADS` | 8 | threads |
| `MIN_DP_WGS` | 12 | minimum WGS depth per disomy line |
| `WGS_AD_THR` | 0.05 | WGS allele-depth gate: max fraction of a copy's private allele tolerated in the deletion line |
| `ASE_MIN_TOTAL` | 20 | minimum pooled reads per gene-class (gene-level summaries) |
| `CALL_DROP` | 0.40 | blind calling: a copy is "reduced" if its dosage falls by more than this fraction vs baseline |
| `CENSAT_BED` | (optional) | satellite/rDNA mask |

`set -a` in the config auto-exports these variables to the Python steps.

---

## Pipeline

Run each step after `source 00_config.sh`.

| Step | Role | Key outputs |
|------|------|-------------|
| `01_call_disomy_variants.sh` | Joint diploid genotyping of dP/dM1/dM2 on chr21 (bcftools mpileup/call, FORMAT/AD,DP); parental line genotyped for QC only. | `vcf/trio.chr21.vcf.gz`, `vcf/T21.qc.vcf.gz` |
| `02_classify_alleles.py` | Assign each biallelic SNV to a private copy from the disomy genotypes; drop multiallelic/ambiguous sites and apply the **WGS allele-depth gate** (the deletion line must carry essentially none of the private allele; symmetric across P/M1/M2; RNA-independent). | `snp/informative_sites.tsv`, `snp/hetsites.vcf.gz` |
| `03_build_resolvability_map.py` | Map informative sites to exons (incl. UTR/non-coding) and biotypes; normalize `NC_060945.1`→`chr21`. | `snp/sites_with_gene.tsv`, `snp/gene_resolvability.tsv` |
| `04_align_count.sh` | Manifest-driven HISAT2 alignment (index built once; existing BAMs reused) and GATK ASEReadCounter at the informative sites (MAPQ ≥ 60, BQ ≥ 20). | `ase/<sid>.ase.table` |
| `05_estimate_and_plot.py` | Per-condition allele fractions and a chromosome-wide **dosage summary** with automatic (incl. blind) calling of the affected copy; per-condition class-shift figures; three-copy stacked composition (affected copy always directly measured, at most one retained copy residual-completed); a two-way *affected-vs-rest* figure maximizing gene coverage. | `fig/dosage_summary.tsv`, `fig/class_shift_<cond>.*`, `fig/gene_allele_fractions.tsv`, `fig/fig_stacked_full.pdf`, `fig/fig_affected_vs_rest_<cond>.pdf` |
| `06_diagnose_residual.py` | For each deleted copy, classify any residual allele signal in the edited RNA against the deletion line's WGS (ground truth): `MISMAP_WGS` / `RNA_ONLY` / low-floor. | `fig/residual_diagnosis_<D>.tsv` |
| `07_absolute_expression.py` | Gene-level expression via featureCounts over chr21 exons, normalized by non-chr21 mapped reads, multiplied by directly measured allele fractions to give **absolute per-copy expression**. | `expr/absolute_expression.tsv`, `expr/fig_absolute_<cond>.pdf`, `expr/fig_absolute_stacked.pdf` |
| `08_statistics.py` | **Replicate-level (n = 3) statistics** as the primary unit: per-replicate chromosome-wide allele fraction, total chr21 dosage and absolute per-copy expression, compared between conditions by Welch's t-test (primary) and Mann-Whitney U (reported; underpowered at n = 3). Secondary: gene-level paired Wilcoxon signed-rank. | `stats/replicate_metrics.tsv`, `stats/stats_summary.tsv`, `stats/fig_stats_*` |

---

## Sample manifest (`rna/manifest.tsv`)

Tab-separated, no header, `#` denotes comments. Columns:

```
sid    R1    R2    condition    deleted_copy(none|P|M1|M2|unknown)
```

- `deleted_copy = none` — baseline (parental trisomy)
- `deleted_copy = P|M1|M2` — known target (deletion or inactivation)
- `deleted_copy = unknown` — bulk sample of unknown outcome (the affected copy
  and its residual dosage are inferred from the data)
- `condition` is a free-text label (e.g. `delM2`, `inactM1`, `editX`) used to
  group replicates.

`04_align_count.sh` writes a default manifest on first run; add rows for new
experiments and re-run.

---

## Running the pipeline

Full run:

```bash
source 00_config.sh
bash   01_call_disomy_variants.sh
python 02_classify_alleles.py
python 03_build_resolvability_map.py
bash   04_align_count.sh
python 05_estimate_and_plot.py
python 06_diagnose_residual.py
python 07_absolute_expression.py
python 08_statistics.py
```

**Adding a new experiment** (deletion, inactivation, or unknown bulk): steps
01–03 are experiment-independent (site assignment and the WGS gate depend only
on the disomy WGS and are already fixed). Add rows to `rna/manifest.tsv` and
re-run **04 → 05 → 06 → 07 → 08** only.

**Re-running step 02** (e.g. after changing a threshold) does not require
re-aligning the RNA: the informative-site set is simply narrowed, and steps
05/06 inner-join on position, so only surviving sites are used.

---

## Outputs and interpretation

- `fig/dosage_summary.tsv` — per-condition normalized P:M1:M2 with an automatic
  call (e.g. *"M2 reduced to 0% of baseline"*; partial loss gives an
  intermediate value).
- `stats/stats_summary.tsv` — replicate-level (Welch/MWU) and gene-level
  (Wilcoxon) results. Expected for a clean elimination: the affected copy's
  fraction and absolute expression fall to ~0; the two retained copies' fractions
  rise toward 1/2 while their absolute expression is unchanged; total chr21
  dosage falls to ~2/3.
- `fig/fig_class_shift_<cond>.pdf`, `expr/fig_absolute_<cond>.pdf` — one point
  per gene, baseline vs condition.
- `stats/fig_stats_*` — three replicate points per group with mean ± SD.

---

## Notes and limitations

- **Replication and power.** With n = 3, the Mann-Whitney U test cannot fall
  below p = 0.1; Welch's t-test is the primary test and the gene-level Wilcoxon
  provides orthogonal, well-powered support. Non-significant absolute expression
  for the retained copies indicates *no evidence of change* (consistent with no
  compensation), not proof of equivalence.
- **Resolution ceiling.** Copy-private informative sites are far rarer for the
  maternal homologs (M1, M2) than for P, so the set of genes resolvable into all
  three copies is intrinsically limited; the *affected-vs-rest* view (only one
  informative site needed) covers more genes.
- **Distal mis-mapping.** Subtelomeric chr21 genes can show distorted allelic
  ratios from RNA-specific mis-mapping; the three-copy stack therefore requires
  the affected copy to be directly measured and never residual-fills it.
- **Reference bias.** WASP-style re-mapping can further refine absolute
  fractions; an optional `CENSAT_BED` mask and ONT read-backed phasing QC are
  supported.

---

## Citation

**Software.** R. Hashizume, *T21_ASE_pipeline*, GitHub (2026),
https://github.com/RayHashizume/T21-allele-specific-expression. A versioned, citable archive is
available via Zenodo: https://doi.org/10.5281/zenodo.20772770.

**Associated study.** A manuscript describing this work (R. Hashizume et al.) is
in preparation; this section will be updated with the full citation and DOI upon
publication.

## License

MIT License

## Contact

Ryotaro Hashizume — Department of Pathology and Matrix Biology / Department of
Genomic Medicine, Mie University Graduate School of Medicine.
