# Methods (draft)

## Allele-specific quantification of the three chromosome 21 homologs

**Reference and annotation.** All sequencing data were aligned to the
telomere-to-telomere human reference T2T-CHM13v2.0, and gene models were taken
from the corresponding NCBI RefSeq annotation (assembly GCF_009914755.1). Only
chromosome 21 was analysed; the CHM13v2.0 chr21 contig (NCBI accession
NC_060945.1) was renamed `chr21` for consistency across tools.

**Definition of the three chr21 copies and the induced-disomy panel.** The
trisomy 21 iPSC line carries three chr21 homologs, designated P (paternal) and
M1 and M2 (the two maternal homologs). Because the trisomy arose through a
meiosis I non-disjunction event, M1 and M2 are non-identical homologs and are
therefore distinguishable along the whole chromosome. Three isogenic induced-
disomy lines, each retaining a defined pair of homologs after elimination of a
single chr21 copy — dP (M1+M2), dM1 (P+M2) and dM2 (P+M1) — were used to phase
chr21 variants; the single-chromosome composition of each disomy line had been
confirmed by G-banding, FISH and STR analysis.

**Variant calling and copy (homolog) assignment.** Whole-genome sequencing
reads from the three disomy lines were used jointly to call chr21 variants with
bcftools (v1.20) `mpileup`/`call` under a diploid model, emitting per-sample
allele depths (FORMAT/AD) and read depths (FORMAT/DP). At every biallelic SNV,
the genotypes of dP, dM1 and dM2 were used to infer, for each of the three
copies, which allele it carries: a "private" allele of a given copy is the
allele that is absent from the disomy line lacking that copy and present in the
other two. SNVs whose three disomy genotypes were jointly consistent with a
unique copy-to-allele assignment were retained as informative sites; multi-
allelic records, monomorphic sites and sites yielding an ambiguous assignment
were discarded. Sites were required to have a read depth ≥12 in each disomy
line.

**WGS allele-depth gate.** Because the disomy lines provide a genomic ground
truth, each informative site was additionally required to be consistent with the
expected genomic absence of the relevant allele: at a site private to copy C,
the disomy line lacking C (dC) must carry essentially none of the C-private
allele in its WGS reads. Sites at which the deleted allele nonetheless accounted
for >5% of WGS reads in dC were removed, as such sites reflect mis-mapping in
segmentally duplicated or repetitive regions or copy-number artefacts. This
filter is symmetric across the three copies and is independent of the RNA-seq
data. The trisomy WGS was genotyped only for quality control, as its triploid
allele ratios are mis-called by a diploid model.

**RNA-seq alignment and allele-specific read counting.** Paired-end RNA-seq
reads from the parental trisomy line and from the chr21-elimination line(s)
(n = 3 biological replicates per condition) were aligned to T2T-CHM13v2.0 with
HISAT2 (v[VERSION]). Allele-specific read counts at the informative
heterozygous sites were obtained with GATK ASEReadCounter (GATK v4.4.0.0) using
a minimum mapping quality of 60 (uniquely mapped reads) and a minimum base
quality of 20.

**Allele-fraction estimation.** For each gene, reads supporting each copy's
allele were summed across the informative sites assigned to that gene's exons
(including UTRs and non-coding exons) and across biological replicates, and the
expression fraction of a copy was estimated as the copy's allele count divided
by the total count at its sites (assuming a uniform allelic ratio across the
gene). A gene was considered fully resolvable when all three copies were
directly measured by their own informative sites; fractions were never inferred
by residual subtraction, so a copy lacking a directly informative site was not
assigned a value for that gene. Chromosome-wide per-copy fractions were computed
analogously by pooling all informative sites of each copy class.

**Absolute per-copy expression.** Gene-level expression was quantified with
featureCounts (Subread v2.1.1) over chr21 exons in paired-end fragment-counting
mode and normalized to the number of mapped reads on chromosomes other than
chr21 (reads per million non-chr21 reads), a denominator that is robust to the
chr21 dosage change itself. The absolute expression attributable to each copy
was obtained as the normalized chr21 (or gene-level) expression multiplied by
the directly measured allele fraction of that copy.

**Statistical analysis.** Biological replicates (n = 3 per condition) were the
unit of replication. For each replicate, the chromosome-wide allele fraction of
each copy, the total chr21 expression dosage, and the absolute per-copy
expression were computed as above. The parental and elimination lines were
compared as independent groups using Welch's two-sample t-test (primary test);
the Mann–Whitney U test is also reported, but with three replicates per group
its minimum two-sided p value is 0.1 and it is therefore underpowered. As an
orthogonal, gene-level test of a consistent shift across chr21 genes, per-gene
allele fractions were compared between conditions with the Wilcoxon signed-rank
test (genes as the unit). Data were analysed in Python (pandas, NumPy, SciPy);
analysis scripts are available at https://github.com/RayHashizume/T21-allele-specific-expression (DOI: 10.5281/zenodo.20772770).

**Note on mechanism.** RNA-seq quantifies expressed allelic dosage and does not
distinguish physical chromosome elimination from transcriptional silencing of an
intact homolog, as both abolish the corresponding allele's transcripts;
mechanism was established by orthogonal DNA/cytogenetic assays.
