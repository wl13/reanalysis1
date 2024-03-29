# reanalysis1
Original data and codes uesd for reanalysis purposed

<br />

## File descriptions

<br />

### Fig1a,b, FigS1a~g

These folders contain source data and codes used for producing each figure and supplementary figures


<br />

### igv_snapshots
Example of BAM alignments (against TAIR10 reference) shown in IGV (http://www.broadinstitute.org/software/igv/), the File A2.example_snapshots_of_homomeric_associated_mutations.v2.pptx gives examples of biased strands, more snapshots for problematic clustered "mutations" are given in **clustered_somatic_mutations**

<br />


### equilibrium
Code for the dinucleotide equilibrium analysis

<br />


### Notes:
(1) TEs are annotated based on "TAIR10_GFF3_genes_transposons.gff", so a few more mutations were annotated as TEs than labeled in original Weng's data. But no apparent influence was found in all analysis.

(2) Raw fastq files were downloaded from EBI (PRJNA434660:SRR6862404-SRR6862510). BAM files were generated using BWA-mem (version 0.7.10-r789), and were sorted and processed with MarkDuplicates from picard-tools (version 1.114), and further realigned with GenomeAnalysisTK (version 3.7).

<br />
