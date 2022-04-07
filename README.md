# reanalysis1
Original data and plots uesd for reanalysis purposed


## File descriptions

### complete_mutation_list_with_labels.v3.xlsx

This excel file contains three sheets:

full_list - A complete list of mutations used in the reanalysis studies. Detailed explanations are given in the header and comments.

counts_of_filtered_SNMs - Counts of numbers in different regions. 

gene_analysis - Analysis of 2339 genes of different essentiality, same content as the input file in essentiality_plot/data/gene_level_data.SNMs.csv


### essentiality_plot
Data used and original plots of CDS/INTRON mutation rates in genes of different essentiality


### TSS_TTS_plot
Data used and original plots regarding TSS/TTS regions, detailed statistics are given in "mutation_compare_in_regions_of_different_essentiality.xlsx"


### example_snapshots_of_homomeric_associated_mutations.v2.pptx
Example of BAM alignments shown in IGV (http://www.broadinstitute.org/software/igv/), more snapshots are available in igv_snapshots_for_clustered_somatic_mutations


### Notes:
(1) TEs are annotated based on "TAIR10_GFF3_genes_transposons.gff", so a few more mutations were annotated as TEs than labeled in original Weng's data. But no apparent influence was found in all analysis.

(2) The TSS/TTS positions are likely directly extracted based on the "gene" features provided in "TAIR10_GFF3_genes.gff" though the "gene" features there included 5'- and 3'- UTRs. This was left as it is to keep consistency with original analyses. 

(3) BAM files are generated using BWA-mem (version 0.7.10-r789), and were sort and MarkDuplicates with picard-tools (version 1.114), and further realigned with GenomeAnalysisTK (version 3.7).
