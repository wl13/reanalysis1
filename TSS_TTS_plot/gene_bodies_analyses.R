library(data.table)

#####

#variants<-fread("data/2000000_random_positions_in_AT_homomeric_regions.csv")
#variants<-fread("data/6320_somatic_muts.csv")
#variants<-fread("data/1865_germline_SNMs.csv")
#variants<-fread("data/1021_germline_nonTE_SNMs.csv")
#variants<-fread("data/760_original_Weng_version_of_unique_TE_SNMs.csv")
#variants<-fread("data/844_germline_TE_SNMs.csv")
variants<-fread("data/841_Weng_TE_SNMs.csv")

# annotaiton file source and cleanup
#download.file("https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff", destfile="data/TAIR10_GFF3_genes.gff")
annotation<-fread("data/TAIR10_GFF3_genes.gff")
setnames(annotation, c("V1","V3","V4","V5"),c("chr", "type", "start", "stop"))
annotation$chr<-gsub("Chr","",annotation$chr)

# polymorphology contains functions to visualize variation in relation to TSS and TTS of gene bodies
# install if needed:
#library(devtools)
#install_github("greymonroe/polymorphology")
library(polymorphology)
tss_view<-tss_tts.variants(gff = annotation, vcf=variants)

## for 2000000 random positions
#tss_tts.variants.plot(tss_view, window=10)

## for germline or somatic mutations
tss_tts.variants.plot(tss_view, window=200)



