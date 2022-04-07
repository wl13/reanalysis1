library(data.table)
library(ggplot2)

# gene level data contains epigenomic, functional, and mutation data for coding regions of each gene
#gene_data<-fread("data/gene_level_data.csv")
gene_data<-fread("data/gene_level_data.mut_counts.csv")

# we can examine for example, mutation rates in essential vs non-essential functional categories of genes
#mutation_essentiality<-gene_data[essentiality!="" ,.(mutations=sum((SNV+InDel)*`CDS length`), length=sum(`CDS length`)), by=.(essentiality)]
#mutation_essentiality<-gene_data[essentiality!="" ,.(mutations=sum(`Monroe CDS Count`), length=sum(`Local CDS length`)), by=.(essentiality)]
#mutation_essentiality<-gene_data[essentiality!="" ,.(mutations=sum(`Monroe intron Count`), length=sum(`Local Intron length`)), by=.(essentiality)]
#mutation_essentiality<-gene_data[essentiality!="" ,.(mutations=sum(`Weng CDS Count`), length=sum(`Local CDS length`)), by=.(essentiality)]
#mutation_essentiality<-gene_data[essentiality!="" ,.(mutations=sum(`Weng intron Count`), length=sum(`Local Intron length`)), by=.(essentiality)]
#mutation_essentiality<-gene_data[essentiality!="" ,.(mutations=sum(`Monroe INDEL Count`), length=sum(`Overall length`)), by=.(essentiality)]
mutation_essentiality<-gene_data[essentiality!="" ,.(mutations=sum(`Weng INDEL Count`), length=sum(`Overall length`)), by=.(essentiality)]



ggplot(mutation_essentiality, aes(x=essentiality, y=mutations/length))+
  geom_bar(stat="identity")
  