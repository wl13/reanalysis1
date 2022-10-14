### plot mutated regions
library(ggplot2)


regions <- read.csv("numbers_in_different_regions.csv", header = TRUE, sep="\t")

cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(regions, aes(fill=Affected_regions, y=SNMs, x=Source)) + 
    geom_bar(position="fill", stat="identity", width = 0.6) +
    scale_fill_manual(values=cbbPalette)

ggsave("numbers_in_different_regions.pdf", width=9.6, height=6.4)

