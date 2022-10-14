#********************************************************
#   ggplot-windows-multi.R -- Plot number of each windows along chromosomes using ggplot, support multi-sample input.
#                          
#
#   Author: Nowind
#   Created: 2013-07-10
#   Updated: 2015-11-26
#   Version: 1.0.1
#
#   Change logs:
#   Version 1.0.0 14/11/29: The initial version.
#   Version 1.0.1 15/11/26: Updated: rename some variables.
#*********************************************************


library(lsr)
library(reshape2)
library(ggplot2)
library(grid)
library(lemon)



#*********************************************************
## plot distributions of markers
infile     <- commandArgs(TRUE)[1]
outfile    <- commandArgs(TRUE)[2]

count_data <- as.data.frame(read.table(infile, header=T, sep="\t"))

chrom_num <- length(unique(count_data$chrom))
out_height <- chrom_num * 20

count_data_melted <- melt(count_data, id.vars=c("chrom", "interval", "bin"))

##pdf(file=outfile, width = 240, height = out_height)

ggplot(count_data_melted, aes(x=bin, y=value, colour=variable, group=variable)) +
    geom_line() +
    xlab("Windows(x100kbp)") +
    ylab("Number of Mutations") + 
    scale_colour_hue(name="variable",
                     breaks=unique(count_data_melted$variable),
                     labels=unique(count_data_melted$variable),
                     l=40) +
    scale_x_continuous(breaks = c(0, 50, 100, 150, 200, 250, 300),
                       labels = c(0, 5,  10,  15,  20,  25,  30)) +
    theme_bw() +
    theme(legend.justification=c(1,0), legend.position=c(1,0)) +
    theme(panel.grid = element_blank(),
          axis.line  = element_line(colour = "black")) +
    facet_rep_wrap(~ chrom, ncol=1, scales="free_y", repeat.tick.labels = 'bottom', strip.position="right")
    
ggsave(file=outfile, width=180, height=2*out_height, units="mm", dpi=600, bg='transparent')

#*********************************************************
dev.off()
warnings()

