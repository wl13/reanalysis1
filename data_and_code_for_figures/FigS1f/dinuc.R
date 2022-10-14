library(stringr)
library(ggrepel)
library(ggplot2)
library(gridExtra)

data1 <- read.csv("dinuc.csv", header = TRUE)
attach(data1)

W_all <- data1$Weng_all
M_all <- data1$Monroe_all
cnts <- data1[,4]
dinucs <- data1[,1]


W_norm <- W_all/cnts
M_norm <- M_all/cnts

W_tot <- sum(W_norm)
M_tot <- sum(M_norm)


W_rel <- W_norm/W_tot
M_rel <- M_norm/M_tot

data1$W_rel <- W_rel
data1$M_rel <- M_rel


dinucs.1 <- dinucs
dinucs.1 <- str_replace_all(dinucs.1, "[GC]A>AA", "blue")
dinucs.1 <- str_replace_all(dinucs.1, "A[GC]>AA", "blue")
dinucs.1 <- str_replace_all(dinucs.1, "[GC]T>TT", "blue")
dinucs.1 <- str_replace_all(dinucs.1, "T[GC]>TT", "blue")

dinucs.1 <- str_replace_all(dinucs.1, "TA>AA", "darkred")
dinucs.1 <- str_replace_all(dinucs.1, "AT>AA", "darkred")
dinucs.1 <- str_replace_all(dinucs.1, "AT>TT", "darkred")
dinucs.1 <- str_replace_all(dinucs.1, "TA>TT", "darkred")
dinucs.1 <- str_replace_all(dinucs.1, "CG>[CGAT][CGAT]", "darkgreen")
dinucs.1 <- str_replace_all(dinucs.1, ".+>.+", "darkgrey")

data1$dinucs.1 <- dinucs.1

#xlim = c(0,0.065), ylim = c(0, 0.065)

plot(log10(W_rel), log10(M_rel),  type = "n", xlab = "Log(relative normalized frequency - Weng et al)", ylab = "Log(relative normalized frequency - Monroe et al)")
text(log10(W_rel), log10(M_rel),labels=dinucs, col=dinucs.1, cex= 0.8, offset = 10)
pdf("Fig2.pdf")
plot(W_rel, M_rel,  type = "n", xlab = "Relative normalized frequency - Weng et al", ylab = "Relative normalized frequency - Monroe et al", xlim = c(0,0.065), ylim = c(0, 0.065))
text(W_rel, M_rel,labels=dinucs, col=dinucs.1, cex= 0.8, offset = 10)

abline(0,1)

dev.off()


#try plot with ggrepel
options(ggrepel.max.overlaps = 10)
p.M <- ggplot(data1, aes(W_rel,M_rel, label = Di_nucl_changes)) +
  geom_point(color = dinucs.1)+ theme(aspect.ratio=1) + geom_abline(slope=1, intercept =0) + xlim(0, 0.065) + ylim(0, 0.05)
  p.M.1 <- p.M + geom_text_repel(colour=dinucs.1) +labs(x = "Relative normalised freq - Weng", y="Relative normalised freq - Monroe")
  
  
  
 pdf(file = "FigED1f.pdf", height = 6, width = 6)
 
 gridExtra::grid.arrange(p.M.1, ncol = 1)
  
dev.off()

#chi squared on raw numbers

df <- data.frame(W_all, M_all)
cols <- c("Weng_all", "Monroe_all")
df.m <- data.matrix(df)

chi.P <- c()
chi.C <- c()

   v1 <- df.m[,1]
   v2 <- df.m[,2]
  prematrix <- c(v1, v2)
  m <- matrix(prematrix, ncol=2, nrow=96)
 chi <-  chisq.test(m)
 
chi.P <- c(chi.P, chi$p.value)
 chi.C <- c(chi.C, chi$statistic)
   


  