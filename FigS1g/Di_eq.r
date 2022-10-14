library(stringr)
library(ggrepel)
library(ggplot2)
library(gridExtra)


data1 <- read.csv("Dinuc_equil.csv", header = TRUE)

attach(data1)
dinucs <- data1$dinucleotide
M_e <- data1[,2]
W_e <- data1[,3]
intron_abs <- data1$intron
interg_abs <- data1$intergenic
geno_abs <- data1$genomic


intron_rel <- intron_abs/(sum(intron_abs))
interg_rel <- interg_abs/(sum(interg_abs))
geno_rel <- geno_abs/(sum(geno_abs))

data1$interg_rel <- interg_rel

dinucs.1 <- dinucs
dinucs.1 <- str_replace_all(dinucs.1, "AA", "blue")

dinucs.1 <- str_replace_all(dinucs.1, "TT", "blue")
dinucs.1 <- str_replace_all(dinucs.1, "[ACTG][ACTG]", "black")

data1$dinucs.1 <- dinucs.1

M_e_notAATT <- M_e[c(2:15)]
W_e_notAATT <- W_e[c(2:15)]
interg_rel_notAATT <- interg_rel[c(2:15)]

pdf("Fig3.pdf")


par(mfrow=c(1,2))
par(pty="s")
plot(W_e, interg_rel,  type = "n", xlab = "Equil. frequency - Weng et al", ylab = "Intergenic observed frequency", cex.axis=0.8, cex.lab=0.8)
text(W_e, interg_rel,labels=dinucs, col=dinucs.1, cex= 0.6, offset = 10)

abline(lm(interg_rel_notAATT~W_e_notAATT))


plot(M_e, interg_rel,  type = "n", xlab = "Equil. frequency - Monroe et al", ylab = "Intergenic observed frequency", cex.axis=0.8,cex.lab=0.8)
text(M_e, interg_rel,labels=dinucs, col=dinucs.1, cex= 0.5, offset = 10)

abline(lm(interg_rel_notAATT~M_e_notAATT))


dev.off()

M_e_AA <- M_e[1]
M_e_TT <- M_e[16]
M_e_AATT <- M_e_AA + M_e_TT

W_e_AA <- W_e[1]
W_e_TT <- W_e[16]
W_e_AATT <-W_e_AA + W_e_TT



interg_rel_AA <- interg_rel[1]
interg_rel_TT <- interg_rel[16]
interg_rel_AATT <-interg_rel_AA + interg_rel_TT 

interg_abs_AA <- interg_abs[1]
interg_abs_TT <- interg_abs[16]
interg_abs_AATT <-interg_abs_AA + interg_abs_TT 

Num_E_M <- sum(interg_abs)*M_e_AATT
Num_E_W <- sum(interg_abs)*W_e_AATT
Num_O_Interg <- interg_abs_AATT

ch.M <- (Num_O_Interg - Num_E_M)^2/Num_E_M
ch.W <- (Num_O_Interg - Num_E_W)^2/Num_E_W

lm_fit_M_notAATT <- lm(interg_rel_notAATT~M_e_notAATT)
sm.M_notAATT <- summary(lm_fit_M_notAATT)
lm_fit_W_notAATT <- lm(interg_rel_notAATT~W_e_notAATT)
sm.W_notAATT <- summary(lm_fit_W_notAATT)  
lm_fit_M_all <- lm(interg_rel~M_e)
sm.M_all <- summary(lm_fit_M_all)
lm_fit_W_all <- lm(interg_rel~W_e)
sm.W_all <- summary(lm_fit_W_all) 


#try plot with ggrepel
options(ggrepel.max.overlaps = Inf)
p.M <- ggplot(data1, aes(Monroe_estimate,interg_rel, label = dinucleotide)) +
  geom_point(color = dinucs.1)+ theme(aspect.ratio=1) + geom_abline(slope=lm_fit_M_notAATT$coef[2], intercept =lm_fit_M_notAATT$coef[1]) + geom_abline(slope=lm_fit_M_all$coef[2], intercept =lm_fit_M_all$coef[1], linetype=2)
  p.M.1 <- p.M + geom_text_repel(colour=dinucs.1) +labs(x = "Monroe equilibrium freq", y="Intergenic frequency")
  
  p.W <- ggplot(data1, aes(Weng_estimate,interg_rel, label = dinucleotide))+
  geom_point(color = dinucs.1)+ theme(aspect.ratio=1) + geom_abline(slope=lm_fit_W_notAATT$coef[2], intercept =lm_fit_W_notAATT$coef[1])+ geom_abline(slope=lm_fit_W_all$coef[2], intercept =lm_fit_W_all$coef[1], linetype=2)
  p.W.1 <- p.W + geom_text_repel(colour=dinucs.1) +labs(x = "Weng Equilibrium freq", y="Intergenic frequency")
  
 pdf(file = "Figs1g.pdf", height = 6, width = 12)
 
  gridExtra::grid.arrange(p.W.1, p.M.1, ncol = 2)
  
  dev.off()
  
  #check for equal slopes
  
lm_fit_M_notAATT <- lm(interg_rel_notAATT~M_e_notAATT)
sm.M_notAATT <- summary(lm_fit_M_notAATT)
lm_fit_W_notAATT <- lm(interg_rel_notAATT~W_e_notAATT)
sm.W_notAATT <- summary(lm_fit_W_notAATT)  
lm_fit_M_all <- lm(interg_rel~M_e)
sm.M_all <- summary(lm_fit_M_all)
lm_fit_W_all <- lm(interg_rel~W_e)
sm.W_all <- summary(lm_fit_W_all)  


coef.M_notAATT <- sm.M_notAATT$coef[,2]
coef.W_notAATT <- sm.W_notAATT$coef[,2]
coef.M_all <- sm.M_all$coef[,2]
coef.W_all <- sm.W_all$coef[,2]

slope.M.nAATT <- sm.M_notAATT$coef[2]
slope.W.nAATT <- sm.W_notAATT$coef[2]
slope.M.all <- sm.M_all$coef[2]
slope.W.all <- sm.W_all$coef[2]


sem.M.nAATT <- sm.M_notAATT$coef[4]
sem.W.nAATT <- sm.W_notAATT$coef[4]
sem.M.all <- sm.M_all$coef[4]
sem.W.all <- sm.W_all$coef[4]


t.M <- (slope.M.nAATT - slope.M.all)/(sqrt(sem.M.nAATT^2 + sem.M.all^2))
t.W <- (slope.W.nAATT - slope.W.all)/(sqrt(sem.W.nAATT^2 + sem.W.all^2))

p.value.M <- dt(t.M, df=26)
  p.value.W <- dt(t.W, df=26)
