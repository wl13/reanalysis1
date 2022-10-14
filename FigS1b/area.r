euclidean <- function(a, b) sqrt(sum((a - b)^2))


data1 <- read.csv("area.csv", header=TRUE)


WC <- data1[,2]
MC <- data1[,3]
WL <- data1[,4]
ML <- data1[,5]
Pol <- data1[,6]

W_full <- WC + WL
M_full <- MC + ML



WC_tot <- sum(WC)
MC_tot <- sum(MC)
WL_tot <- sum(WL)
ML_tot <- sum(ML)
Pol_tot <- sum(Pol)

W_tot <- sum(W_full)
M_tot <- sum(M_full)

WC_f <- WC/WC_tot
MC_f <- MC/MC_tot
WL_f <- WL/WL_tot
ML_f <- ML/ML_tot
Pol_f <- Pol/Pol_tot
W_f <- W_full/W_tot
M_f <- M_full/M_tot



df <- data.frame(WC_f, MC_f,WL_f, ML_f, W_f, M_f, Pol_f)
cols <- c("Weng_C", "Monroe_C", "Weng_L", "Monroe_L", "Weng_all","Monroe_all", "Pol")

dists <- c()
comp <- c()

for (i in 1:7) {
	for (j in 1:7) {
   
  
 ed <-  euclidean(df[i], df[j])
 
 dists <- c(dists, ed)
 comp <- c(comp, c(i,j))
   
   	}		
}


dists.mat <- matrix(dists, nrow=7, ncol=7, byrow=TRUE)

colnames(dists.mat) <- cols
rownames(dists.mat) <- cols


write.table(format(dists.mat, digits = 4), file = "FigS1b_top_euclidean.csv", sep = ",", row.names = TRUE, col.names = TRUE)

df <- data.frame(WC, MC,WL, ML, W_full, M_full, Pol)
cols <- c("Weng_C", "Monroe_C", "Weng_L", "Monroe_L", "Weng_all","Monroe_all", "Pol")
df.m <- data.matrix(df)

chi.P <- c()
chi.C <- c()
for (i in 1:7) {
	for (j in 1:7) {
   v1 <- df.m[,i]
   v2 <- df.m[,j]
  prematrix <- c(v1, v2)
  m <- matrix(prematrix, ncol=2, nrow=8)
 chi <-  chisq.test(m)
 
chi.P <- c(chi.P, chi$p.value)
 chi.C <- c(chi.C, chi$statistic)
   
   	}		
}


chiP.mat <- matrix(chi.P, nrow=7, ncol=7, byrow=TRUE)
chiC.mat <- matrix(chi.C, nrow=7, ncol=7, byrow=TRUE)

colnames(chiP.mat) <- cols
rownames(chiP.mat) <- cols
colnames(chiC.mat) <- cols
rownames(chiC.mat) <- cols


write.table(format(chiP.mat, digits = 4), file = "FigS1b_low_chi_P_area.csv", sep = ",", row.names = TRUE, col.names = TRUE)
write.table(format(chiC.mat, digits = 4), file = "chi_C_area.csv", sep = ",", row.names = TRUE, col.names = TRUE)
