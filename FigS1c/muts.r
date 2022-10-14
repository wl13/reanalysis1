euclidean <- function(a, b) sqrt(sum((a - b)^2))


data1 <- read.csv("Muts.csv", header=TRUE)


WC <- data1[,2]
MC <- data1[,5]
WL <- data1[,4]
ML <- data1[,3]
#sites <- data1[,6]



W_full <- WC + WL
M_full <- MC + ML



WC_tot <- sum(WC)
MC_tot <- sum(MC)
WL_tot <- sum(WL)
ML_tot <- sum(ML)


W_tot <- sum(W_full)
M_tot <- sum(M_full)

WC_f <- WC/WC_tot
MC_f <- MC/MC_tot
WL_f <- WL/WL_tot
ML_f <- ML/ML_tot

W_f <- W_full/W_tot
M_f <- M_full/M_tot



df <- data.frame(WC_f, MC_f,WL_f, ML_f, W_f, M_f)
cols <- c("Weng_C", "Monroe_C", "Weng_L", "Monroe_L", "Weng_all","Monroe_all")

dists <- c()
comp <- c()

for (i in 1:6) {
	for (j in 1:6) {
   
  
 ed <-  euclidean(df[i], df[j])
 
 dists <- c(dists, ed)
 comp <- c(comp, c(i,j))
   
   	}		
}


dists.mat <- matrix(dists, nrow=6, ncol=6, byrow=TRUE)

colnames(dists.mat) <- cols
rownames(dists.mat) <- cols


write.table(format(dists.mat, digits = 4), file = "FigS1c_top_euclidean_muts.csv", sep = ",", row.names = TRUE, col.names = TRUE)




WC_ML <- c(WC, ML)
WCML_mat <- matrix(WC_ML, ncol=2, nrow=12)
chi.WC_ML <- chisq.test(WCML_mat)

df <- data.frame(WC, MC,WL, ML, W_full, M_full)
df.m <- data.matrix(df)
cols <- c("Weng_C", "Monroe_C", "Weng_L", "Monroe_L", "Weng_all","Monroe_all")
chi.P <- c()
chi.C <- c()
for (i in 1:6) {
	for (j in 1:6) {
   v1 <- df.m[,i]
   v2 <- df.m[,j]
  prematrix <- c(v1, v2)
  m <- matrix(prematrix, ncol=2, nrow=12)
 chi <-  chisq.test(m)
 
chi.P <- c(chi.P, chi$p.value)
 chi.C <- c(chi.C, chi$statistic)
   
   	}		
}


chiP.mat <- matrix(chi.P, nrow=6, ncol=6, byrow=TRUE)
chiC.mat <- matrix(chi.C, nrow=6, ncol=6, byrow=TRUE)

colnames(chiP.mat) <- cols
rownames(chiP.mat) <- cols
colnames(chiC.mat) <- cols
rownames(chiC.mat) <- cols


write.table(format(chiP.mat, digits = 4), file = "FigS1c_low_chi_P_muts.csv", sep = ",", row.names = TRUE, col.names = TRUE)
write.table(format(chiC.mat, digits = 4), file = "chi_C_muts.csv", sep = ",", row.names = TRUE, col.names = TRUE)
