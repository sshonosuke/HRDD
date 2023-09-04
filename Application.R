rm(list=ls())
library(rdrobust)

# load dataset
set.seed(1)
load("Dataset.RData")

Y <- Data[,1]   # response
X <- Data[,2]   # covariate (running variable)
G <- max(ID)    # number of subgroups



## separate RDD (local bandwidth)
Est <- matrix(NA, G, 4) 
CI1 <- CI2 <- CI3 <- matrix(NA, G, 2)

hg <- c()     # vector of subgroup-wise bandwidth
for(g in 1:G){
  try({
    rd <- rdrobust(y=Y[ID==g], x=X[ID==g], c=0)
    CI1[g,] <- rd$ci[1,]
    CI2[g,] <- rd$ci[2,]
    CI3[g,] <- rd$ci[3,]
    Est[g, 1:2] <- rd$coef[1:2]
    hg[g] <- rd$bws
  })
}

hg[is.na(hg)] <- mean(na.omit(hg))    # imputation of subgroups that "rdrobust" cannot apply


## HRDD with normal prior
HRDD <- HRDD_bin(Y=Y, x=X, ID=ID, c=0, band="local", h=hg, order=1, SS=F, print=T, MC=2000, bn=500)
Est[,3] <- apply(HRDD$TE, 2, mean)
CI4 <- t( apply(HRDD$TE, 2, quantile, prob=c(0.025, 0.975)) )




## Figures of results
ng <- as.vector(table(ID))

# fig1
par(mfcol=c(1,2))
plot(Est[,3], Est[,1], xlab="HRDD-N", ylab="sRDD-c")
abline(0, 1)
plot(ng, Est[,1]-Est[,3], ylab="Difference", xlab="Sample size")
abline(0, 0)


# fig2
par(mfcol=c(1,1))
id <- order(Est[,3])
plot(Est[id, 3], ylim=range(0, CI4), pch=20, xlab="Group index", ylab="Treatment effect")
for(g in 1:G){
  lines(c(g,g), CI4[id[g],])
}


