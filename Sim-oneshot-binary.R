rm(list=ls())

## load functions for HRDD
source("HRDD-Binary.R")

## scenario 
tau_sc <- 1       # 1-3
error_sc <- 1     # 1-3


## settings
set.seed(1)
G <- 100
ng <- rep(c(200, 300, 400, 500), rep(G/4, 4))
ID <- rep(1:G, ng)
N <- sum(ng)


## data generation 
x <- 2*rbeta(N, 2, 4) - 1 

## true treatment effects
Tau.true <- c() 
if(tau_sc==1){
  Tau.true <- rgamma(G, 3, 1) - 3
}
if(tau_sc==2){
  sel <- sample(1:3, G, prob=c(0.4, 0.2, 0.4), replace=T)
  Tau.true <- (-2)*(sel==1) + 2*(sel==3)
}
if(tau_sc==3){
  sel <- sample(1:3, G, prob=c(0.4, 0.2, 0.4), replace=T)
  Tau.true <- runif(G, -3, -1)*(sel==1) + runif(G, 1, 3)*(sel==3)
}

## mean and variance 
Mu1 <- runif(G, 0.4, 1.4)[ID]*x + runif(G, 3, 7)[ID]*x^2 + runif(G, 9, 11)[ID]*x^3  
Mu2 <- Tau.true[ID] + runif(G, 0.4, 1.4)[ID]*x + runif(G, 5, 9)[ID]*x^2 + runif(G, 3, 5)[ID]*x^3  
Sig <- runif(G, 0.5, 1.2)[ID]

## generation of Y
if(error_sc==1){  Ep <- rnorm(N) }
if(error_sc==2){  Ep <- rt(N, 3) }
if(error_sc==3){  Ep <- rgamma(N, 4, 2) - 2 }
Y_ast <- Mu1*ifelse(x<0, 1, 0) + Mu2*ifelse(x>0, 1, 0) + Sig*rnorm(N)
Y <- ifelse(Y_ast>0, 1, 0)
TE <- pnorm(Tau.true) - pnorm(0)    # true treatment effect


## separate RDD
Est <- matrix(NA, G, 3) 
dimnames(Est)[[2]] <- c("sRDD-c", "sRDD-r", "HRDD-N")
CI1 <- CI2 <- CI3 <- matrix(NA, G, 2)
hg <- c()
for(g in 1:G){
  rd <- rdrobust(y=Y[ID==g], x=x[ID==g], c=0)
  CI1[g,] <- rd$ci[1,]
  CI2[g,] <- rd$ci[2,]
  CI3[g,] <- rd$ci[3,]
  Est[g, 1:2] <- rd$coef[1:2]
  hg[g] <- rd$bws[1,1]
}


## HRDD
HRDD <- HRDD_bin(Y=Y, x=x, ID=ID, band="local", h=hg, c=0, order=1, print=T)
Est[,3] <- apply(HRDD$TE, 2, mean)
CI4 <- t( apply(HRDD$TE, 2, quantile, prob=c(0.025, 0.975)) )


## RMSE
sqrt( apply((Est-TE)^2, 2, mean) )


## Figure 
bb <- 0.2
ind <- order(TE)
plot(TE[ind], ylim=range(CI3), ylab="Treatment effect", xlab="Group", pch=20, col=1)
points((1:G)+bb, Est[ind, 3], pch=1, col=2)
points((1:G)-bb, Est[ind, 2], pch=2, col=4)
for(g in 1:G){
  lines(c(g,g)+bb, CI4[ind[g],], col=2)
  lines(c(g,g)-bb, CI3[ind[g],], col=4)
}
abline(h=0)
legend("bottomright", c("True", "HRDD-N", "sRDD-r"), col=c(1,2,4), lty=c(NA, 1, 1), pch=c(20,1,2))











# coverage 
CI <- list(CI1, CI2, CI3, CI4)
for(l in 1:4){
  print( mean(CI[[l]][,1]<TE & CI[[l]][,2]>TE) )
}




# plot
plot(TE, Est[,3], ylim=range(TE, Est))
points(TE, Est[,2], col=2)
abline(0, 1)
for(i in 1:m){
  lines(c(TE[i],TE[i]), CI4[i,])
}

