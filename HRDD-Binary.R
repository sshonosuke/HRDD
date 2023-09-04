library(MCMCpack)
library(mvtnorm)
library(rdrobust)
library(BayesLogit)

### INPUT
## Y: study variable (binary)
## x: running variable
## ID: grouping information (from 1 to G)
## c: threshold value
## SS: Use spike-and-slab prior (if "True")

HRDD_bin <- function(Y, x, ID, c, band="local", h=NULL, order=1, kernel="Triangular", SS=F, 
                     MC=1500, bn=500, print=F){
  ## preparation 
  # basic information
  G <- max(ID)
  ng <- as.vector( table(ID) )
  
  # bandwidth selection 
  if(is.null(h)){
    h <- rep(NA, G)
    if(band=="local"){
      for(g in 1:G){
        h[g] <- rdrobust(y=Y[ID==g], x=x[ID==g], c=c)$bws[1,1] 
      }
    }
    if(band=="global"){ 
      h <- rep(rdrobust(y=Y, x=x, c=c)$bws[1,1], G)   
    }
  }
  
  logit <- function(x){ 1/(1+exp(-x)) }   # logistic function 
  
  # kernel function
  if(kernel=="Triangular"){
    Kernel <- function(x){
      u <- 1-abs(x)
      return((u>0)*u)
    }
  }
  if(kernel=="Window"){
    Kernel <- function(x){
      return( ifelse(x<1, 1, 0) )
    }
  }

  # covariate matrix and local weight
  z <- x - c
  K <- Kernel(abs(z)/h)    # weight
  W <- 1*(z>0)      # treatment indicator
  Z <- cbind(1, z*(z<0), z*(z>0))   # local linear
  if(order>1){
    for(k in 2:order){
      Z <- cbind(Z, z^k*(z<0), z^k*(z>0))
    }
  }
  
  ## sub-samples
  Y <- Y[K>0]
  W <- W[K>0]
  Z <- Z[K>0,]
  ID <- ID[K>0]
  K <- K[K>0]
  N <- length(Y)
  G <- max(ID)
  ng <- as.vector( table(ID) )
  Kap <- K*(Y-0.5)
  
  ## MCMC 
  p <- dim(Z)[2] + 1   # dimension
  
  # initial values
  om <- rep(1, N)
  Tau <- rep(0, G)
  Beta <- matrix(0, G, p-1)
  Psi <- rep(1, p)    # prior variance of theta (estimated)
  M <- rep(0, p)    # prior mean of theta
  Sel <- rep(0, G)    # representing null (1) and non-null (0) 
  null_prob <- 0.5

  # objects to store posterior samples
  Tau.pos <- matrix(NA, MC, G)
  Beta.pos <- array(NA, c(MC, G, p-1))
  Psi.pos <- matrix(NA, MC, p)
  M.pos <- matrix(NA, MC, p)
  Sel.pos <- matrix(NA, MC, G)
  null_prob.pos <- rep(NA, MC) 
  
  # MCMC iteration 
  for(r in 1:MC){
    # omega (Polya-gamma latent variable)
    mu <- Tau[ID]*W + apply(Z*Beta[ID,], 1, sum)
    om <- rpg.gamma(num=N, h=K, z=mu)
    
    # Tau (spike-slab) and Sel 
    Kap_ast <- Kap - om*apply(Z*Beta[ID,], 1, sum)
    for(g in 1:G){
      sW <- W[ID==g]
      sOm <- om[ID==g]
      
      # Tau (treatment effect)
      ep <- 0.01*sqrt(Psi[1])
      if(Sel[g]==1){
        A <- (sum(sOm*sW) + 1/ep^2)^(-1)
        B <- sum(sW*Kap_ast[ID==g]) 
      }
      if(Sel[g]==0){
        A <- (sum(sOm*sW) + 1/Psi[1])^(-1)
        B <- sum(sW*Kap_ast[ID==g]) + M[1]/Psi[1]
      }
      Tau[g] <- rnorm(1, A*B, sqrt(A))
      
      # Sel (null or non-null)
      if(SS){
        dens1 <- dnorm(Tau[g], 0, ep, log=T)
        dens2 <- dnorm(Tau[g], M[1], sqrt(Psi[1]), log=T)
        odds <- (1-null_prob)/null_prob*exp(dens2-dens1)
        Sel[g] <- ifelse(runif(1)<1/(1+odds), 1, 0)
      }
    }
    Tau.pos[r,] <- Tau
    Sel.pos[r,] <- Sel 
   
    # null_prob
    if(SS){
      null_prob <- rbeta(1, 1+sum(Sel), 1+G-sum(Sel))
      null_prob.pos[r] <- null_prob
    }
    
    # Beta
    Kap_ast <- Kap - om*Tau[ID]*W
    for(g in 1:G){
      sW <- W[ID==g]
      sOm <- om[ID==g]
      sZ <- Z[ID==g,]
      A <- solve( t(sZ)%*%(sOm*sZ) + diag(1/Psi[-1]) )
      B <- as.vector( t(sZ)%*%Kap_ast[ID==g] ) + M[-1]/Psi[-1]
      Beta[g,] <- mvrnorm(1, as.vector(A%*%B), A)
    }
    Beta.pos[r,,] <- Beta
    
    # Psi (prior variance)
    ss <- ifelse(Sel==0, 1, 100)
    Psi[1] <- rinvgamma(1, 1+G/2, 1+sum(ss*(Tau-M[1])^2)/2)
    if(Psi[1]>10){ Psi[1] <- 10 }
    vv <- apply((t(Beta)-M[-1])^2, 1, sum)
    Psi[-1] <- rinvgamma(p-1, 1+G/2, 1+vv/2)
    Psi.pos[r,] <- Psi
    
    # M (prior mean)
    if(sum(Sel==0)>0){
      M[1] <- rnorm(1, mean(Tau), sqrt(Psi[1]/sum(Sel==0)))
    }
    M[-1] <- rnorm(p-1, apply(Beta, 2, mean), sqrt(Psi[-1]/G))
    M.pos[r,] <- M
    
    # print
    if(print & round(10*r/MC)==(10*r/MC)){ 
      print( paste0(100*r/MC, "% completeted") ) 
    }
  }
  
  # omit burn-in samples
  Tau.pos <- Tau.pos[-(1:bn),]
  Beta.pos <- Beta.pos[-(1:bn),,]
  Psi.pos <- Psi.pos[-(1:bn),]
  M.pos <- M.pos[-(1:bn),]
  Sel.pos <- Sel.pos[-(1:bn),]
  null_prob.pos <- null_prob.pos[-(1:bn)]
  
  # treatment effect
  TE.pos <- logit(Beta.pos[,,1]+Tau.pos) - logit(Beta.pos[,,1])
  
  ## Result
  if(SS){
    Result <- list(Tau=Tau.pos, Beta=Beta.pos, Psi=Psi.pos, M=M.pos, TE=TE.pos, band=h,
                   Sel=Sel.pos, null_prob=null_prob.pos)
  }else{
    Result <- list(Tau=Tau.pos, Beta=Beta.pos, Psi=Psi.pos, M=M.pos, TE=TE.pos, band=h)
  }
  return(Result)
}

