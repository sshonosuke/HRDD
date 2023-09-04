library(MCMCpack)
library(mvtnorm)
library(rdrobust)

### INPUT
## Y: study variable (continuous)
## x: running variable
## ID: grouping information (from 1 to G)
## c: threshold value
## SS: Use spike-and-slab prior (if "True")

HRDD_cont <- function(Y, x, ID, c, band="local", h=NULL, order=1, kernel="Window", SS=F,
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
  
  # kernel function
  if(kernel=="Window"){
    Kernel <- function(x){
      u <- 1-abs(x)
      return((u>0)*u)
    }
  }
  if(kernel=="Triangular"){
    Kernel <- function(x){
      return( ifelse(x<1, 1, 0) )
    }
  }
  
  # covariate matrix and local weight
  z <- x - c
  K <- Kernel(abs(z)/h[ID])    # weight
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

  ## MCMC 
  p <- dim(Z)[2] + 1   # dimension
  N_om <- sum(K)
  
  # initial values
  Om <- 1
  Tau <- rep(0, G)
  Beta <- matrix(0, G, p-1)
  Psi <- rep(1, p)    # prior variance of theta (estimated)
  Mu <- rep(0, p)    # prior mean of theta
  Sel <- rep(0, G)    # representing null (1) and non-null (0) 
  null_prob <- 0.5
  U <- U_latent <- rep(1, N)    # latent scale variables (in error term)
  w <- 0.03     # contamination ratio
  Out <- rep(0, N)
  nu <- 2     # degrees of freedom for outlying observations
  
  # objects to store posterior samples
  Tau.pos <- matrix(NA, MC, G)
  Beta.pos <- array(NA, c(MC, G, p-1))
  Psi.pos <- matrix(NA, MC, p)
  Mu.pos <- matrix(NA, MC, p)
  Om.pos <- rep(NA, MC)
  Sel.pos <- matrix(NA, MC, G)
  null_prob.pos <- rep(NA, MC) 
  U.pos <- matrix(NA, MC, N)
  w.pos <- rep(NA, MC)
  
  # MCMC iteration 
  for(r in 1:MC){
    # Tau (spike-slab) and Sel
    dY <- Y - apply(Z*Beta[ID,], 1, sum)
    for(g in 1:G){
      sK <- K[ID==g]
      sW <- W[ID==g]
      sY <- dY[ID==g]
      sU <- U[ID==g]
      
      # Tau (treatment effect)
      ep <- 0.01*sqrt(Psi[1])
      if(Sel[g]==1){
        A <- (Om*sum(sU*sK*sW^2) + 1/ep^2)^(-1)
        B <- Om*sum(sU*sK*sW*sY) 
      }
      if(Sel[g]==0){
        A <- (Om*sum(sU*sK*sW^2) + 1/Psi[1])^(-1)
        B <- Om*sum(sU*sK*sW*sY) + Mu[1]/Psi[1]
      }
      Tau[g] <- rnorm(1, A*B, sqrt(A))
      
      # Sel (null or non-null)
      if(SS){
        dens1 <- dnorm(Tau[g], 0, ep, log=T)
        dens2 <- dnorm(Tau[g], Mu[1], sqrt(Psi[1]), log=T)
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
    dY <- Y - Tau[ID]*W
    for(g in 1:G){
      sZ <- Z[ID==g,]
      sY <- dY[ID==g]
      sK <- K[ID==g]
      sU <- U[ID==g]
      A <- solve( Om*t(sZ)%*%(sU*sK*sZ) + diag(1/Psi[-1]) )
      B <- as.vector( Om*t(sZ)%*%(sU*sK*sY) ) + Mu[-1]/Psi[-1]
      Beta[g,] <- mvrnorm(1, as.vector(A%*%B), A)
    }
    Beta.pos[r,,] <- Beta
    
    # Psi
    Psi[1] <- rinvgamma(1, 1+sum(Sel==0)/2, 1+sum((Tau[Sel==0]-Mu[1])^2)/2)
    vv <- apply((t(Beta)-Mu[-1])^2, 1, sum)
    Psi[-1] <- rinvgamma(p-1, 1+G/2, 1+vv/2)
    Psi.pos[r,] <- Psi
    
    # Mu
    Mu[1] <- rnorm(1, mean(Tau), sqrt(Psi[1]/sum(Sel==0)))
    Mu[-1] <- rnorm(p-1, apply(Beta, 2, mean), sqrt(Psi[-1]/G))
    Mu.pos[r,] <- Mu
    
    # Om
    resid <- Y - Tau[ID]*W - apply(Z*Beta[ID,], 1, sum)
    Om <- rinvgamma(1, 1+N_om/2, 1+sum(U*K*resid^2)/2) 
    Om.pos[r] <- Om
    
    # U & Out
    U_latent <- rgamma(N, nu+Out*K/2, nu+Out*K*Om*resid^2/2)
    L1 <- K*dnorm(resid, 0, 1/sqrt(Om*U_latent), log=T)
    L0 <- K*dnorm(resid, 0, 1/sqrt(Om), log=T)
    prob <- 1/(1+exp(log(1-w)+L0-log(w)-L1))
    Out <- rbinom(N, 1, prob)
    U <- (1-Out) + Out*U_latent
    U.pos[r,] <- U
    
    # w
    w <- rbeta(1, 1+sum(Out), 1+N-sum(Out))
    w.pos[r] <- w
    
    # print
    if(print & round(10*r/MC)==(10*r/MC)){ 
      print( paste0(100*r/MC, "% completeted") ) 
    }
  }
  
  # omit burn-in samples
  Tau.pos <- Tau.pos[-(1:bn),]
  Beta.pos <- Beta.pos[-(1:bn),,]
  Psi.pos <- Psi.pos[-(1:bn),]
  Mu.pos <- Mu.pos[-(1:bn),]
  Om.pos <- Om.pos[-(1:bn)]
  Sel.pos <- Sel.pos[-(1:bn),]
  null_prob.pos <- null_prob.pos[-(1:bn)]
  U.pos <- U.pos[-(1:bn),]
  w.pos <- w.pos[-(1:bn)]
  
  ## Result
  if(SS){
    Result <- list(Tau=Tau.pos, Beta=Beta.pos, Psi=Psi.pos, Mu=Mu.pos, Om=Om.pos,
                   Sel=Sel.pos, null_prob=null_prob.pos, w=w.pos, U=U.pos, band=h)
  }else{
    Result <- list(Tau=Tau.pos, Beta=Beta.pos, Psi=Psi.pos, Mu=Mu.pos, Om=Om.pos,
                   w=w.pos, U=U.pos, band=h)
  }
  return(Result)
}

