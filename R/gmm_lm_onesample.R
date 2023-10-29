### Multivariable GMM and Kleibergen's Lagrange multplier (LM) statistics
### with one-sample individual-level data

## by Haodong Tian, Ashish Patel, Stephen Burgess

## Inputs
# Z = n x J instrument matrix
# X = n x K exposure matrix
# Y = n x 1 outcome vector
# beta0 = the tested null of the causal parameter value (for LM test only)

## Output
# gmm_est = K vector of causal effect estimates using GMM
# gmm_se = K vector of standard errors corresponding to gmm_est  (没有covariance matrix of GMM estimators吗？)
# gmm_pval = K vector of p-values corresponding to gmm_est
# Q_stat = overidentification test statistic (即，IV pleiotropy test; or IV validity test)
# Q_pval = overidentification test p-value
# lm_pval = p-value of the LM test of the null hypothesis H0: beta=beta0

gmm_lm_onesample <- function(X,Y,Z,beta0=NA){
  n = nrow(Z); J = ncol(Z); K = ncol(X)

  # mean-centered variables  (Z X Y都centred了)
  Y <- Y-mean(Y)
  X0 <- function(k){X[,k]-mean(X[,k])}; X <- sapply(1:K,X0); rm(X0)
  Z0 <- function(k){Z[,k]-mean(Z[,k])}; Z <- sapply(1:J,Z0); rm(Z0)

  # GMM functions----------------------------------------------------------------------------------------
  g <- function(bet){as.vector(t(Z)%*%(Y-(X%*%bet)))/n} #moment condition vector: checked
  Om <- function(bet){(t(Z)%*%Z/n)*as.numeric(t(Y-(X%*%bet))%*%(Y-(X%*%bet))/n)} #Omega matrix: checked; the estimated variance matrix with homoskedasticity assumption
  G <- -(t(Z)%*%X)/n #checked , the estimated G vector #记得有个负号
  Q <- function(bet){as.numeric(t(g(bet))%*%solve(Om(bet))%*%g(bet))}  #Q statistic; also the argmin -> GMM estimator

  # preliminary estimate (主要是给一个initial value for optimization)
  Q.gg <- function(bet){as.numeric(t(g(bet))%*%g(bet))}  #weighting matrix == identify matrix I
  bet.gg <- nlminb(rep(0,K),objective=Q.gg)$par  #optimization
  #这种preliminary设置其实和我的iterative Q的preliminary设置有点像
  #即,先把weight设置成最naive的情况，然后起始点为c(0,0),nlminb得到的极值点其实就是Weighted Linear Regression的解

  # GMM estimate
  DQ <- function(bet){2*as.matrix(t(G)%*%solve(Om(bet))%*%g(bet))}#gradient辅助函数；其实就是Q关于beta的导数

  gmm <- nlminb(bet.gg,objective=Q,gradient=DQ)$par #optimization

  var.gmm <- as.matrix(solve(t(G)%*%solve(Om(gmm))%*%G))#estimated variance matrix of the GMM estimator

  #Q statistics---------
  Qstat <- Q(gmm)#Q statistics - used for later over-idenitifcation testing (i.e. IV pleiotropy test)
  Q.pval <- pchisq(n*Qstat, df = J-K, lower.tail = F) #p value; df=J-K #checked

  gmm.pval <- 2*pnorm(  -abs(  gmm/ sqrt(    diag(var.gmm)/n   )   )    ) #check; the margin parameter z-test

  LM<-NA; LM.pval<-NA
  if(!is.na(as.matrix(beta0)[1,1] )  ){
    # LM statistic--------------------------------------------------------------------------------------------
    g0 <- g(beta0)
    Om0 <- Om(beta0)
    evec0 <- eigen(Om0)$vectors; eval0 <- eigen(Om0)$values #eigenvectors and eigenvalues
    sqrt.Om0 <- evec0%*%diag(sqrt(eval0))%*%t(evec0) #Diagonalization #方便算matrix^{-1/2}呀
    delG <- list()#Delta matrix (Delta_1 Delta_2 ... Delta_K) where Delta_k is a J*J matrix
    for (k in 1:K){
      delG[[k]] <- matrix(NA,J,J)
      for (j1 in 1:J){
        for (j2 in 1:J){
          delG[[k]][j1,j2] <- mean(-X[,k]*Z[,j1]*Z[,j2]*(Y-X%*%beta0))#X[,k]*Z[,j1]*Z[,j2]*(Y-X%*%beta0) 没问题
        }
      }
    }

    D <- matrix(NA,J,K)
    for (k in 1:K){
      D[,k] <- G[,k] - (delG[[k]]%*%solve(Om0)%*%g0)  #D_k sub-matrix
    }
    P0 = solve(sqrt.Om0)%*%D #checked
    U0 = solve(sqrt.Om0)%*%g0 #checked
    P = P0%*%solve(t(P0)%*%P0)%*%t(P0) #checked
    LM = as.numeric(t(U0)%*%P%*%U0)*n #checked
    LM.pval <- 1-pchisq(LM,K)
  }

  res.list <- list("gmm_est"=gmm,
                   "gmm_se"=sqrt(diag(var.gmm)/n),
                   'variance_matrix'=var.gmm/n,
                   "gmm_pval"=2*pnorm(-abs(gmm/sqrt(diag(var.gmm)/n))),
                   "Q_stat"=Qstat, "Q_pval"=Q.pval,
                   "lm_stat"=LM,
                   "lm_pval"=LM.pval)
  return(res.list)
}



getLM<-function(X,Y,Z,beta0){
  n = nrow(Z); J = ncol(Z); K = ncol(X)

  # mean-centered variables  (Z X Y都centred了)
  Y <- Y-mean(Y)
  X0 <- function(k){X[,k]-mean(X[,k])}; X <- sapply(1:K,X0); rm(X0)
  Z0 <- function(k){Z[,k]-mean(Z[,k])}; Z <- sapply(1:J,Z0); rm(Z0)

  g <- function(bet){as.vector(t(Z)%*%(Y-(X%*%bet)))/n} #moment condition vector: checked
  Om <- function(bet){(t(Z)%*%Z/n)*as.numeric(t(Y-(X%*%bet))%*%(Y-(X%*%bet))/n)} #Omega matrix: checked; the estimated variance matrix with homoskedasticity assumption
  G <- -(t(Z)%*%X)/n #checked , the estimated G vector #记得有个负号

  # LM statistic--------------------------------------------------------------------------------------------
  g0 <- g(beta0)
  Om0 <- Om(beta0)
  evec0 <- eigen(Om0)$vectors; eval0 <- eigen(Om0)$values #eigenvectors and eigenvalues
  sqrt.Om0 <- evec0%*%diag(sqrt(eval0))%*%t(evec0) #Diagonalization #方便算matrix^{-1/2}呀

  #Step1: 算Delta matrix
  delG <- list()#Delta matrix (Delta_1 Delta_2 ... Delta_K) where Delta_k is a J*J matrix
  for (k in 1:K){
    delG[[k]] <- matrix(NA,J,J)
    for (j1 in 1:J){
      for (j2 in 1:J){
        delG[[k]][j1,j2] <- mean(-X[,k]*Z[,j1]*Z[,j2]*(Y-X%*%beta0))#X[,k]*Z[,j1]*Z[,j2]*(Y-X%*%beta0) 没问题
      }
    }
  }

  #Step2: 算D
  D <- matrix(NA,J,K)
  for (k in 1:K){
    D[,k] <- G[,k] - (delG[[k]]%*%solve(Om0)%*%g0)  #D_k sub-matrix
  }

  #Step3：算LM statistic
  P0 = solve(sqrt.Om0)%*%D #checked
  U0 = solve(sqrt.Om0)%*%g0 #checked

  P = P0%*%solve(t(P0)%*%P0)%*%t(P0) #checked

  LM = as.numeric(t(U0)%*%P%*%U0)*n #checked

  LM.pval <- 1-pchisq(LM,K)

  res.list <- list("lm_stat"=LM,
                   "lm_pval"=LM.pval)
  return(res.list)
}

# ###examples
# set.seed(1123)
# N<-10000 #sample size
# J<-30 #30 instruments
# K<-3 #3 exposures
# Z<-matrix(  rbinom(N*J,2,0.3), N,J)  #three-valued SNPs
# U<-matrix(  rnorm(N*K,0,1), N,K )
#
# Alphas<-matrix( runif(J*K  , 0, 0.5 ) , J, K   )  #Z-X effects matrix
# X<- Z%*%Alphas + U #exposure matrix
# Beta<-runif(K,0.5,1) #causal effects
# Beta
# Y<- X%*%Beta + apply( U, 1, sum )
#
# gmm_res<-gmm_lm_onesample( X=X,Y=Y,Z=Z,beta0=NA)
# gmm_res$gmm_est;gmm_res$gmm_se
#
#
# gmm_res<-gmm_lm_onesample( X=X,Y=Y,Z=Z,beta0=c(1,2,3) )


#
# mm_res$variance_matrix
# gmm_res$gmm_est
# gmm_res_<-gmm_lm_onesample_( X=X,Y=Y,Z=Z,beta0=c(1,2,3) )
# gmm_res_$gmm_est
#
#
# LMres<-getLM( X=X,Y=Y,Z=Z,beta0=c(1,2,3) )

