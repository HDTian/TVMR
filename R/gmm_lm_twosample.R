### Multivariable GMM and Kleibergen's Lagrange multplier (LM) statistics
### for two-sample summary data Mendelian randomization

## by Haodong Tian, Ashish Patel, Stephen Burgess

## Inputs
# bx = J x K matrix of genetic variant associations with K exposures
# sx = J x K matrix of corresponding standard errors from variant-exposure regressions
# by = J vector of genetic variant associations with outcome
# sy = J vector of corresponding standard errors from variant-outcome regressions
# nx = sample size for genetic variant-exposure regressions
# ny = sample size for genetic variant-outcome regressions
# ld = J x J genetic variant correlation matrix   #很重要的matrix；传统的IVW几乎都要保证SNP之间的covariance = 0才能让likelihood简化
# cor.x = K x K matrix exposure correlation matrix  #即，phenotype之间来自于同一sample的情况
# beta0 = the tested null of the causal parameter value (for LM test only)

## Output
# condFstat = K vector of conditional F statistics (只和G X data有关；即one-sample)
# gmm_est = K vector of causal effect estimates using GMM
# gmm_se = K vector of standard errors corresponding to gmm_est
# gmm_pval = K vector of p-values corresponding to gmm_est
# Q_stat = overidentification test statistic
# Q_pval = overidentification test p-value
# lm_pval = p-value of the LM test of the null hypothesis H0: beta=beta0


gmm_lm_twosample <- function(bx,by,sx,sy,nx,ny,ld,cor.x,beta0=NA){
  J = nrow(bx); K = ncol(bx)

  # exposures quantities of interest
  #(即用univariate regression的summary statistics转变成multi regression的summary statistics)
  ax <- matrix(NA,nrow=J,ncol=K); Ax <- list()
  for (k in 1:K){ax[,k] <- 1/((nx*sx[,k]^2)+bx[,k]^2)}
  for (k in 1:K){Ax[[k]] <- (sqrt(ax[,k])%*%t(sqrt(ax[,k])))*ld}
  Bx <- function(k){ax[,k]*bx[,k]}; Bx <- sapply(1:K,Bx)

  sqrt.Ax <- function(k){                                      #i.e. 对角化然后再求^{1/2}
    evec <- eigen(Ax[[k]])$vectors; eval <- eigen(Ax[[k]])$values
    return((evec%*%diag(sqrt(eval))%*%t(evec)))
  }
  sqrt.Ax <- lapply(1:K,sqrt.Ax)

  SigX <- function(k,l){  #这些SigX的操作主要是为了算论文中的V这个stacking matrix
    solve(sqrt.Ax[[k]]%*%t(sqrt.Ax[[l]]))*(   cor.x[k,l]-(  as.numeric(t(Bx[,k])%*%solve(sqrt.Ax[[k]]%*%t(sqrt.Ax[[l]]))%*%Bx[,l])  )   )*(1/(nx-J+1))
  }
  SigX2 <- function(m){
    SigX2a <- list()
    for (m1 in 1:K){SigX2a[[m1]] <- SigX(m1,m)}
    SigX2a <- do.call(rbind, SigX2a)
    return(SigX2a)
  }
  SigX3 <- list()
  for (m1 in 1:K){SigX3[[m1]] <- SigX2(m1)}
  SigX3 <- do.call(cbind, SigX3)
  # check: SigX3[(2*J+1):(3*J),(1*J+1):(2*J)] == SigX(3,2)   #这里其实是stacking matrix操作
  SigX <- SigX3; rm(SigX2, SigX3)
  gamX_est <- function(k){as.vector(solve(Ax[[k]])%*%Bx[,k])}  #Ax 和 Bx结合就是multi regression下的effect estimator这个summary statistic
  gamX_est <- sapply(1:K,gamX_est)  #dim(gamX_est) #J  K  #没问题



  # outcome quantities of interest
  ay <- 1/((ny*sy^2)+by^2) #为什么要加上一个by^2再做倒数？ -已证；没问题
  Ay <- (sqrt(ay)%*%t(sqrt(ay)))*ld
  By <- ay*by
  SigY <- solve(Ay)*(1-as.numeric(t(By)%*%solve(Ay)%*%By))*(1/(ny-J+1))#后半段感觉像是SSE/(n-p),但是为什么是1- ，并且df为什么是ny-J+1
  #注意SigY是stable的吗？错！是->0的！；所以对于传统的multivariate regression的estimator variance matrix, 需要乘以n用来stablize！
  gamY_est <- as.vector(solve(Ay)%*%By)#这个是estimated G-Y assocation

  #没问题: gamY_est和SigY就是multi regression下的GY effect estimator 和 variance matrix (注意：SigY是->0的)



  ## compute conditional F-statistics
  if(K>2){
    condF <- function(j){
      SigXX <- SigX[c(((((j-1)*J)+1):(j*J)),((1:(J*K))[-((((j-1)*J)+1):(j*J))])),c(((((j-1)*J)+1):(j*J)),((1:(J*K))[-((((j-1)*J)+1):(j*J))]))]
      g <- function(bet){as.vector(gamX_est[,j] - (gamX_est[,-j]%*%bet))}
      Q.gg <- function(bet){as.numeric(t(g(bet))%*%g(bet))}
      bet.gg <- nlminb(rep(0,(K-1)),objective=Q.gg)$par
      Om.nr <- function(bet){as.matrix(cbind(diag(J),kronecker(t(-bet),diag(J)))%*%SigXX%*%t(cbind(diag(J),kronecker(t(-bet),diag(J)))))}
      Q.nr <- function(bet){as.numeric(t(g(bet))%*%solve(Om.nr(bet))%*%g(bet))}
      G <- -gamX_est[,-j]
      DQ.nr <- function(bet){2*as.matrix(t(G)%*%solve(Om.nr(bet))%*%g(bet))}
      condF <- (nlminb(bet.gg,objective=Q.nr,gradient=DQ.nr)$objective/(J-K+1))
      return(condF)
    }
    condF <- sapply(1:K,condF)
  }

  if(K==2){
    condF <- function(j){
      SigXX <- SigX[c(((((j-1)*J)+1):(j*J)),((1:(J*K))[-((((j-1)*J)+1):(j*J))])),c(((((j-1)*J)+1):(j*J)),((1:(J*K))[-((((j-1)*J)+1):(j*J))]))]
      g <- function(bet){as.vector(gamX_est[,j] - (gamX_est[,-j]*bet))}
      Q.gg <- function(bet){as.numeric(t(g(bet))%*%g(bet))}
      bet.gg <- nlminb(rep(0,(K-1)),objective=Q.gg)$par
      Om.nr <- function(bet){as.matrix(cbind(diag(J),kronecker(t(-bet),diag(J)))%*%SigXX%*%t(cbind(diag(J),kronecker(t(-bet),diag(J)))))}
      Q.nr <- function(bet){as.numeric(t(g(bet))%*%solve(Om.nr(bet))%*%g(bet))}
      G <- -gamX_est[,-j]
      DQ.nr <- function(bet){as.numeric(2*as.matrix(t(G)%*%solve(Om.nr(bet))%*%g(bet)))}
      condF <- (nlminb(bet.gg,objective=Q.nr,gradient=DQ.nr)$objective/(J-K+1))
      return(condF)
    }
    condF <- sapply(1:K,condF)
  }


  ## continuously-updating GMM inference

  # preliminary estimate for initial value
  g <- function(bet){as.vector(gamY_est - (gamX_est%*%bet))}

  Q.gg <- function(bet){as.numeric(t(g(bet))%*%g(bet))}
  bet.gg <- nlminb(rep(0,K),objective=Q.gg)$par    #get a initial starting point

  # GMM estimate
  phi <- function(bet){kronecker(t(bet),diag(J))} #checked 只是为了stacking成一个matrix而已 和论文上的kronecker操作一模一样
  Om <- function(bet){as.matrix(SigY+(phi(bet)%*%SigX%*%t(phi(bet))))}  #checked! 注意是two-sample的情况； Om都是stable的 #SigX就是论文里的V
  Q <- function(bet){as.numeric(t(g(bet))%*%solve(Om(bet))%*%g(bet))} #checked #Q -> 0 as n -> infty

  G <- -gamX_est#checked

  DQ <- function(bet){2*as.matrix(t(G)%*%solve(Om(bet))%*%g(bet))} #辅助函数而已

  gmm <- nlminb(bet.gg,objective=Q,gradient=DQ)$par #optimization
  var.gmm <- as.matrix(   solve( t(G)%*%solve(Om(gmm))%*%G )   )  #all elements use the estimates
  Qstat <- Q(gmm)
  Q.pval <- pchisq(Qstat, df = J-K, lower.tail = F)  #Q statistics and Q testing #不用乘以n吗？ 专这样Q(gmm)不会随着样本增大然后趋近于0吗？
  #不用乘，因为Omega是带来1/n*stable_elements的，那么1/n通过^{-1}被提出去了；也就是n*Qstat的情况
  gmm.pval <- 2*pnorm(-abs(gmm/sqrt(diag(var.gmm))))


  LM<-NA; LM.pval<-NA
  if(   !is.na(as.matrix(beta0)[1,1])      ){
    ###LM----------------------------------
    ## identification-robust inference with LM statistics
    Om0 <- Om(beta0) #Omega matrix

    evec0 <- eigen(Om0)$vectors; eval0 <- eigen(Om0)$values
    sqrt.Om0 <- evec0%*%diag(sqrt(eval0))%*%t(evec0) #Omega^{1/2}

    g0 <- g(beta0)

    #Step1: 算Delta
    del <- SigX%*%t(phi(beta0))  #由于是two-sample 所以cov( \hat{G} , \hat{g}(\beta_0) )可以只看\hat{\gamma}部分，极大的简化了运算！
    #dim(del) JK * J #转置了一下 和之后的操作相匹配

    #Step2: 算D
    ind <- matrix(1:(J*K),ncol=J,byrow=TRUE)#只是方便定位matrix位置的复制位置index矩阵而已
    D <- list()
    for (k in 1:K){
      D[[k]] <- G[,k]-(del[ind[k,],]%*%solve(Om0)%*%g0) #checked
    }

    #Step3：算LM statistic
    D0 <- do.call(cbind, D); rm(D)

    P0 = solve(sqrt.Om0)%*%D0 #checked

    U0 = solve(sqrt.Om0)%*%g0 #checked

    P = P0%*%solve(t(P0)%*%P0)%*%t(P0) #checked

    LM = as.numeric(t(U0)%*%P%*%U0) #checked
    LM.pval <- 1-pchisq(LM,K)
    }

  res.list <- list("condFstat"=condF,
                 "gmm_est"=gmm,
                 "gmm_se"=sqrt(diag(var.gmm)),
                 'variance_matrix'=var.gmm,
                 "gmm_pval"=2*pnorm(-abs(gmm/sqrt(diag(var.gmm)))),
                 "Q_stat"=Qstat,
                 "Q_pval"=Q.pval,
                 "lm_stat"=LM,
                 "lm_pval"=LM.pval)
  return(res.list)
}


getLM_twosample<-function(bx,by,sx,sy,nx,ny,ld,cor.x,beta0){

  J = nrow(bx); K = ncol(bx)

  # exposures quantities of interest
  #(即用univariate regression的summary statistics转变成multi regression的summary statistics)
  ax <- matrix(NA,nrow=J,ncol=K); Ax <- list()
  for (k in 1:K){ax[,k] <- 1/((nx*sx[,k]^2)+bx[,k]^2)}
  for (k in 1:K){Ax[[k]] <- (sqrt(ax[,k])%*%t(sqrt(ax[,k])))*ld}
  Bx <- function(k){ax[,k]*bx[,k]}; Bx <- sapply(1:K,Bx)

  sqrt.Ax <- function(k){                                      #i.e. 对角化然后再求^{1/2}
    evec <- eigen(Ax[[k]])$vectors; eval <- eigen(Ax[[k]])$values
    return((evec%*%diag(sqrt(eval))%*%t(evec)))
  }
  sqrt.Ax <- lapply(1:K,sqrt.Ax)

  SigX <- function(k,l){  #这些SigX的操作主要是为了算论文中的V这个stacking matrix
    solve(sqrt.Ax[[k]]%*%t(sqrt.Ax[[l]]))*(   cor.x[k,l]-(  as.numeric(t(Bx[,k])%*%solve(sqrt.Ax[[k]]%*%t(sqrt.Ax[[l]]))%*%Bx[,l])  )   )*(1/(nx-J+1))
  }
  SigX2 <- function(m){
    SigX2a <- list()
    for (m1 in 1:K){SigX2a[[m1]] <- SigX(m1,m)}
    SigX2a <- do.call(rbind, SigX2a)
    return(SigX2a)
  }
  SigX3 <- list()
  for (m1 in 1:K){SigX3[[m1]] <- SigX2(m1)}
  SigX3 <- do.call(cbind, SigX3)
  # check: SigX3[(2*J+1):(3*J),(1*J+1):(2*J)] == SigX(3,2)   #这里其实是stacking matrix操作
  SigX <- SigX3; rm(SigX2, SigX3)
  gamX_est <- function(k){as.vector(solve(Ax[[k]])%*%Bx[,k])}  #Ax 和 Bx结合就是multi regression下的effect estimator这个summary statistic
  gamX_est <- sapply(1:K,gamX_est)  #dim(gamX_est) #J  K  #没问题


  # outcome quantities of interest
  ay <- 1/((ny*sy^2)+by^2) #为什么要加上一个by^2再做倒数？ -已证；没问题
  Ay <- (sqrt(ay)%*%t(sqrt(ay)))*ld
  By <- ay*by
  SigY <- solve(Ay)*(1-as.numeric(t(By)%*%solve(Ay)%*%By))*(1/(ny-J+1))#后半段感觉像是SSE/(n-p),但是为什么是1- ，并且df为什么是ny-J+1
  #注意SigY是stable的吗？错！是->0的！；所以对于传统的multivariate regression的estimator variance matrix, 需要乘以n用来stablize！
  gamY_est <- as.vector(solve(Ay)%*%By)#这个是estimated G-Y assocation


  g <- function(bet){as.vector(gamY_est - (gamX_est%*%bet))}
  phi <- function(bet){kronecker(t(bet),diag(J))} #checked 只是为了stacking成一个matrix而已 和论文上的kronecker操作一模一样
  Om <- function(bet){as.matrix(SigY+(phi(bet)%*%SigX%*%t(phi(bet))))}  #checked! 注意是two-sample的情况； Om都是stable的 #SigX就是论文里的V

  G <- -gamX_est#checked

  ###LM-------------------------------------------------------------------------
  ## identification-robust inference with LM statistics
  Om0 <- Om(beta0) #Omega matrix

  evec0 <- eigen(Om0)$vectors; eval0 <- eigen(Om0)$values
  sqrt.Om0 <- evec0%*%diag(sqrt(eval0))%*%t(evec0) #Omega^{1/2}

  g0 <- g(beta0)

  #Step1: 算Delta
  del <- SigX%*%t(phi(beta0))  #由于是two-sample 所以cov( \hat{G} , \hat{g}(\beta_0) )可以只看\hat{\gamma}部分，极大的简化了运算！
  #dim(del) JK * J #转置了一下 和之后的操作相匹配

  #Step2: 算D
  ind <- matrix(1:(J*K),ncol=J,byrow=TRUE)#只是方便定位matrix位置的复制位置index矩阵而已
  D <- list()
  for (k in 1:K){
    D[[k]] <- G[,k]-(del[ind[k,],]%*%solve(Om0)%*%g0) #checked
  }

  #Step3：算LM statistic
  D0 <- do.call(cbind, D); rm(D)

  P0 = solve(sqrt.Om0)%*%D0 #checked

  U0 = solve(sqrt.Om0)%*%g0 #checked

  P = P0%*%solve(t(P0)%*%P0)%*%t(P0) #checked

  LM = as.numeric(t(U0)%*%P%*%U0) #checked
  LM.pval <- 1-pchisq(LM,K)

  res.list <- list("lm_stat"=LM,
                   "lm_pval"=LM.pval)
  return(res.list)

}


# # ###examples
# set.seed(1123)
# N<-10000 #sample size  #用N<-50000 ?
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
# #prepare: Z X Y J K
#
# #bx,by,sx,sy,nx,ny,ld,cor.x,beta0=NA
#
# ###summary information
#
# ###one sample (current sample) for (G X); another for (G Y)
#
#
# #bx sx nx
# bx<-matrix( NA, J,K ); sx<-matrix( NA, J,K )
# for(j in 1:J){
#   for(k in 1:K){
#     GXfit<-lm( X[,k] ~ Z[,j] )
#     bx[j,k]<-as.numeric(  summary( GXfit )$coef[2,1]  )
#     sx[j,k]<-as.numeric(  summary( GXfit )$coef[2,2]  )
#   }
# }
# nx<-nrow(X)
#
#
# #ld: J x J genetic variant correlation matrix
#
# ld<-cor(  Z    ) #checked; dim(ld) #30 30
#
# #cor.x: K x K matrix exposure correlation matrix
# cor.x<-cor( X )  #dim(cor.x  ) #3 3
#
#
# #another sample: (Z2 Y2)
# set.seed(1123+1)
#
# Z2<-matrix(  rbinom(N*J,2,0.3), N,J)  #three-valued SNPs
# U2<-matrix(  rnorm(N*K,0,1), N,K )
# X2<- Z%*%Alphas + U #exposure matrix
# Y2<- X%*%Beta + apply( U, 1, sum )
#
# #by sy ny
# by<-rep(NA,J); sy<-rep(NA,J)
# for(j in 1:J){
#   GYfit<-lm(  Y2 ~  Z2[,j] )
#   by[j]<-as.numeric(  summary( GYfit )$coef[2,1]  )
#   sy[j]<-as.numeric(  summary( GYfit )$coef[2,2]  )
# }
# ny<-length(Y2)
#
# #GMM
# gmm_res<-gmm_lm_twosample( bx,by,sx,sy,nx,ny,ld,cor.x,beta0=NA )
# Beta
# gmm_res$gmm_est
# gmm_res$gmm_se
# gmm_res$Q_pval #为什么很小于0.05? 试试用自己的Qfunction?
#
#
# LMres<-getLM_twosample(  bx,by,sx,sy,nx,ny,ld,cor.x,beta0=Beta )
# LMres$lm_pval
#
# LMres<-getLM_twosample(  bx,by,sx,sy,nx,ny,ld,cor.x,beta0=c(10,20,30) )
# LMres$lm_pval
