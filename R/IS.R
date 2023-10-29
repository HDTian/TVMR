###weak instrument strength for MVMR with conditional F statistics (individual data + one-sample - i.e. conditional F)



###one-sample -> conditional F statistic

IS<-function(J,# # of SNPs
             K,# number of measured timpoints
             timepoints,  #timepoints vector
             datafull){  # cbind( SNPs , Xs  )
  #Instrument strength check (conditional F statistic value) with individual data
  X<-data.matrix(datafull)  #dim(datafull)[1] = 30 + 51 + 1  (30=J)  timepoints<-1:51
  FF<-rep(NA,K)  #F statistic value
  RR<-rep(NA,K)  #R^2     i.e. proportion of variance in the exposure explained by the genetic variants
  cFF<-rep(NA,K) #conditional F statistic value
  for(k in 1:K){ #generate the conventional weak instrument F statistic value
    t<-timepoints[k]  #当前的这个k-th variable就是作为'outcome'的variable
    fitX<-lm( X[,J+t ]~matrix(   as.numeric(X[,1:J ]) ,  dim(X)[1], J    )   )
    FF[k]<-as.numeric(summary(fitX)$fstatistic[1])
    RR[k]<-as.numeric(summary(fitX)$r.squared[1])
  }
  for(k in 1:K){  #conditional F 其实就是2SLS版本的F statistics
    XXm<-c()
    for(kk in (1:K)[-k]){#当前的这个k-th variable就是作为'outcome'的variable
      t<-timepoints[kk]  
      fitX<-lm( X[,J+t ]~matrix(   as.numeric(X[,1:J ]) ,  dim(X)[1], J    )   )
      Xhat<- X[,J+t ] - summary(fitX)$residuals  #checked
      XXm<-cbind(XXm,as.numeric(Xhat))
    }
    t<-timepoints[k]
    fitX<-lm( X[,J+t ]~XXm     )  #其实就是2SLS呀！!!!!!! #当前的这个k-th variable就是作为'outcome'的variable
    resi<-as.numeric( summary(fitX)$residuals )
    fit<-lm(  resi  ~ matrix(   as.numeric(X[,1:J ]) ,  dim(X)[1], J    )  )
    cFF[k]<-as.numeric(summary(fit)$fstatistic[1])*summary(fit)$fstatistic[2]/(summary(fit)$fstatistic[2]-(K-1)    )     #conditional F statistic value  #adjusted!!!
  }
  return( cbind(   timepoints, RR ,FF , cFF    ) )
}

# #examples
# PC_<-PC[,1:nPC]
# J ; K<- nPC  ;  timepoints<- 1:nPC  ; datafull<- cbind(G,PC_)  #datafull 得是Gx+Ex的matrix
# IS(J,K,timepoints, datafull )#threshold
