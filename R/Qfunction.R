


#Parametric testing function (with B must be inputted)
#if let B is a square matrix (i.e. K=L), then it is just a nonparametric testing (equivalent to test the IV validity)

#Qfunction also give robust MR estimator and CIs (based on weight-changed IVW regression)


#Qfunction: according to the inserted effect vector to return the Q statistics and (BTW) the robust MR estimate and s.e./CI at the current v effect vector 
#this inserted style is used for iteration style for geting the final Q statistics (and BTW the robust MR inference results) - like Jack Bowden's style 

#Gam: the covariance between the G-X estimate and G-Y estimate for each SNPs
Qfunction_old<-function(v,#v: vector #the initial values used in the weight
                    by, #G-Y est
                    byse, #s.e.(G-Y est)
                    B,#B matrix; i.e. \int basisfunction(t)*eigenfunction(t) dt matrix   dim(B)=K*L
                    BX, #G-X est
                    BXse, #s.e.(G-X est)
                    Gam #Gamma matrix  #dim(Gam): J K  #无所谓是否为over-lapping data；反正都会体现在Gam中
){
  ###let possible vector into matrix
  BX<-as.matrix(BX);BXse<-as.matrix(BXse)

  RES<-list()
  RES$original_input_v<-v
  ###first get the estimated vector
  ##dummy weighted regression #because linear regression fitting = OLS = MLE = GMM
  J<-nrow(BX )
  weights<-c()  #用在依次每个SNPs的weight中; weight就是Q formula里的每个SNP的分母部分呀
  for(j in 1:J){
    weight<-byse[j]^2+t(  B%*%v )%*%diag(  BXse[j,]^2, nrow=length( BXse[j,]) )%*%( B%*%v ) - 2*t(v)%*%t(B)%*%Gam[j,]
    weights<-c(weights,weight)
  }

  S<-diag(    weights ) #dummy Sigma matrix
  BBX<-BX%*%B #BBX: BX+B #Alpha%*%B    #dim(BBX) #30 2

  Est<-solve(    t(BBX)%*%solve(S)%*%BBX )%*%t(BBX)%*%solve(S)%*%by  #this is interpreated as the robust MR ests

  #################################Robust statistical inference#######################################
  ##random-effect setting
  tau2<-max(1,   t(by- BBX%*%Est)%*%solve( S )%*%(by- BBX%*%Est)/(J-length(Est))      )

  #variance matrix of the IVW regression for MPCMR
  Estvar<-solve(    t(BBX)%*%solve( S )%*%BBX            )*tau2

  #estimate and s.e.s
  inference_results<-cbind( as.numeric(Est) , sqrt( diag( Estvar  ) )    )
  colnames(inference_results)<-c('estimate','s.e.')
  RES$inference_results<-inference_results
  RES$Estvar<-Estvar  #variance matrix
  ####################################################################################################

  ###then get the updated Q values
  Qvalue<-sum( (by-BBX%*%Est)^2/weights  )


  RES$Est<-as.vector(Est)  #kind of robust estimates
  RES$Qvalue<-Qvalue
  return(RES)
}


###for simple IV pleiotropy effect/weak instrument testing for PCs:
#B==identify matrix


#已证：无论是用之前的B分离出G-X effect variance matrix的形式还是直接用transformed exposure的形式；Q statistics的result是一模一样的



#Gam: the covariance between the G-X estimate and G-Y estimate for each SNPs
Qfunction<-function(v,#v: vector #the initial values used in the weight #length(v)=L
                    by, #G-Y est #length(by)=J
                    byse, #s.e.(G-Y est) #length(byse)=J
                    B,#B matrix; i.e. \int basisfunction(t)*eigenfunction(t) dt matrix   #dim(B)=K*L
                    BX, #G-X est #dim(BX) #J K
                    Sigma,#a array (K*K*J) contains the covariance matrix for G-PCs assocation estimators
                    Gam #Gamma matrix  #dim(Gam): J K  #无所谓是否为over-lapping data；反正都会体现在Gam中
){
  ###let possible vector into matrix
  BX<-as.matrix(BX)
  
  RES<-list()
  RES$original_input_v<-v
  ###first get the estimated vector
  ##dummy weighted regression #because linear regression fitting = OLS = MLE = GMM
  J<-nrow(BX )
  weights<-c()  #用在依次每个SNPs的weight中; weight就是Q formula里的每个SNP的分母部分呀
  for(j in 1:J){
    Sigma_j<-Sigma[,,j]
    weight<-byse[j]^2+t(  B%*%v )%*%Sigma_j%*%( B%*%v ) - 2*t(v)%*%t(B)%*%Gam[j,]
    weights<-c(weights,weight)
  }
  
  S<-diag(    weights ) #dummy Sigma matrix
  BBX<-BX%*%B #BBX: BX+B #Alpha%*%B    #dim(BBX) #30 2
  
  Est<-solve(    t(BBX)%*%solve(S)%*%BBX )%*%t(BBX)%*%solve(S)%*%by  #this is interpreted as the robust MR ests
  
  #################################Robust statistical inference#######################################
  ##random-effect setting
  tau2<-max(1,   t(by- BBX%*%Est)%*%solve( S )%*%(by- BBX%*%Est)/(J-length(Est))      )
  
  #variance matrix of the IVW regression for MPCMR
  Estvar<-solve(    t(BBX)%*%solve( S )%*%BBX            )*tau2
  
  #estimate and s.e.s
  inference_results<-cbind( as.numeric(Est) , sqrt( diag( Estvar  ) )    )
  colnames(inference_results)<-c('estimate','s.e.')
  RES$inference_results<-inference_results
  RES$Estvar<-Estvar  #variance matrix
  ####################################################################################################
  
  ###then get the updated Q values
  Qvalue<-sum( (by-BBX%*%Est)^2/weights  )
  
  
  RES$Est<-as.vector(Est)  #kind of robust estimates
  RES$Qvalue<-Qvalue
  return(RES)
}
