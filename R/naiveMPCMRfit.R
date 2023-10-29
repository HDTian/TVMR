
#naive association fit其实只能是One-sample;且Gmatrix不重要
naiveMPCMRfit<-function(Gmatrix,
                        Yvector,
                        res, #res<-FPCA(...)
                        nPC=NA, #how many nPC to be retained? default is the the just-exceed-to-95%
                        XYmodel=NA #if labelled as '0' '1' '2', ...; use this as the true effect (real data cannot have this)
){
  fitRES<-list()
  if(is.na(nPC)){   nPC<-sum(res$cumFVE<0.95)+1    }
  if( nPC>ncol( res$xiEst ) ){stop('The exposure curves cannot support much information; please use smaller nPC')  }
  fitRES$nPC_used<-nPC #nPC_used就是nPC=K
  PCmatrix<-res$xiEst[,1:nPC]
  PC_<-PCmatrix[!is.na(  PCmatrix[,1] ),] #removing missing values
  if( nrow(PCmatrix)!=nrow(PC_)){ 
    cat('Not that there exists NA values for some individuals in fPCA results; check the measured timepoints for all individuals', '\n')}
  if( sum(is.na(Gmatrix))>0   ){stop('there exists missing data for the Gmatrix') }
  if(   nrow( Gmatrix)!=nrow(PC_  )    ){
    stop('the nrow of the Gmatrix is not equal to the nrow of the principal components matrix (after removing possible missing-PC individuals)')
  }
  J<-ncol(Gmatrix)
  K<- ncol(PCmatrix)#nPC<-K
  
  PC<-res$xiEst[,1:nPC] #PC matrix  #PC=PC_=PCmatrix
  #cor(PC[,1],RES$details$UUU[,ncol(RES$details$UUU)])
  #fit<-lm( DAT$Y ~ PC ) #naive regression  #注意：SLR = OLS = MLE = GMM 
  fit<-lm( Yvector ~ PC ) #naive regression  #注意：SLR = OLS = MLE = GMM 
  
  MPCMRest<-coef(fit)[-1] #此时的MPCMR就是MPCregression；并不是MR；也不存在random-effect 这一说
  fitRES$Est<-MPCMRest
  MPCMRvar<-vcov(fit)[-1,-1]#vcov(): Variance-Covariance Matrix
  fitRES$VM<-MPCMRvar
  pointwise_shape_var<-diag(   ( (res$phi)[,1:nPC]   )%*%MPCMRvar%*%t( (res$phi)[,1:nPC])  )#取diagonal;易理解
  ggdata<-data.frame(time=res$workGrid, effect=(res$phi)[,1:nPC]%*%MPCMRest,
                     effect_low=(res$phi)[,1:nPC]%*%MPCMRest-1.96* sqrt(pointwise_shape_var  ),
                     effect_up=(res$phi)[,1:nPC]%*%MPCMRest+1.96* sqrt(pointwise_shape_var  ),
                     LM_low=NA,
                     LM_up=NA,
                     true_shape=NaN)
  
  if(!is.na(XYmodel)){
    if( !XYmodel%in%c('0','1','2','3','4','5','3.5') ){ stop('please use correct X-Y model') }
    
    if(XYmodel=='0'){
      fun<-function(t){  0*(t<Inf)   }   #Lifetime null effect
    }
    if(XYmodel=='1'){
      fun<-function(t){  0.1*(t<Inf)   }  #Lifetime constant effect
    }
    if(XYmodel=='2'){
      fun<-function(t){  0.02*t   }  #Lifetime increasing effect
    }
    if(XYmodel=='3'){
      fun<-function(t){  0.5-0.02*t   }  #Lifetime sign-change effect
    }
    if(XYmodel=='4'){
      fun<-function(t){  0.1*(t<20)   } #Early-age-only effect
    }
    if(XYmodel=='5'){
      fun<-function(t){  0.1*(t>30)   }  #Later-age-only effect
    }
    if(XYmodel=='3.5'){
      fun<-function(t){  0.4-0.02*t   }  #Lifetime sign-change effect
    }
    if(XYmodel=='6'){
      fun<-function(t){  0.05*(-t+20)*(t<20)   } #Early-age-only effect (continuous threshold effect)
    }
    if(XYmodel=='7'){
      fun<-function(t){  0.05*(t-30)*(t>30)   } #Later-age-only effect (continuous threshold effect)
    }
    
    ggdata$true_shape<-fun( res$workGrid   )
    
  }
  fitRES$ggdata<-ggdata
  p<- ggplot(ggdata, aes(time, effect))+
    geom_hline(yintercept = 0,linewidth=0.5,linetype = 2,col='grey' )+
    geom_line(ggdata, mapping =aes(time, true_shape), alpha=1,linewidth=1,col='blue'  )+
    geom_line(ggdata, mapping =aes(time, effect), alpha=1,linewidth=1  )+
    geom_line(ggdata, mapping =aes(time, effect_low), alpha=1,linewidth=1,linetype = 2  )+
    geom_line(ggdata, mapping =aes(time, effect_up), alpha=1,linewidth=1,linetype = 2  )+
    labs(x='Age',y='Time-varying effect')+
    theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    coord_cartesian( ylim = c(min(ggdata$effect)-0.5  ,   max(ggdata$effect)+0.5 ) )
  p
  
  fitRES$p<-p
  
  
  ###MSE and coverage rate (最好是返回一个vector over res$workGrid; 方便后续局部分析)
  SE<- ( (res$phi)[,1:nPC]%*%MPCMRest - fun(  res$workGrid ) )^2 #Squared errors #mean(SE)->MSE
  
  ###coverage analysis
  if(!is.na(XYmodel)){
    Co<-abs( (res$phi)[,1:nPC]%*%MPCMRest  - fun( res$workGrid)   )< 1.96* sqrt(pointwise_shape_var  )
    #covaer-or-not vector #mean(Co)->coverage rate
  }else{
    Co<-NA  #no true shape, no coverage
  }
  fitRES$SE<-SE;fitRES$MSE<-mean(SE)
  fitRES$Co<-Co;fitRES$Coverage_rate<-mean(Co)
  
  ###significant timepoints information
  fitRES$time_points<-res$workGrid
  fitRES$sig_points<-as.numeric(   (ggdata$effect_up*ggdata$effect_low)>0   )*(2*as.numeric( ggdata$effect_low>0  )-1 )
  #前者返回 0或1 (显著与否)； 后者返回-1或1
  
  
  return(fitRES)
}
