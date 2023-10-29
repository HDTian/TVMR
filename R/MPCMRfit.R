#MPCMR: MVMR + functional back estimation ;in the middle: weak instrument analysis (conditional F + Q)


#MPCMRfit: 1. nonparametric fitting +2.  polynomial parametric fitting + 3. two-step test
#res: usually FPCA result - a list must contains: $xiEst $cumFVE  $phi  $workGrid

#ParMPCMRfit: (这个更好用；且更general)： get the 1. weak IV results 2. curve estimates, and 3. the Q test based on current L (L=nPC is just the nonparametric testing result)
#已证：用L=nPC的ParMPCMRfit得到的Q statistics result和直接用MVMR(无视任何basisfunction和eigenfunctions)的Q statistics结果是一样的；所以直接用作IV valisity test

#note: for comparasion; the weak IV analysis in MPCMRfit and ParMPCMRfit will only use naive Q style


MPCMRfit<-function(   Gmatrix, #G matrix for PCs
                      Gymatrix=NA,  #G matrix for Y
                      Yvector, #Y 未来其实可以带时间戳，方便做time-varying outcome analysis
                      IDmatch=NA, #以Gmatrix为先手主导的match vector； 默认为     1: nrow(Gmatrix)
                      res, #res<-FPCA(...);如果都导入FPCA res了，那就没必要搞PCmatrix了？错！PCmatrix可以人工选择！
                      nPC=NA, #how many nPC to be retained? default is the the just-exceed-to-95%
                      XYmodel=NA #if labelled as '1' '2', ...; use this as the true effect (real data cannot have this)
){
  fitRES<-list()
  if(is.na(nPC)){   nPC<-sum(res$cumFVE<0.95)+1    }
  if( nPC>ncol( res$xiEst ) ){stop('The exposure curves cannot support much information; please use smaller nPC')  }
  fitRES$nPC_used<-nPC #nPC_used就是nPC=K

  PCmatrix<-res$xiEst[,1:nPC];PCmatrix<-as.matrix(PCmatrix   ) #以防nrol=1的matrix退化成vector
  PC_<-PCmatrix[!is.na(  PCmatrix[,1] ),] #removing possible missing values
  PC_<-as.matrix(PC_   ) #以防nrol=1的matrix退化成vector

  if( nrow(PCmatrix)!=nrow(PC_)){
    cat('Not that there exists NA values for some individuals in fPCA results; check the measured timepoints for all individuals', '\n')}
  if( sum(is.na(Gmatrix))>0   ){stop('there exists missing data for the Gmatrix') }
  if(   nrow( Gmatrix)!=nrow(PC_  )    ){
    stop('the nrow of the Gmatrix is not equal to the nrow of the principal components matrix (after removing possible missing-PC individuals)')
  }


  is.md<-function(x){ return(   is.matrix(x)|is.data.frame(x)     )  } #is matrix  or data.frame?


  if((length(IDmatch)==1)&(is.md(Gymatrix))){
    stop('Carefully!: you decide to use different Gmatrix and Gymatrix but not define the IDmatch. Define the IDmatch first.')
  }
  if(length(IDmatch)==1){IDmatch<-1:nrow(Gmatrix)}
  if(length(IDmatch)!=nrow(Gmatrix)){
    stop('The length of IDmatch is not equal to the nrow of Gmatrix')
  }



  if(!is.md(Gymatrix)){ Gymatrix<-Gmatrix }  #complete one sample case
  if( ncol(Gmatrix)!= ncol( Gymatrix )    ){
    stop('Gmatrix and Gymatrix do not have the same columns' )
  }
  if( nrow( Gymatrix)!=length(Yvector  )  ){
    stop('the nrow of the Gymatrix is not equal to the length of the Y vector')
  }

  J<-ncol(Gmatrix)
  K<- ncol(PCmatrix)#nPC<-K

  ###Weak instrument strength analysis
  if(K==1){#only one PC case - UVMR
    fitGX<-lm(PC_~  as.matrix( G )  )
    ISres<-c(   as.numeric(summary(fitGX)$r.squared[1]) ,   as.numeric(summary(fitGX)$fstatistic[1])      )
    names( ISres  )<-c('RR','F')
    fitRES$ISres<-ISres
  }else{#IS只有K>=2才能用
    timepoints<- 1:K
    datafull<- cbind(Gmatrix,PC_)  #datafull 得是Gx+Ex的matrix for IS
    ISres<-IS(J,K,timepoints, datafull ); rownames(ISres)<-paste0('PC',1:K)
    fitRES$ISres<-ISres[,-1]
  }






  ###nonparametric MVMR fitting (ignore the G-PC estimate uncertainty)------------------
  ###-----------------------------------------------------------------------------------
  G<-Gmatrix ; Y<-Yvector  ; Gy<-Gymatrix

  #G-X fitting
  BX<-c(); BXse<-c()
  nPC<-K
  for(p in 1:nPC){  #只用前nPC个
    fitGX<-lm(PC_[,p]~  as.matrix( G )  )
    bx<-as.numeric(summary(fitGX)$coef[-1,1]) ;  bxse<-as.numeric(summary(fitGX)$coef[-1,2])
    BX<-cbind(BX, bx  );   BXse<-cbind( BXse, bxse )
  }


  #IV-PC1 IV-PV2 scatterplot
  if( nPC>=2  ){
    ggdata<-data.frame(   beta1=BX[,1], beta2=BX[,2],
                          beta1_error_up=BX[,1]+1.96*BXse[,1], beta1_error_low=BX[,1]-1.96*BXse[,1] ,
                          beta2_error_up=BX[,2]+1.96*BXse[,2], beta2_error_low=BX[,2]-1.96*BXse[,2])
    scatterp<-     ggplot(ggdata, aes(beta1, beta2))+
      geom_errorbar(data=ggdata,mapping=aes(x=beta1,ymin=beta2_error_low,ymax=beta2_error_up),width = 0,color='#666666')+
      geom_errorbar(data=ggdata,mapping=aes(y=beta2,xmin=beta1_error_low,xmax=beta1_error_up),width = 0,color='#666666')+
      geom_point(data=ggdata,mapping=aes(x=beta1,y=beta2),size=1.8)+
      geom_hline(yintercept = 0,linetype='dashed',col='#999999')+
      geom_vline(xintercept = 0,linetype='dashed',col='#999999')+
      geom_abline(intercept=0,slope=1,linetype='dashed',col='#999999')+
      labs(x = 'G-X ests for PC1', y = 'G-X ests for PC2')+
      theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
    fitRES$scatterp<-scatterp
  }


  ###IS: instrument Strength with Q style value----------------------------
  if(nPC>=2){
    Gam<-matrix(0,nrow=J,ncol=nPC-1) #Gam直接固定成a vector of 0，因为完美fPCA下，不同的phenotype (PCs) 之间是covariance-zero的
    IS_Q<-c()
    for(p in 1:nPC){ #每个PCs都跑一遍Q test #当前的这个k-th variable就是作为'outcome'的variable
      B_<-diag(1, nPC-1  )
      by_<-BX[,p]; byse_<-BXse[,p]; BX_<-as.matrix(BX[,-p]); BXse_<-as.matrix(BXse[,-p])
      #initial starting estimate
      Qfunction_res<-Qfunction_old(  v=rep(0,nPC-1  ),by=by_,byse=byse_, B=B_,BX=BX_,BXse=BXse_,Gam=Gam)
      v0<-Qfunction_res$Est
      #iteratiing to get stable robust estimate
      iter_time<-0;total_dif<-Inf
      while( abs(total_dif)>1e-10  ){
        iter_time<-iter_time+1
        v_updated<-Qfunction_old( v=v0,by=by_, byse=byse_,B=B_,BX=BX_,BXse=BXse_,Gam=Gam)$Est
        total_dif<-sum( abs(v_updated-v0) )
        v0<-v_updated
      }
      #get the <Q values>, <df> and <p-value> for the p-th PC, and storage into fitRES$ISres
      Qvalue<-Qfunction_old( v=v0,by=by_, byse=byse_,B=B_,BX=BX_,BXse=BXse_,Gam=Gam)$Qvalue  #qchisq(0.95,J-(nPC-1))
      df<-J-(nPC-1)
      pvalue<-1-pchisq(Qvalue,J-(nPC-1))
      IS_Q<-rbind(IS_Q,  c(Qvalue ,df,  pvalue ) )
    }
    #storage into fitRES$ISres
    fitRES$ISres<-cbind( fitRES$ISres, IS_Q   )
    colnames(fitRES$ISres)<-c('RR','F','cF','Qvalue','df','pvalue') #pcalue越小越好，越应该拒绝！说明instrument strength没问题！
  }


  ###-----------------------------------------------------------------------

  #G-Y fitting
  fitGY<-lm( Y ~  as.matrix( Gy )  )
  by<-as.numeric(summary(fitGY)$coef[-1,1]) ;  byse<-as.numeric(summary(fitGY)$coef[-1,2])

  S<-vcov(fitGY  )[-1,-1]  #variance-covariance matrix



  ###MVMR fitting------------------------------
  #point estimates
  MPCMRest<-solve(    t(BX)%*%solve( S )%*%BX            )%*%t(BX)%*%solve(S)%*%by
  fitRES$MPCMRest<-MPCMRest
  ##random-effect setting
  tau2<-max(1,   t(by- BX%*%MPCMRest)%*%solve( S )%*%(by- BX%*%MPCMRest)/(J-length(MPCMRest))      )
  fitRES$tau<-sqrt(  tau2 )
  #variance matrix of the IVW regression for MPCMR
  MPCMRvar<-solve(    t(BX)%*%solve( S )%*%BX            )*tau2
  fitRES$MPCMRvar<-MPCMRvar

  ##nonparametric fitting visualization results
  pointwise_shape_var<-diag(   ( (res$phi)[,1:nPC]   )%*%MPCMRvar%*%t( (res$phi)[,1:nPC])  )#取diagonal;易理解
  ggdata<-data.frame(time=res$workGrid, effect=(res$phi)[,1:nPC]%*%MPCMRest,
                     effect_low=(res$phi)[,1:nPC]%*%MPCMRest-1.96* sqrt(pointwise_shape_var  ),
                     effect_up=(res$phi)[,1:nPC]%*%MPCMRest+1.96* sqrt(pointwise_shape_var  ),
                     true_shape=NaN)
  if(!is.na(XYmodel)){
    #if( !XYmodel%in%c('0','1','2','3','4','5','3.5') ){ stop('please use correct X-Y model') }

    if(XYmodel=='0'){
      fun<-function(t){  0*(t<Inf)   }
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
  plotdif<- max(ggdata$effect_up)-min(ggdata$effect_low)
  p1<- ggplot(ggdata, aes(time, effect))+
    geom_hline(yintercept = 0,linewidth=0.5,linetype = 2,col='grey' )+
    geom_line(ggdata, mapping =aes(time, true_shape), alpha=1,linewidth=1,col='blue'  )+
    geom_line(ggdata, mapping =aes(time, effect), alpha=1,linewidth=1  )+
    geom_line(ggdata, mapping =aes(time, effect_low), alpha=1,linewidth=1,linetype = 2  )+
    geom_line(ggdata, mapping =aes(time, effect_up), alpha=1,linewidth=1,linetype = 2  )+
    labs(x='Age',y='Time-varying effect')+
    theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    coord_cartesian( ylim = c(min(ggdata$effect_low)-plotdif  ,   max(ggdata$effect_up)+plotdif ) )

  fitRES$p1<-p1 #Removed 51 rows containing missing values (`geom_line()`). #这些都是NaN


  ###MSE and coverage rate (最好是返回一个vector over res$workGrid; 方便后续局部分析)
  if(!is.na(XYmodel)){
    SE<- ( (res$phi)[,1:nPC]%*%MPCMRest - fun(  res$workGrid ) )^2 #Squared errors #mean(SE)->MSE
  }else{
    SE<-NA
  }
  fitRES$SE<-SE;fitRES$MSE<-mean(SE)

  ###coverage analysis
  if(!is.na(XYmodel)){
    Co<-abs( (res$phi)[,1:nPC]%*%MPCMRest  - fun( res$workGrid)   )< 1.96* sqrt(pointwise_shape_var  )
    #covaer-or-not vector #mean(Co)->coverage rate
  }else{
    Co<-NA  #no true shape, no coverage
  }
  fitRES$Co<-Co;fitRES$Coverage_rate<-mean(Co)

  ###significant timepoints information
  fitRES$time_points<-res$workGrid
  fitRES$sig_points<-as.numeric(   (ggdata$effect_up*ggdata$effect_low)>0   )*(2*as.numeric( ggdata$effect_low>0  )-1 )
  #前者返回 0或1 (显著与否)； 后者返回-1或1




  ###Two-step testing--------------------------------------------
  ###------------------------------------------------------------

  maxL<-K #K=nPC #complete model

  ##now decide the best L by backward stepwise LRT
  Qvalues<-c()
  for(l in maxL:1){  #l=maxL时，test为instrument validity test
    ParMPCMRFiting_res<-ParMPCMRfit( Gmatrix=Gmatrix,
                                     Gymatrix=Gymatrix,
                                     Yvector=Yvector,
                                     IDmatch=IDmatch,
                                     res=res,
                                     L=l,
                                     XYmodel=XYmodel,
                                     Fit=FALSE )
    Qvalues<-c(Qvalues,ParMPCMRFiting_res$Qvalue)
    if(l!=maxL){
      Qdif<-Qvalues[length(Qvalues)]-Qvalues[length(Qvalues)-1 ]
    }else{ Qdif<- -Inf    }
    ll<-l
    if(Qdif>qchisq(0.95,1  )){ break  }  #一旦break 就是ll=l；如果一路跑完到l=1都不break，则ll=l-1=0
    ll<-l-1
  }
  fitRES$Qvalues<-Qvalues
  L_used<- min(maxL, ll+1)
  fitRES$L_used<-L_used


  ##first-step testing results (first_step_test就是IV validity test)
  fitRES$first_step_test<-c( Qvalue=Qvalues[1], Qdf=J-maxL, pvalue=1-pchisq(Qvalues[1],J-maxL)    )

  ##second-step testing results
  ParMPCMRFiting_res<-ParMPCMRfit(Gmatrix=Gmatrix,
                                  Gymatrix=Gymatrix,
                                  Yvector=Yvector,
                                  IDmatch=IDmatch,
                                  res=res,
                                  L=L_used,
                                  XYmodel=XYmodel,
                                  Fit=TRUE )
  fitRES$second_step_test<-list(
    Qvalue=ParMPCMRFiting_res$Qvalue-Qvalues[1],
    Qdf=maxL-L_used,
    pvalue=1-pchisq(ParMPCMRFiting_res$Qvalue-Qvalues[1],maxL-L_used)    )

  ##testing results


  ###parametric fitting results with L_used
  fitRES$par_fit_res<-ParMPCMRFiting_res

  return(fitRES)

}













#Parametric MPCMR fitting results and Q test statistic values
ParMPCMRfit<-function(Gmatrix,  #G matrix for PCs
                      Gymatrix=NA,  #G matrix for Y
                      Yvector,
                      IDmatch=NA, #以Gmatrix为先手主导的match vector； 默认为1: nrow(Gmatrix)
                      res, #res<-FPCA(...)
                      nPC=NA, #how many nPC to be retained? default is the the just-exceed-to-95%
                      L=NA, #L<=nPC; otherwise ill-identification #default setting is NA, means L=nPC, corresponding to the nonparametric case
                      XYmodel=NA, #if labelled as '0' '1' '2', ...; use this as the true effect (real data cannot have this)
                      Fit=TRUE  #whether fit? or just get the test value
){
  fitRES<-list()


  PC<-res$xiEst[!is.na(  res$xiEst[,1] ),] #removing missing values
  PC<-as.matrix(PC   ) #以防nrol=1的matrix退化成vector
  if( nrow( res$xiEst)!=nrow(PC)){
    cat('Not that there exists NA values for some individuals in fPCA results; check the measured timepoints for all individuals', '\n')}
  if( sum(is.na(Gmatrix))>0   ){stop('there exists missing data for the Gmatrix') }
  if(   nrow( Gmatrix)!=nrow(PC  )    ){
    stop('the nrow of the Gmatrix is not equal to the nrow of the principal components matrix (after removing possible missing-PC individuals)')
  }

  is.md<-function(x){ return(   is.matrix(x)|is.data.frame(x)     )  } #is matrix  or data.frame?

  if((length(IDmatch)==1)&(is.md(Gymatrix))){
    stop('Carefully!: you used different Gmatrix and Gymatrix but not define the IDmatch. Define the IDmatch first.')
  }
  if(length(IDmatch)==1){IDmatch<-1:nrow(Gmatrix)}
  if(length(IDmatch)!=nrow(Gmatrix)){
    stop('The length of IDmatch is not equal to the nrow of Gmatrix')
  }


  if(!is.md(Gymatrix)){ Gymatrix<-Gmatrix }  #complete one sample case
  if( ncol(Gmatrix)!= ncol( Gymatrix )    ){
    stop('Gmatrix and Gymatrix do not have the same columns' )
  }
  if( nrow( Gymatrix)!=length(Yvector  )  ){
    stop('the nrow of the Gymatrix is not equal to the length of the Y vector')
  }




  J<-ncol(Gmatrix)
  N<-nrow(Gmatrix)
  if(is.na(nPC)){   nPC<-sum(res$cumFVE<0.95)+1    }

  #nPC<-sum(res$cumFVE<0.95)+1  #默认采取这个nPC,没有人为调整的余地 #更新：增加人为调整的余地
  fitRES$nPC<-nPC
  K<-nPC
  if(is.na(L)){  L<- K  } #the default setting
  if(L>K){stop(paste0('L>',nPC,': ill-identification, please use smaller L' ) )}


  ###BX BXse by byse
  G<-Gmatrix ; Y<-Yvector; Gy<-Gymatrix
  #G-X fitting
  BX<-c(); BXse<-c()
  nPC<-K #K<-nPC 其实绕了一圈而已
  for(p in 1:nPC){  #只用前nPC个
    fitGX<-lm(PC[,p]~  as.matrix( G )  )
    bx<-as.numeric(summary(fitGX)$coef[-1,1]) ;  bxse<-as.numeric(summary(fitGX)$coef[-1,2])
    BX<-cbind(BX, bx  );   BXse<-cbind( BXse, bxse )
  }

  #G-Y fitting
  fitGY<-lm( Y ~  as.matrix( Gy )  )
  by<-as.numeric(summary(fitGY)$coef[-1,1]) ;  byse<-as.numeric(summary(fitGY)$coef[-1,2])


  ###B matrix #i.e. \int basisfunction(t)*eigenfunction(t) dt matrix  ##dim(B) #K L
  bbb<-c()
  for(l in 1:L){
    bbb<-cbind(bbb, res$workGrid^(l-1))  #polynomial basisfunction
  }
  bb<-bbb[-1,]

  #dim(bb); dim( res$phi[-1,1:K]  )

  B<-t(res$phi[-1,1:K])%*%bb*( res$workGrid[3] - res$workGrid[2]  ) #numeric integration



  if(Fit ){  #whether fit to get est/CIs, ggplot, MSE an coverage rate (based on true XYmodel)?
    ###Parametric fitting with the current L-------------------------------------------------
    ###--------------------------------------------------------------------------------------
    BBX<-BX%*%B #BBX: BX+B #Alpha%*%B

    S<-vcov(fitGY  )[-1,-1] #checked！ 和之前的对角线元素一样： Sigma<-diag(  byse^2 )
    Est<-solve(    t(BBX)%*%solve(S)%*%BBX )%*%t(BBX)%*%solve(S)%*%by  #确实一样的
    fitRES$Est<-Est
    tau2<-max(1,   t(by- BBX%*%Est)%*%solve( S )%*%(by- BBX%*%Est)/(J-length(Est))      )
    fitRES$tau<-sqrt(tau2)
    VM<-solve(    t(BBX)%*%solve(S)%*%BBX )*tau2  #estimator variance matrix
    fitRES$VM<-VM

    pointwise_shape_var<-diag(   bbb%*%VM%*%t( bbb ) )#取diagonal;易理解
    ggdata<-data.frame(time=res$workGrid, effect=bbb%*%Est,
                       effect_low=bbb%*%Est-1.96* sqrt(pointwise_shape_var  ),
                       effect_up=bbb%*%Est+1.96* sqrt(pointwise_shape_var  ),
                       true_shape=NaN)
    if(!is.na(XYmodel)){
      #if( !XYmodel%in%c('0','1','2','3','4','5','3.5') ){ stop('please use correct X-Y model') }

      if(XYmodel=='0'){
        fun<-function(t){  0*(t<Inf)   }
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
    plotdif<- max(ggdata$effect_up)-min(ggdata$effect_low)
    p<- ggplot(ggdata, aes(time, effect))+
      geom_hline(yintercept = 0,linewidth=0.5,linetype = 2,col='grey' )+
      geom_line(ggdata, mapping =aes(time, true_shape), alpha=1,linewidth=1,col='blue'  )+
      geom_line(ggdata, mapping =aes(time, effect), alpha=1,linewidth=1  )+
      geom_line(ggdata, mapping =aes(time, effect_low), alpha=1,linewidth=1,linetype = 2  )+
      geom_line(ggdata, mapping =aes(time, effect_up), alpha=1,linewidth=1,linetype = 2  )+
      labs(x='Age',y='Time-varying effect')+
      theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
      coord_cartesian( ylim = c(min(ggdata$effect_low)-plotdif  ,   max(ggdata$effect_up)+plotdif ) )

    fitRES$p<-p



    ###MSE and coverage rate (最好是返回一个vector over res$workGrid; 方便后续局部分析)
    if(!is.na(XYmodel)){
      SE<- ( bbb%*%Est - fun(  res$workGrid ) )^2 #Squared errors #mean(SE)->MSE
    }else{
      SE<-NA
    }
    fitRES$SE<-SE;fitRES$MSE<-mean(SE)


    ###coverage analysis
    if(!is.na(XYmodel)){
      Co<-abs( bbb%*%Est  - fun( res$workGrid)   )< 1.96* sqrt(pointwise_shape_var  )
      #cover-or-not vector #mean(Co)->coverage rate
    }else{
      Co<-NA  #no true shape, no coverage
    }
    fitRES$Co<-Co;fitRES$Coverage_rate<-mean(Co)

    ###significant timepoints information
    fitRES$time_points<-res$workGrid
    fitRES$sig_points<-as.numeric(   (ggdata$effect_up*ggdata$effect_low)>0   )*(2*as.numeric( ggdata$effect_low>0  )-1 )
    #前者返回 0或1 (显著与否)； 后者返回-1或1

  }



  ###testing--------------------------------------------------------------------------
  ###this is the advanced IVW Q test (also nonparametric MVMR test is just let L=nPC)-
  ###---------------------------------------------------------------------------------

  ##overlapping-sample instrument-phenotype estimates covariance (Q statistic中的Gamma vector)
  fitGY<-lm(Y~  as.matrix( Gy )  ); Yresiduals<-resid( fitGY )[IDmatch] #IDmatch 很关键； IDmatch: c(7,NA,4,NA,NA,...)
  Cov<-c()  #这个对所有SNPs都是共用的
  for(p in 1:nPC){  #只用前nPC个
    fitGX<-lm(PC[,p]~  as.matrix( G )  ); PCresiduals<-resid( fitGX )
    Cov<-c( Cov ,  cov( PCresiduals[!is.na(IDmatch)],Yresiduals[!is.na(IDmatch)] )      )   #complete overlap sample size so 1/N
  }
  Gam<-c()#Gam是个J*K的matrix
  for(j in 1:J){
    Gam<-rbind( Gam, sum(!is.na(IDmatch)  )/( as.numeric(nrow(G))*as.numeric(nrow(Gy))  )*1/( var(G[,j]) )*Cov)  #sum(!is.na(IDmatch)  )就是Ns
  }
  Gam[,is.na(Gam[1,])]<-0 #把NA的col的项换成0，以防Qfunction结果出问题； 注意，出现NA那么肯定整个一列都是NA


  #30/8/2023 updates:  also take into account the covariance of G-PCs covariance
  #we need Sigma: a K*K*J array, each j contains the covariance matrix of jSNP-PCs association estimators
  Cov_X<-matrix(NA,K,K)  #这个对所有SNPs都是共用的  #dim( Cov_X   ) #K K
  ref_tab<-cbind( rep(   1:K, each=K  )  ,  rep(1:K   )  )
  for(p in 1:(K*K)){  #只用前nPC个
    p_1<-ref_tab[p,1]; p_2<-ref_tab[p,2]
    fitGX1<-lm(PC[,p_1]~  as.matrix( G )  ); PCresiduals1<-resid( fitGX1 )
    fitGX2<-lm(PC[,p_2]~  as.matrix( G )  ); PCresiduals2<-resid( fitGX2 )
    Cov_X[p_1,p_2  ]<- cov( PCresiduals1,PCresiduals2   )
  }
  Sigma<-array(NA, dim=c(K,K,J))#Gam是个K*K*J的array
  for(j in 1:J){
    Sigma[,,j]<-1/as.numeric(nrow(G))*1/( var(G[,j]) )*Cov_X  #surely one sample; #Cov_X是共用的 #dim( Cov_X   ) #K K
  }

  ##iterative algorithm to get the final robust estimate and the Q final value
  #Qfunction_res<-Qfunction(  v=rep(0,L  ),by=by,byse=byse, B=B,BX=BX,BXse=BXse,Gam=Gam) #Gamma vector #initial input vector
  Qfunction_res<-Qfunction(  v=rep(0,L  ),by=by,byse=byse, B=B,BX=BX,Sigma=Sigma,Gam=Gam)
  v0<-Qfunction_res$Est
  iter_time<-0;total_dif<-Inf
  while( abs(total_dif)>1e-10  ){
    iter_time<-iter_time+1
    #v_updated<-Qfunction( v=v0,by=by, byse=byse,B=B,BX=BX,BXse=BXse,Gam=Gam)$Est
    v_updated<-Qfunction( v=v0,by=by, byse=byse,B=B,BX=BX,Sigma=Sigma,Gam=Gam)$Est
    total_dif<-sum( abs(v_updated-v0) )
    v0<-v_updated
  }
  #iter_time;v0
  #Qvalue<-Qfunction( v=v0,by=by, byse=byse,B=B,BX=BX,BXse=BXse,Gam=Gam)$Qvalue  #qchisq(0.95,J-L)
  Qvalue<-Qfunction( v=v0,by=by, byse=byse,B=B,BX=BX,Sigma=Sigma,Gam=Gam)$Qvalue  #qchisq(0.95,J-L)

  #15/sep/2023新增一个保护机制?：
  #如果迭代后的Q value比super naive的Q value(即，完全不考虑GXeffect的uncertainty和GX GY之间的covariance)还要大
  #那么干脆直接用super naice的Q value
  # if(    Qvalue >Qfunction_res$Qvalue ){
  #   Qvalue<-Qfunction_res$Qvalue
  #   print( 'warning: the iterative converged Q value is larger than the super naive Q value; may occur optimization problem!so we use the latter Q value' )
  #   }

  pvalue<-1-pchisq(Qvalue,J-L)

  fitRES$Qfunction_results<-Qfunction( v=v0,by=by, byse=byse,B=B,BX=BX,Sigma=Sigma,Gam=Gam) #it contains the robust inference results!
  #注意： fitRES$Qfunction_results用的是最终迭代完(即趋于稳定)的Q information；不一定是我们所需要的(可能super naive Q的值更小)
  fitRES$Qvalue<-Qvalue
  fitRES$Qdf<-J-L
  fitRES$pvalue<-pvalue

  return(fitRES)
}
#Parametric fitting test返回的Q test statistic就是一步到位的test (df = J-L)



#
# seed<-8;XYmodel_used<-'3'
# ZXmodel_used<-'B'
# #get res
# res<-readRDS( paste0('D:\\files\\R new\\MPCMR\\fPCA_res\\sim_scenario',ZXmodel_used,'_',seed,'.RData') )
# plotEifun(res)
#
# #get RES (seed track)
# set.seed(seed)
# RES<-getX(J=30,ZXmodel=ZXmodel_used) #和上方readRDS中文件名保持一致
# indi_plot(RES,123)
# DAT<-getY(RES, XYmodel=XYmodel_used)#DAT: complete data
#
#
# ##nonparametric MPCMR fitting--------------------------------------
# MPCMRFiting_res<-MPCMRfit( Gmatrix=RES$DAT[,1:RES$details$J],
#                            Yvector=DAT$Y,
#                            res=res,
#                            XYmodel=XYmodel_used )
# MPCMRFiting_res$ISres
# MPCMRFiting_res$p1
# MPCMRFiting_res$MSE
# MPCMRFiting_res$par_fit_res$p
# MPCMRFiting_res$par_fit_res$MSE
# MPCMRFiting_res$first_step_test#IV validity testing
#
# ##naive MPC fiting----------------------------------------------
# ##directly PC-Y regression---------------------------
# MPCMRFiting_res<-naiveMPCMRfit( Gmatrix=RES$DAT[,1:RES$details$J],
#                                 Yvector=DAT$Y,
#                                 res=res,
#                                 XYmodel=XYmodel_used )
# MPCMRFiting_res$p
# MPCMRFiting_res$MSE
#
# ##parametric MPCMR fitting-------------------------------------------
# ParMPCMRFiting_res<-ParMPCMRfit( Gmatrix=RES$DAT[,1:RES$details$J],
#                                  Yvector=DAT$Y,
#                                  res=res,
#                                  L=1,  #人为选择
#                                  XYmodel=XYmodel_used,
#                                  Fit=TRUE )
# ParMPCMRFiting_res$p
# ParMPCMRFiting_res$Qvalue
# ParMPCMRFiting_res$pvalue
# ParMPCMRFiting_res$Qdf
#
#
#
# #
# #许久之前的测试
# #
#
# ParMPCMRFiting_res<-ParMPCMRfit( Gmatrix=RES$DAT[,1:RES$details$J],
#                                  Yvector=DAT$Y,
#                                  res=res,
#                                  L=2,
#                                  XYmodel='3',
#                                  Fit=FALSE )
# QvalueL_2<-ParMPCMRFiting_res$Qvalue
# QdfL_2<-ParMPCMRFiting_res$Qdf
# ParMPCMRFiting_res<-ParMPCMRfit( Gmatrix=RES$DAT[,1:RES$details$J],
#                                  Yvector=DAT$Y,
#                                  res=res,
#                                  L=1,
#                                  XYmodel='3',
#                                  Fit=FALSE )
# QvalueL_1<-ParMPCMRFiting_res$Qvalue
# QdfL_1<-ParMPCMRFiting_res$Qdf
#
#
# QvalueL_1-QvalueL_2  ;qchisq(0.95,QdfL_1- QdfL_2  )
#
#
#
#
# #weak instrument analysis-----------------------------------------------------
# cFF<-c()
# QQ<-c()
# for(ZXmodel_used in c('A','B','C')){
#   cF<-c()#conditional F values
#   Q<-c()#Q p-values
#   for(ii in 1:50){
#     seed<-ii;XYmodel_used<-'1'#XYmodel_used不重要；IS都会是一样的
#     #get res
#     res<-readRDS( paste0('D:\\files\\R new\\MPCMR\\fPCA_res\\sim_scenario',ZXmodel_used,'_',seed,'.RData') )
#     #get RES (seed track)
#     set.seed(seed)
#     RES<-getX(J=30,ZXmodel=ZXmodel_used) #和上方readRDS中文件名保持一致
#     DAT<-getY(RES, XYmodel=XYmodel_used)#DAT: complete data
#     ##nonparametric MPCMR fitting--------------------------------------
#     MPCMRfit_res<-MPCMRfit( Gmatrix=RES$DAT[,1:RES$details$J],
#                             Yvector=DAT$Y,
#                             res=res,
#                             XYmodel=XYmodel_used )
#
#     cF<-cbind(cF,   as.vector( MPCMRfit_res$ISres[,3]) )
#     Q<-cbind(Q,   as.vector( MPCMRfit_res$ISres[,6]) )
#   }
#   cFF<-rbind(cFF,cF)
#   QQ<-rbind(QQ,round(Q,3))
# }
# View(cFF);View(QQ)
#
# #one-value weak instrument strength definition:   minimal cF value and maximal Q p-values
#
# cFF_min<-c()
# QQ_max<-c()
# for(i in 1:3){ # 'A' 'B' 'C' 共3种ZX model scenarios
#   cF_sub<-cFF[2*(i-1)+(1:2),]
#   Q_sub<-QQ[2*(i-1)+(1:2),]
#   cF_min<-apply(  cF_sub, 2, min )
#   Q_max<-apply(  Q_sub, 2, max )
#   cFF_min<-rbind(cFF_min,cF_min )
#   QQ_max<-rbind(QQ_max,Q_max )
# }
# View(cFF_min);View(QQ_max)
#
#
# ggdata<-data.frame(  Values=c( as.vector( t(cFF_min) ) ,  as.vector( t(QQ_max) )  ),
#                      ZXmodel= rep(c('A','B','C') , each=50  ) ,
#                      IStype=rep(c('Conditional F values' , 'Q p-values'  ) , each=3*50  ))
#
# p<- ggplot(ggdata, aes(Values))+
#   geom_histogram(aes(y = after_stat(ncount)),bins = 100) +
#   facet_grid( rows=vars(ZXmodel), cols=vars(IStype) , scales ='free' )+
#   theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
#   scale_x_continuous( n.breaks=10)
# p
#
# ggsave(paste0('IShist.eps' ),   #.eps
#        plot = p ,  #非指定aste0('C:\\Users\\Haodong Tian\\Desktop\\All_TVMR\\plots\\'),
#        height = 6, width = 8, units = "in",limitsize=TRUE)
#
#
# #not-reject (i.e. believed weak IS) rate under ZXmodel='A' (the known weak instrument case)
# sum(as.vector( t(QQ_max) )[1:50]>0.05)  #46
#
#





