###MPCMR_GMM_twosample: fit the MPCMR with GMM
#(two-sample setting; one is from individual level; another can be from summary-level; GMM fitting: all summary-level)

#return: the IV strength results (only relevant to GX data) + parametric fitting curve
#(IV validity test is now removed); if use: simple use Qfunction



#该function逻辑是：只适用于indiivudal data：如果是one-sample 直接用用gmm_lm_onesample;
#                                           如果是overlapping data，直接prun成one-sample再gmm_lm_onesample；
#                                           如果是two-sample 直接用gmm_lm_onesample



get_true_shape_values<-function( workGrid, XYmodel    ){
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
  
  return(   fun( workGrid   )  )
  
}



MPCMR_GMM_twosample<-function(
    Gmatrix, #G matrix for PCs
    res, #res<-FPCA(...); (Gmatrix res) is the individual-level GX data
    by_used, #instrument-outcome association vector
    sy_used,#standard error vector of by  #(Gymatrix Yvector) is the indiivudal-level GY data for MPCMR_GMM; but now store as the summary information
    ny_used, #sample size for instrument-outcome regressions
    nPC=NA, #how many nPC to be retained? default is the the just-exceed-to-95%
    nL=NA, #how many polynomial basisfunctions used for fitting? only available for <= nPC; the default is = nPC
    LMCI=TRUE, #whether to calculate the CI with LM for nonparametric fit? this could be time-consuming
    LMCI2=TRUE, #whether to calculate the CI with LM for semiparametric fit? this could be time-consuming
    nLM=20, #how many grid used for LM CIs?
    Parallel=TRUE, #whether to use parallel? this is to fasten the CI calculation with LM
    XYmodel=NA #if labelled as '1' '2', ...; use this as the true effect (real data cannot have this)
){
  ###result list
  fitRES<-list()
  
  ###logic check and warnings
  if(is.na(nPC)){   nPC<-sum(res$cumFVE<0.95)+1    } #选用FVE 恰好大于95%的principal components
  if( nPC>ncol( res$xiEst ) ){stop('The exposure curves cannot support much information; please use smaller nPC')  }
  fitRES$nPC_used<-nPC #nPC_used就是nPC=K
  
  PCmatrix<-res$xiEst[,1:nPC];PCmatrix<-as.matrix(PCmatrix   ) #以防nrol=1的matrix退化成vector
  PC_<-PCmatrix[!is.na(  PCmatrix[,1] ),] #removing possible missing values(如果PACE内有null的list element,$xiEst会有missing data)
  PC_<-as.matrix(PC_   ) #以防nrol=1的matrix退化成vector
  
  if( nrow(PCmatrix)!=nrow(PC_)){
    cat('Note that there exists NA values for some individuals in fPCA results; check the measured timepoints for all individuals when running FPCA', '\n')}
  if( sum(is.na(Gmatrix))>0   ){stop('there exists missing data for the Gmatrix') }
  if(   nrow( Gmatrix)!=nrow(PC_  )    ){
    stop('the nrow of the Gmatrix is not equal to the nrow of the principal components matrix (after removing possible missing-PC individuals)')
  }
  
  
  is.md<-function(x){ return(   is.matrix(x)|is.data.frame(x)     )  } #is matrix  or data.frame?
  
  
  if( length(by_used)!=length(sy_used)   ){  stop('the G-Y assocation and its standard error are not of the same length')       }
  
  J<-ncol(Gmatrix)            #敲定J
  K<- ncol(PCmatrix)#nPC<-K   #敲定K
  
  
  if( length(by_used)!= J   ){ stop( 'be careful about the number of instruments used for the GY summary data'  ) }
  
  if(is.na(nL)){  L<- K  }else{ L<-nL }  #敲定L
  
  if(  L> K   ){stop('the number of basisfunctions should not be more than the number of exposures (scores)' )  }
  
  
  fitRES$L<-L
  fitRES$K<-K
  
  ###now start---------------------------------------------------------
  
  
  ###Weak instrument strength analysis (只与Z X data有关,和Y无关;所以不用担心overlapping data或者GY summary data)
  ###################################################################################################
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
  
  G<-Gmatrix
  
  #G-X fitting (only for *naive* Q-based weak instrument test)    *naive意味着不考虑PCs之间的相关性
  BX<-c(); BXse<-c()
  nPC<-K
  for(p in 1:nPC){  #只用前nPC个
    fitGX<-lm(PC_[,p]~  as.matrix( G )  )
    bx<-as.numeric(summary(fitGX)$coef[-1,1]) ;  bxse<-as.numeric(summary(fitGX)$coef[-1,2])
    BX<-cbind(BX, bx  );   BXse<-cbind( BXse, bxse )
  }
  
  
  #IV-PC1 IV-PC2 scatterplot
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
  
  ###IS: instrument Strength with Q style value (now updates with advanced Q test statistic)----------------------------
  if(nPC>=2){
    
    
    #Gam<-matrix(0,nrow=J,ncol=nPC-1) #Gam直接固定成a vector of 0，因为完美fPCA下，不同的phenotype (PCs) 之间是covariance-zero的 #错！
    IS_Q<-c() #dim(IS_Q) #nPC-1  3
    
    for(p in 1:nPC){ #每个PCs都跑一遍Q test #当前的这个k-th variable就是作为'outcome'的variable
      K_<-nPC-1 #K已经定义成nPC了
      #把第p个PC当作Y
      fitGY<-lm(PC_[,p]~  as.matrix( G )  ); Yresiduals<-resid( fitGY )
      by<-as.numeric(summary(fitGY)$coef[-1,1]) ;  byse<-as.numeric(summary(fitGY)$coef[-1,2])
      #把1:nPC减去第p个的对象们当作X
      PC__<-PC_[,  (1:nPC)[-p]  ]              ;PC__<-as.matrix( PC__ ) #以防matrix退化成vector
      nPC_<-nPC-1#nPC_<-K_
      #BX (i.e. G-X effect) calculation  BXse用后续的算法得到
      BX<-c();BXse<-c()
      for(p_ in 1:nPC_){  #只用前nPC个
        fitGX<-lm(PC__[,p_]~  as.matrix( G )  )
        bx<-as.numeric(summary(fitGX)$coef[-1,1]) ;  bxse<-as.numeric(summary(fitGX)$coef[-1,2])  #bxse和BXse其实不重要；可忽略
        BX<-cbind(BX, bx  ) ;   BXse<-cbind( BXse, bxse )
      }
      
      ##GX GY covariance---------
      Cov<-c()  #这个对所有SNPs都是共用的
      for(p_ in (1:nPC_)   ){  #只用前1:nPC减去第p个
        fitGX<-lm(PC__[,p_]~  as.matrix( G )  ); PCresiduals<-resid( fitGX )
        Cov<-c( Cov ,  cov( PCresiduals,Yresiduals )      )   #complete overlap sample size so 1/N
      }
      Gam<-c()#Gam是个J*(nPC-1)的matrix
      for(j in 1:J){
        Gam<-rbind( Gam, as.numeric(nrow(G))/( as.numeric(nrow(G))*as.numeric(nrow(G))  )*1/( var(G[,j]) )*Cov)  #sum(!is.na(IDmatch)  )就是Ns
      }
      Gam[,is.na(Gam[1,])]<-0 #把NA的col的项换成0，以防Qfunction结果出问题； 注意，出现NA那么肯定整个一列都是NA
      
      ##GX GX covariance-----------
      Cov_X<-matrix(NA,K_,K_)  #这个对所有SNPs都是共用的  #dim( Cov_X   ) #K_ K_
      ref_tab<-cbind( rep(   1:K_, each=K_  )  ,  rep(1:K_   )  )
      for(p in 1:(K_*K_)){  #只用前nPC个
        p_1<-ref_tab[p,1]; p_2<-ref_tab[p,2]
        fitGX1<-lm(PC__[,p_1]~  as.matrix( G )  ); PCresiduals1<-resid( fitGX1 )
        fitGX2<-lm(PC__[,p_2]~  as.matrix( G )  ); PCresiduals2<-resid( fitGX2 )
        Cov_X[p_1,p_2  ]<- cov( PCresiduals1,PCresiduals2   )
      }
      Sigma<-array(NA, dim=c(K_,K_,J))#Sigma是个K*K*J的array
      for(j in 1:J){
        Sigma[,,j]<-1/as.numeric(nrow(G))*1/( var(G[,j]) )*Cov_X  #surely one sample; #Cov_X是共用的 #dim( Cov_X   ) #K K
      }
      
      
      ##iterative algorithm to get the final robust estimate and the Q final value
      B<-diag(1, nPC-1  )#经典操作，B为I matrix表示最精简的Q statistics
      
      Qfunction_res<-Qfunction(  v=rep(0,nPC_  ),by=by,byse=byse, B=B,BX=BX,Sigma=Sigma,Gam=Gam)
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
      df<-J-(nPC-1)
      pvalue<-1-pchisq(Qvalue,J-(nPC-1))
      
      #result storage
      IS_Q<-rbind(IS_Q,  c(Qvalue ,df,  pvalue ) )
      
    }
    
    #storage into fitRES$ISres
    fitRES$ISres<-cbind( fitRES$ISres, IS_Q   )
    colnames(fitRES$ISres)<-c('RR','F','cF','Qvalue','df','pvalue') #pvalue越小越好，越应该拒绝！说明instrument strength没问题！
    
    
  }
  
  
  ###IS analysis end
  
  #############################################################################
  
  
  
  ###GMM fit (Ash's GMM is similar to MR-RAPS/GRAPPLE - already consider the uncertainty of the first stage)
  ##################################################################################################
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
  
  
  bx_used<-matrix( NA, J,K ); sx_used<-matrix( NA, J,K )
  for(j in 1:J){
    for(k in 1:K){
      GXfit<-lm( PC_[,k] ~ Gmatrix[,j] )
      bx_used[j,k]<-as.numeric(  summary( GXfit )$coef[2,1]  )
      sx_used[j,k]<-as.numeric(  summary( GXfit )$coef[2,2]  )
    }
  }
  nx_used<-nrow(Gmatrix)
  
  ld_used<-cor(Gmatrix)
  cor.x_used<-cor(PC_)
  
  
  
  gmm_res<-gmm_lm_twosample( bx=bx_used,by=by_used,sx=sx_used,sy=sy_used,
                             nx=nx_used,ny=ny_used,
                             ld=ld_used,cor.x=cor.x_used,
                             beta0=NA )         #GMM inference + Q statistics(这是GMM的Q,弃用) + LM statistics
  
  
  
  LM_CI_low<-NaN
  LM_CI_up<-NaN
  if( LMCI ){ #whether to calculate the CI with LM statistics?
    ##LM confidence interval------------
    #(直接算beta的CI region；然后套进functional curve，在Grid上算)
    #先准备candidate points
    nn<-nLM
    beta_candidates<-c()  #nn^K * K matrix
    for(k in 1:K){
      seqs<-seq( gmm_res$gmm_est[k] -4*gmm_res$gmm_se[k], gmm_res$gmm_est[k] +4*gmm_res$gmm_se[k],length=nn)
      beta_candidates<-cbind(  beta_candidates ,    rep(  seqs , each=nn^(  K-k ) ) )
    }
    vector_to_LM<-function(vector){ LMres<-getLM_twosample( bx=bx_used,by=by_used,sx=sx_used,sy=sy_used,
                                                            nx=nx_used,ny=ny_used,
                                                            ld=ld_used,cor.x=cor.x_used,
                                                            beta0=vector );return( LMres$lm_pval>0.05  )  }
    #TRUE or 1:代表着pvalue>0.05;即不拒绝H0;即在CI内部
    if(!Parallel){
      LMres_vector<-apply( beta_candidates, 1,  vector_to_LM   )#比较耗时，可以用parallel
    }else{
      #parallel
      if(detectCores()-1 <1    ){ cl<-makeCluster(1) }else{   cl<-makeCluster(detectCores()-1)   }
      clusterExport(  cl=cl ,  varlist=c( 'bx_used', 'by_used', 'sx_used' , 'sy_used' , 'nx_used' , 'ny_used' , 'ld_used' , 'cor.x_used',
                                          'getLM_twosample','vector_to_LM'),envir=environment()  )
      LMres_vector<-parApply(cl, beta_candidates, 1, vector_to_LM)
      stopCluster(cl)
      #matrix( LMres_vector, nn,nn  )
    }
    
    
    iterative_points<-c(  1, nn-1) #when K=1  #用来找到boundary point上LMrestult的LMres_vector对应位置
    if(K>=2){ for(k in 2:K){iterative_points<-c(iterative_points, nn^k- nn^(k-1) + iterative_points )}    }
    boundaryLMres<-cbind(beta_candidates[iterative_points,]   , LMres_vector[iterative_points]  )   #边缘坐标点上的LM情况；共2^K个values (最好都是0或者FALSE)
    colnames(boundaryLMres)<-c( paste0('expsoure', 1:K) , 'LMresult' )
    fitRES$boundaryLMres<-boundaryLMres
    
    LMres_used<-beta_candidates[LMres_vector,]#the CI inside points
    pointwise_LM_range<-(res$phi)[,1:nPC]%*%t( LMres_used  )  #dim(   t( LMres_used  )   )[1] #nPC #i.e. K
    LM_CI_low<-apply( pointwise_LM_range, 1, min  )
    LM_CI_up<-apply( pointwise_LM_range, 1, max  )
  }
  
  ###Q test results (i.e. IV validity test) (valid IV情况下不应该拒绝)
  #########################################################
  fitRES$IV_validity_test<-c(     'Q_statistic'= gmm_res$Q_stat ,'df'= J-K, 'p_value' = gmm_res$Q_pval  )
  
  
  ###semi-parametric (之前的nonparametric其实属于一种特定的semi-parametric)
  #point estimates
  MPCMRest<-gmm_res$gmm_est
  fitRES$MPCMRest<-MPCMRest
  #variance matrix using GMM for MPCMR
  MPCMRvar<-gmm_res$variance_matrix
  fitRES$MPCMRvar<-MPCMRvar
  
  #when basisfunctions are eigenfunctions---------------------------------
  pointwise_shape_var<-diag(   ( (res$phi)[,1:nPC]   )%*%MPCMRvar%*%t( (res$phi)[,1:nPC])  )#取diagonal;易理解
  ggdata<-data.frame(time=res$workGrid, effect=(res$phi)[,1:nPC]%*%MPCMRest,
                     effect_low=(res$phi)[,1:nPC]%*%MPCMRest-1.96* sqrt(pointwise_shape_var  ),
                     effect_up=(res$phi)[,1:nPC]%*%MPCMRest+1.96* sqrt(pointwise_shape_var  ),
                     LM_low=LM_CI_low,
                     LM_up=LM_CI_up,
                     true_shape=NaN)#true_shape=NaN 是为了没有已知XYmodel的情况
  if(!is.na(XYmodel)){
    #if( !XYmodel%in%c('0','1','2','3','4','5','3.5') ){ stop('please use correct X-Y model') }
    
    ggdata$true_shape<-get_true_shape_values( res$workGrid,XYmodel  )
    
  }
  
  fitRES$ggdata1<-ggdata
  plotdif<- max(ggdata$effect_up)-min(ggdata$effect_low)
  p1<- ggplot(ggdata, aes(time, effect))+
    geom_hline(yintercept = 0,linewidth=0.5,linetype = 2,col='grey' )+
    geom_line(ggdata, mapping =aes(time, true_shape), alpha=1,linewidth=1,col='blue'  )+
    geom_line(ggdata, mapping =aes(time, effect), alpha=1,linewidth=1  )+
    geom_line(ggdata, mapping =aes(time, effect_low), alpha=1,linewidth=1,linetype = 2  )+
    geom_line(ggdata, mapping =aes(time, effect_up), alpha=1,linewidth=1,linetype = 2  )+
    geom_line(ggdata, mapping =aes(time, LM_low), alpha=1,linewidth=1,linetype = 2 , col='#666666' )+
    geom_line(ggdata, mapping =aes(time, LM_up), alpha=1,linewidth=1,linetype = 2  , col='#666666' )+
    labs(x='Age',y='Time-varying effect')+
    theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    coord_cartesian( ylim = c(min(ggdata$effect_low)-0.5*plotdif  ,   max(ggdata$effect_up)+0.5*plotdif ) )
  
  fitRES$p1<-p1 #Removed 51 rows containing missing values (`geom_line()`). #不用担心 这些都是NaN
  
  ###MSE and COV results
  if(!is.na(XYmodel)){
    SE<- ( (res$phi)[,1:nPC]%*%MPCMRest -  ggdata$true_shape )^2 #Squared errors #mean(SE)->MSE
  }else{
    SE<-NA
  }
  fitRES$SE<-SE;fitRES$MSE<-mean(SE)
  
  ###coverage analysis
  if(!is.na(XYmodel)){
    Co<-abs( (res$phi)[,1:nPC]%*%MPCMRest  - ggdata$true_shape  )< 1.96* sqrt(pointwise_shape_var  )  #GMM的CIs;并不是LM的！
    Co_LM<-( ggdata$LM_up > ggdata$true_shape   )&(   ggdata$LM_low < ggdata$true_shape     ) #LM的coverage
    #covaer-or-not vector #mean(Co)->coverage rate
  }else{
    Co<-NA  #no true shape, no coverage
    Co_LM<-NA
  }
  fitRES$Co<-Co;fitRES$Coverage_rate<-mean(Co)
  fitRES$Co_LM<-Co_LM;fitRES$Coverage_rate_LM<-mean(Co_LM)
  
  ###significant timepoints information
  fitRES$time_points<-res$workGrid
  fitRES$sig_points<-as.numeric(   (ggdata$effect_up*ggdata$effect_low)>0   )*(2*as.numeric( ggdata$effect_low>0  )-1 )
  fitRES$sig_points_LM<-as.numeric(   (ggdata$LM_up*ggdata$LM_low)>0   )*(2*as.numeric( ggdata$LM_low>0  )-1 )
  #前者返回 0或1 (显著与否)； 后者返回-1或1
  
  
  #when basisfunctions are polynomial---------------------------------
  #-------------------------------------------------------------------
  #-------------------------------------------------------------------
  
  ###B matrix #i.e. \int basisfunction(t)*eigenfunction(t) dt matrix  ##dim(B) #K L #used for transforming exposures
  bbb<-c()
  for(l in 1:L){
    bbb<-cbind(bbb, res$workGrid^(l-1))  #polynomial basisfunction
  }
  bb<-bbb[-1,]
  #dim(bb); dim( res$phi[-1,1:K]  )
  
  B<-t(res$phi[-1,1:K])%*%bb*( res$workGrid[3] - res$workGrid[2]  ) #numeric integration
  #dim(B) #K L
  
  PC_transformed<-PC_%*%B  #transformed PC
  
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
  
  
  bx_used<-matrix( NA, J,K ); sx_used<-matrix( NA, J,K )
  for(j in 1:J){
    for(k in 1:K){
      GXfit<-lm( PC_transformed[,k] ~ Gmatrix[,j] )
      bx_used[j,k]<-as.numeric(  summary( GXfit )$coef[2,1]  )
      sx_used[j,k]<-as.numeric(  summary( GXfit )$coef[2,2]  )
    }
  }
  nx_used<-nrow(Gmatrix)
  
  ld_used<-cor(Gmatrix)
  cor.x_used<-cor(PC_transformed)
  
  
  
  gmm_res<-gmm_lm_twosample( bx=bx_used,by=by_used,sx=sx_used,sy=sy_used,
                             nx=nx_used,ny=ny_used,
                             ld=ld_used,cor.x=cor.x_used,
                             beta0=NA )         #GMM inference + Q statistics(这是GMM的Q,弃用) + LM statistics
  
  
  
  ###Q test results (i.e. IV validity + basisfunction assumption test) (valid IV + correct basiafunction assumption情况下不应该拒绝)
  #########################################################
  fitRES$IV_validity_and_basisfunction_test<-c(     'Q_statistic'= gmm_res$Q_stat ,'df'= J-L, 'p_value' = gmm_res$Q_pval  )
  
  
  LM_CI_low<-NaN
  LM_CI_up<-NaN
  if( LMCI2 ){ #whether to calculate the CI with LM statistics?
    ##LM confidence interval------------
    #(直接算beta的CI region；然后套进functional curve，在Grid上算)
    #先准备candidate points
    nn<-nLM
    beta_candidates<-c()  #nn^L * L matrix
    for(l in 1:L){
      seqs<-seq( gmm_res$gmm_est[l] -4*gmm_res$gmm_se[l], gmm_res$gmm_est[l] +4*gmm_res$gmm_se[l],length=nn)
      beta_candidates<-cbind(  beta_candidates ,    rep(  seqs , each=nn^(  L-l ) ) )
    }
    vector_to_LM_p<-function(vector){ LMres<-getLM_twosample( bx=bx_used,by=by_used,sx=sx_used,sy=sy_used,
                                                              nx=nx_used,ny=ny_used,
                                                              ld=ld_used,cor.x=cor.x_used,
                                                              beta0=vector );return( LMres$lm_pval>0.05  )  }
    #TRUE or 1:代表着pvalue>0.05;即不拒绝H0;即在CI内部
    if(!Parallel){
      LMres_vector<-apply( beta_candidates, 1,  vector_to_LM_p   )
    }else{
      #parallel
      if(detectCores()-1 <1    ){ cl<-makeCluster(1) }else{   cl<-makeCluster(detectCores()-1)   }
      clusterExport(  cl=cl ,  varlist=c( 'bx_used', 'by_used', 'sx_used' , 'sy_used' , 'nx_used' , 'ny_used' , 'ld_used' , 'cor.x_used',
                                          'getLM_twosample','vector_to_LM_p'),envir=environment()  )
      LMres_vector<-parApply(cl, beta_candidates, 1, vector_to_LM_p)
      stopCluster(cl)
      #matrix( LMres_vector, nn,nn  )
    }
    
    
    iterative_points<-c(  1, nn-1) #when L=1  #用来找到boundary point上LMrestult的LMres_vector对应位置
    if(L>=2){ for(l in 2:L){iterative_points<-c(iterative_points, nn^l- nn^(l-1) + iterative_points )}    }
    boundaryLMres<-cbind(beta_candidates[iterative_points,]   , LMres_vector[iterative_points]  )   #边缘坐标点上的LM情况；共2^K个values (最好都是0或者FALSE)
    colnames(boundaryLMres)<-c( paste0('transformed_expsoure', 1:L) , 'LMresult' )
    fitRES$boundaryLMres_p<-boundaryLMres
    
    LMres_used<-beta_candidates[LMres_vector,]#the CI inside points
    pointwise_LM_range<-bbb%*%t( LMres_used  )  #dim(   t( LMres_used  )   )[1]  #L
    LM_CI_low<-apply( pointwise_LM_range, 1, min  )
    LM_CI_up<-apply( pointwise_LM_range, 1, max  )
  }
  
  Est<-gmm_res$gmm_est
  fitRES$MPCMRest_p<-Est
  #variance matrix using GMM for MPCMR
  VM<-gmm_res$variance_matrix
  fitRES$MPCMRvar_p<-VM
  
  #ggplot
  pointwise_shape_var<-diag(   bbb%*%VM%*%t( bbb ) )#取diagonal;易理解
  ggdata<-data.frame(time=res$workGrid, effect=bbb%*%Est,
                     effect_low=bbb%*%Est-1.96* sqrt(pointwise_shape_var  ),
                     effect_up=bbb%*%Est+1.96* sqrt(pointwise_shape_var  ),
                     LM_low=LM_CI_low,
                     LM_up=LM_CI_up,
                     true_shape=NaN)
  if(!is.na(XYmodel)){
    #if( !XYmodel%in%c('0','1','2','3','4','5','3.5') ){ stop('please use correct X-Y model') }
    
    ggdata$true_shape<-get_true_shape_values( res$workGrid,XYmodel  )
    
  }
  
  fitRES$ggdata2<-ggdata
  
  plotdif<- max(ggdata$effect_up)-min(ggdata$effect_low)
  p2<- ggplot(ggdata, aes(time, effect))+
    geom_hline(yintercept = 0,linewidth=0.5,linetype = 2,col='grey' )+
    geom_line(ggdata, mapping =aes(time, true_shape), alpha=1,linewidth=1,col='blue'  )+
    geom_line(ggdata, mapping =aes(time, effect), alpha=1,linewidth=1  )+
    geom_line(ggdata, mapping =aes(time, effect_low), alpha=1,linewidth=1,linetype = 2  )+
    geom_line(ggdata, mapping =aes(time, effect_up), alpha=1,linewidth=1,linetype = 2  )+
    geom_line(ggdata, mapping =aes(time, LM_low), alpha=1,linewidth=1,linetype = 2 , col='#666666' )+
    geom_line(ggdata, mapping =aes(time, LM_up), alpha=1,linewidth=1,linetype = 2  , col='#666666' )+
    labs(x='Age',y='Time-varying effect')+
    theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    coord_cartesian( ylim = c(min(ggdata$effect_low)-0.5*plotdif  ,   max(ggdata$effect_up)+0.5*plotdif ) )
  
  fitRES$p2<-p2 #Removed 51 rows containing missing values (`geom_line()`). #不用担心 这些都是NaN
  
  ###MSE and coverage rate (最好是返回一个vector over res$workGrid; 方便后续局部分析)
  if(!is.na(XYmodel)){
    SE<- ( bbb%*%Est -ggdata$true_shape )^2 #Squared errors #mean(SE)->MSE
  }else{
    SE<-NA
  }
  fitRES$SE_p<-SE;fitRES$MSE_p<-mean(SE)
  
  
  ###coverage analysis
  if(!is.na(XYmodel)){
    Co<-abs( bbb%*%Est  - ggdata$true_shape )< 1.96* sqrt(pointwise_shape_var  )
    Co_LM<-( ggdata$LM_up > ggdata$true_shape   )&(   ggdata$LM_low < ggdata$true_shape     )
    #cover-or-not vector #mean(Co)->coverage rate
  }else{
    Co<-NA  #no true shape, no coverage
    Co_LM<-NA
  }
  fitRES$Co_p<-Co;fitRES$Coverage_rate_p<-mean(Co)
  fitRES$Co_p_LM<-Co_LM;fitRES$Coverage_rate_p_LM<-mean(Co_LM)
  
  ###significant timepoints information
  fitRES$sig_points_p<-as.numeric(   (ggdata$effect_up*ggdata$effect_low)>0   )*(2*as.numeric( ggdata$effect_low>0  )-1 )
  fitRES$sig_points_p_LM<-as.numeric(   (ggdata$LM_up*ggdata$LM_low)>0   )*(2*as.numeric( ggdata$LM_low>0  )-1 )
  #前者返回 0或1 (显著与否)； 后者返回-1或1
  
  return(fitRES)
  
}
#$L $K $ISres $scatterp $one_sample_size
#$IV_validity_test                     $MPCMRest $MPCMRvar $p1 $boundaryLMres $ggdata1 $SE $MSE $Co $Coverage_rate $Co_LM $Coverage_rate_LM $timepoints $sig_points $sig_points_LM
#$IV_validity_and_basisfunction_test $MPCMRest_p $MPCMRvar_p $p2 $boundaryLMres_p $ggdata2 $SE_p $MSE_p $Co_p $Coverage_rate_p $Co_p_LM $Coverage_rate_p_LM  $sig_points_p $sig_points_p_LM




# ###example
# 
# 
# 
# 
# seed<-20
# ZXmodel_used<-'D';XYmodel_used<-'2'
# N_used<-50000
# #get res
# res<-readRDS( paste0('D:\\files\\R new\\MPCMR\\fPCA_res\\sim_scenario',ZXmodel_used,'_',seed,'.RData') )
# plotEifun(res)
# 
# #or FPCA but with more sample size
# set.seed(seed)
# RES<-getX(N=N_used,J=30,ZXmodel=ZXmodel_used)
# res <- FPCA(RES$Ly_sim, RES$Lt_sim,list(dataType='Sparse', error=TRUE, verbose=TRUE))  #10 timepoints per individual: takes 5 mins
# saveRDS(res,file=paste0('C:\\Users\\Haodong\\Desktop\\sim_scenarioD_',seed,'_50000.RData') )
# res<-readRDS( paste0('C:\\Users\\Haodong\\Desktop\\sim_scenarioD_',seed,'_50000.RData') )
# 
# 
# #get the original data (note the seed track)
# set.seed(seed)
# RES<-getX(N=N_used,J=30,ZXmodel=ZXmodel_used) #和上方readRDS中 文件名保持一致
# indi_plot(RES,123)
# DAT<-getY(RES, XYmodel=XYmodel_used)
# 
# #get another samples for GY data
# set.seed(seed+10000)
# RES2<-getX(N=N_used,J=30,ZXmodel=ZXmodel_used, MGX_used=RES$details$MGX) #使用之前同一个original data中的genetic effect matrix
# indi_plot(RES2,123)
# DAT2<-getY(RES2, XYmodel=XYmodel_used)
# 
# 
# #and prepare the summary information of the GY data
# by<-rep(NA,RES2$details$J); sy<-rep(NA,RES2$details$J)
# for(j in 1:RES2$details$J){
#   GYfit<-lm(  DAT2$Y ~  DAT2[,j] )
#   by[j]<-as.numeric(  summary( GYfit )$coef[2,1]  )
#   sy[j]<-as.numeric(  summary( GYfit )$coef[2,2]  )
# }
# ny<-nrow(DAT2)
# 
# Gmatrix=RES$DAT[,1:RES$details$J];res=res;by_used=by;sy_used=sy;ny_used=ny;nPC=NA;nL=NA;LMCI=TRUE;LMCI2=TRUE;nLM=20; Parallel=TRUE;XYmodel=XYmodel_used
# 
# 
# 
# #MPCMR
# #start time
# cat('begin time :' ,Sys.time(), '\n')
# time1<-as.numeric(Sys.time())
# 
# 
# MPCMRres<-MPCMR_GMM_twosample(  Gmatrix=RES$DAT[,1:RES$details$J],
#                                 res=res,                                  #Gmatrix res 来自于同一个indiivudual data set
#                                 by_used=by, #instrument-outcome association vector
#                                 sy_used=sy, #standard error vector of by  #(Gymatrix Yvector) is the indiivudal-level GY data for MPCMR_GMM; but now store as the summary information
#                                 ny_used=ny, #sample size for instrument-outcome regressions
#                                 nPC=NA, #how many nPC to be retained? default is the the just-exceed-to-95%
#                                 nL=NA, #how many polynomial basisfunctions used for fitting? only available for <= nPC; the default is = nPC
#                                 LMCI=TRUE, #whether to calculate the CI with LM for nonparametric fit? this could be time-consuming
#                                 LMCI2=TRUE, #whether to calculate the CI with LM for semiparametric fit? this could be time-consumin
#                                 nLM=20, #how many grid used for LM CIs?
#                                 Parallel=TRUE, #whether to use parallel? this is to fasten the CI calculation with LM
#                                 XYmodel=XYmodel_used #if labelled as '1' '2', ...; use this as the true effect (real data cannot have this)
# )
# 
# 
# #record the time - end
# cat('finish time :' ,Sys.time(), '\n')
# time2<-as.numeric(Sys.time())
# ime1)/60, '(mins)', '\n')


