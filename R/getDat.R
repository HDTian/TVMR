


#getDat in MPCMR project


#both the complete exposure and the sparse exposure
getX<-function(N=10000, #sample size
               J=30,# number of SNPs (three-valued SNP)
               ZXmodel='A',
               MGX_used=NA, #if user define MGX, then will use this MGX; this is important for two-sample case where certain SNP should have the same effect
               nSparse=10 #the number of sparse measured timepoints for each individual #之前是20
){
  NT<-1000  #total number of points

  #
  TT<-50 #end timepoint

  # J<-30 # # of SNPs


  times<-seq(0, TT, len=(NT+1  )    );Times<-times[-1]
  times_squ<-times^2

  if( !ZXmodel%in%c('A','B','C','D','E','F') ){ stop('please use correct Z-X model') }

  if(ZXmodel=='A'){#constant effect
    a<-runif( J, -0.1,0.1 )
    b<-runif( J, -0,0)
    c<-runif( J, -0,0 )
    MGX<-t(   t(    matrix( rep(times_squ[-1],J), NT,J )%*%diag( c  )  )+
                t(    matrix( rep(times[-1],J), NT,J )%*%diag( b  )  )+
                a  )
  }

  if(ZXmodel=='B'){#linear effect
    a<-runif( J, -0.1,0.1 )
    b<-runif( J, -0.004 , 0.004)
    c<-runif( J, -0,0 )
    MGX<-t(   t(    matrix( rep(times_squ[-1],J), NT,J )%*%diag( c  )  )+
                t(    matrix( rep(times[-1],J), NT,J )%*%diag( b  )  )+
                a  )
  }

  if(ZXmodel=='C'){
    a<-runif( J , 0,50 )
    b<-runif( J, -0.1,0.1 )
    MGX<-c()
    for(j in 1:J){
      MGX<-cbind( MGX,  b[j]*(times>a[j])*0^{1+(-1)^j}+b[j]*(times<a[j])*0^{1+(-1)^(j+1)} )
    }
    MGX<-MGX[-1,]
  }

  if(ZXmodel=='D'){ #sin effect
    a<-runif( J, -1 , 1 ) #a=0 means constant effect  #之前是a<-runif( J, -1 , 1 )
    b<-runif( J, -0.1,0.1 )
    MGX<-t(   t(    0.05*sin( matrix( rep(times[-1],J), NT,J )%*%diag( a  )  )  )+  b ) #一直是0.05
  }
  # dim( MGX )  #1000   30  #the 30 SNPs effect at 1000 timepoints

  if(ZXmodel=='E'){ #sin effect
    a<-runif( J, -0.1 , 0.1 ) #a=0 means constant effect  #之前是a<-runif( J, -1 , 1 )
    b<-runif( J, -0.1,0.1 )
    MGX<-t(   t(    0.05*sin( matrix( rep(times[-1],J), NT,J )%*%diag( a  )  )  )+  b ) #一直是0.05
  }
  # dim( MGX )  #1000   30  #the 30 SNPs effect at 1000 timepoints

  if(ZXmodel=='F'){
    a<-runif( J, -0.1,0.1 )
    b<-runif( J, -0.002 , 0.002)
    #b<-c(runif( J-5, -0,0), runif( 5,-0.004 , 0.004) )
    c<-runif( J, -0,0 )
    MGX<-t(   t(    matrix( rep(times_squ[-1],J), NT,J )%*%diag( c  )  )+
                t(    matrix( rep(times[-1],J), NT,J )%*%diag( b  )  )+
                a  )
  }



  if( !is.na(  as.matrix(MGX_used)[1,1] )  ){  MGX<-MGX_used  }  #如果指定了MGX，那么就用它

  G<-matrix(  rbinom(J*N,2,0.3) , N , J  ) #J=30

  MBM<-matrix(     rnorm(  NT * N, 0, sqrt(  ( 1 ) /NT)    ) ,N, NT     ) #Matrix for Brownian motion
  UU<-apply(   MBM, 1 , cumsum) #保留第一个维度(即行保留)；注意最后apply的形式可能需要转秩  #dim(UU) #1000 N
  UUU<-rnorm(N,0,1)  + t(UU)#注意 加了一个baseline N(0,1^2)的confounding value
  MEX<-matrix(     rnorm(  NT * N, 0, sqrt(  ( 1 ) /NT)    ) ,N, NT     ) #Matrix for Epsilon_X Brownian motion
  EX<-apply(   MEX, 1 , cumsum)  #dim(EX) #1000 N

  X<- G%*%t(MGX)  +  UUU +t(EX)

  #Y分开产生；为了更好的节省fPCA时间 in simulation
  #Y<- (X%*% effect_vec   )*TT/NT  + UUU[,NT] + rnorm(N,0,1  )

  DAT<-    cbind( G, X[ ,(1:50)*NT/50  ]    )#不包含Y
  DAT<-as.data.frame(DAT)
  names(DAT)<-c(paste0( 'G', 1:J ), paste0( 'X',1:50 ))
  ###Sparse
  ###
  Ly_sim<-list(); Lt_sim<-list()  #Ly_sim 就是指的X的值 at the sparse measured timepoints

  for(i in 1:nrow(  X  ) ){
    index_sparse<-(1:length(Times))[  sort(sample(1:length(Times),nSparse)) ] #nParse=20
    time_sparse<-Times[index_sparse]
    Ly_sim[[i]]<-X[i,index_sparse ]
    Lt_sim[[i]]<-time_sparse
  }

  ###result storage
  RES<-list()
  RES$DAT<-DAT
  RES$details<-list( J=J,
                     X=X,TT=TT,NT=NT,UUU=UUU, #used for further-generate Y
                     G=G,times=times,MGX=MGX)  #extra details for generating individual curve
  RES$Ly_sim<-Ly_sim
  RES$Lt_sim<-Lt_sim

  return(RES)
}
### $DAT: (SNP Xs)-dataframe  $details: elements to generate the Y $Ly_sim: sparse exposures list $Lt_sim: measured timepoints list
### examples:
#


#further get Y
#注意：getY依旧具有随机性，体现在epsilon_Y上
getY<-function(RES, #RES<-getX(...)
               XYmodel='1', #X-Y model scenarios
               b=0,  #the parameter user determined for XYmodel='000' '111' '222' '333'
               plei=FALSE, #whether simulate the outcome with the pleiptrppic effects?
               pleiN=0 #the number of invalid SNPs (only works when plei=TRUE)
){

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

  if(XYmodel=='000'){ #null XY effect
    fun<-function(t){  0*(t<Inf)   }
  }
  if(XYmodel=='111'){ #constant effect  #'000' is a special case of '111'
    fun<-function(t){  b*(t<Inf)   }   #b~0.1
  }
  if(XYmodel=='222'){ #linear effect
    fun<-function(t){  b*t    }     #b~0.02
  }
  if(XYmodel=='333'){ #sin effect
    fun<-function(t){  b*sin(0.1*t)   }  #b~0.1
  }
  # if(XYmodel=='333'){ #threshold effect
  #   fun<-function(t){  b*(t>30)   }  #b~0.1
  # }
  if(XYmodel=='444'){ #continuous threshold effect
    fun<-function(t){  b*(t-30)*(t>30)   }  #b~0.1
  }




  #elements recover
  X<-RES$details$X;TT<-RES$details$TT;NT<-RES$details$NT
  UUU<-RES$details$UUU

  times<-seq(0, TT, len=(NT+1  )    )
  effect_vec<-as.numeric(fun(  times[-1] ))

  Y<- (X%*% effect_vec   )*TT/NT  + 10*UUU[,NT] + rnorm(nrow(X),0,1  )


  ###pleiotropic effect case----------------
  if(plei){
    G<-RES$details$G  #dim(G) 10000 30
    pleiEffect<-  c(rep(0.5, pleiN  ),rep(0, 30-pleiN  ) )  #pleiotropic effect vector
    Y=Y+  G%*%pleiEffect  #length(Y) 10000
  }

  DAT<-RES$DAT
  DAT$Y<-as.numeric(Y)

  return(DAT)
}
#return DAT: the complete data.set (SNPs Xs Y)




#get the complete data
#这个函数用处不大；主要还是靠getX() FPCA() getY() 之后再构造MPC data.frame for MR fitting
getDat<-function(N=10000, #sample size
                 J=30,# number of SNPs (three-valued SNP)
                 ZXmodel='A', # scenario A B C #A: constant B:linear C:threshold
                 nParse=20,#the number of parse measured timepoints for each individual
                 XYmodel='1' # scenario 1 2 3 #0: null effect 1:constant 2:linear 3:threshold
){
  res1<-getX(N=N,J=J,ZXmodel=ZXmodel,nParse=nParse)
  res2<-getY(res1,XYmodel=XYmodel)


  ###result storage
  RES<-list()
  RES$DAT<-res2$DAT
  #RES$details<-list(G=G,times=times, X=X,MGX=MGX,UUU=UUU) #MGX: genetic-effect-matrix #UUU: confounding
  return(RES)  #the complete measured dataset
}


###draw an individual variable (Xt,Zt,Ut) plot over time in simulation
indi_plot<-function(RES, #RES<-getX()
                    i=1, #which individual
                    main=FALSE #if show the main (individual index)
){
  G<-RES$details$G; times<-RES$details$times; X<-RES$details$X; MGX<-RES$details$MGX; UUU<-RES$details$UUU
  if( !main ){
    plot( times[-1],   X[i,],type='l',ylim=c(min(X[i,])-1,max(0, max(X[i,])+1)   ),
          xlab='Time',ylab='Variable level',lwd=1.5 )
  }else{
    plot( times[-1],   X[i,],type='l',ylim=c(min(X[i,])-1,max(0, max(X[i,])+1)   ),
          main=paste0('Individual: ',i),xlab='Time',ylab='Variable level',lwd=1.5 )
  }
  abline( 0,0 ,lty=3,col='#999999')
  lines(  times[-1],    (G%*%t(MGX))[i,], lty=1,col='green'  )
  lines(  times[-1],   UUU[i,], lty=1,col='red'  )
  legend('topright', legend=c('Exposure', 'Gene score','Confouding'),
         lty=c(1, 1,1), col=c('black','green','red'),lwd=c(2,2,2),cex=0.8)
}



###examples:
# indi_plot(getX(),123)
#
# RES<-getX(J=30,ZXmodel='D')
#
# indi_plot(RES,223) #Zt Xt Ut curve
# DAT<-getY(RES,XYmodel='2' ) #complete data

#SNP specific effect
#for(i in 1:30){   plot( 1:1000 , RES$details$MGX[,i] ,type='l',main=i  )     }

#
# RES<-getX(J=35,ZXmodel='C')
#
# #这个函数单独运行
# res <- FPCA(RES$Ly_sim, RES$Lt_sim,list(dataType='Sparse', error=TRUE, verbose=TRUE))#很费时
#
# res$cumFVE




