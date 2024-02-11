#formal_real


library(data.table)#as.data.table()
library(fdapace)#PACE
library(tidyverse) #一些data.table %>% 操作
library(corrplot)#correlated SNPs
library(parallel)



Dat<-read.csv(paste0('/Users/haodongtian/Documents/files/R new/MPCMR/sim_data_new2/Dat.csv'))
Dat_n<-read.csv(paste0('/Users/haodongtian/Documents/files/R new/MPCMR/sim_data_new2/Dat_n.csv'))
DatY<-read.csv(paste0('/Users/haodongtian/Documents/files/R new/MPCMR/sim_data_new2/DatY.csv'))

###clarification
#Dat1: reference data (no any initial selection)
#Dat2: genotype data for SBP (258 nearly independent SNPs)
#Dat3: longitudinal SBP data (保留一个小数点)
#Dat4: urea data with measured timepoints

#Dat_n是Dat1 2的结合：即 （G X）data （for non-time-varying UVMR）
#DatY是Dat1 2 4的结合： 即 （G Y）data

#Dat是Dat1 2 3的结合：即 （G PCs）data

#Dat + DatY for MPCMR fitting (either one-sample individual level or two-sample summary level)



###UVMR---------------------------------------------------------------
dim(Dat_n)#367283    260  # ID + 258 SNPs + SBP(measured at study baseline) ; ID is called 'V1' ;SBP is called 'V2'
dim(DatY)#350223    461  # ID + 258SNPs + urea + timepoint + nuisance phenotypes

#individual-level data
#can use either the one-sample fitting (2SLS or IVW; IVW==2SLS) or overlapping.two--data fitting (IVW)
#their results should be very similar

#only use naive IVW (or 2SLS) fitting: do not consider the uncertainty of the G-X association estimates


#combine Dat_n and DatY into one-sample individual data



ID<-match(  Dat_n$V1, DatY$V1 )  #按照前者的顺序
DatY_<-DatY[ID,]
Dat_uvmr<-cbind(  Dat_n,  DatY_[,260]  )  #names(DatY_)[260]=="Y"
dim(Dat_uvmr) #367283    261
ID<-  apply(  is.na(Dat_uvmr[,1:261]  )  ,  1, sum  )==0  #remove NA individuals
Dat_UVMR<-Dat_uvmr[ID,]
dim(Dat_UVMR) #349883    261 # ID + 258 SNPs + SBP (V2) + urea (Y)
names( Dat_UVMR)[261]<-'Y'

#release some store spaces in the global environment
rm( Dat_uvmr, Dat12, Dat123, Dat124 )


###naive IVW---------------
fitGX<-lm( Dat_UVMR$V2 ~  as.matrix( Dat_UVMR[,2:259])   )
bx<-as.numeric(summary(fitGX)$coef[-1,1]) ;  bxse<-as.numeric(summary(fitGX)$coef[-1,2])
BX<-bx; BXse<-bxse #UVMR  so one-dimensional
fitGY<-lm( Dat_UVMR$Y ~  as.matrix( Dat_UVMR[,2:259] )  )
by<-as.numeric(summary(fitGY)$coef[-1,1]) ;  byse<-as.numeric(summary(fitGY)$coef[-1,2])


###my WLR(weigthed linear regression) code
J<-258 # number of SNPs
S<-vcov(fitGY  )[-1,-1]  #variance-covariance matrix
#point estimates
UVMRest<-solve(    t(BX)%*%solve( S )%*%BX            )%*%t(BX)%*%solve(S)%*%by
##random-effect setting
tau2<-max(1,   t(by- BX%*%UVMRest)%*%solve( S )%*%(by- BX%*%UVMRest)/(J-length(UVMRest))      )

#variance matrix of the IVW-regression estimator
UVMRvar<-solve(    t(BX)%*%solve( S )%*%BX            )*tau2

#estimate and 95% CI
as.numeric(UVMRest);  c(UVMRest)+c(-1,1)*1.96*sqrt( c(UVMRvar) )

# 0.00815089
# 0.004409627 0.011892153 (urea)

###Steve's package
library(MendelianRandomization)
mr_ivw(mr_input(bx, bxse, by, byse))

#Inverse-variance weighted method
#(variants uncorrelated, random-effect model)

#Number of Variants : 258

#------------------------------------------------------------------
#  Method Estimate Std Error 95% CI       p-value
#IVW    0.008     0.002 0.004, 0.012   0.000
#------------------------------------------------------------------
#  Residual standard error =  2.175
#Heterogeneity test statistic (Cochran's Q) = 1215.9872 on 257 degrees of freedom, (p-value = 0.0000). I^2 = 78.9%.
#F statistic = 28.6.


###MR scatterplot
BX<-cbind(bx,by)
BXse<-cbind(bxse,byse)
ggdata<-data.frame(   beta1=BX[,1], beta2=BX[,2],
                      beta1_error_up=BX[,1]+1.96*BXse[,1], beta1_error_low=BX[,1]-1.96*BXse[,1] ,
                      beta2_error_up=BX[,2]+1.96*BXse[,2], beta2_error_low=BX[,2]-1.96*BXse[,2])
scatterp<-     ggplot(ggdata, aes(beta1, beta2))+
  geom_errorbar(data=ggdata,mapping=aes(x=beta1,ymin=beta2_error_low,ymax=beta2_error_up),width = 0,color='#666666')+
  geom_errorbar(data=ggdata,mapping=aes(y=beta2,xmin=beta1_error_low,xmax=beta1_error_up),width = 0,color='#666666')+
  geom_point(data=ggdata,mapping=aes(x=beta1,y=beta2),size=1.8)+
  geom_hline(yintercept = 0,linetype='dashed',col='#999999')+
  geom_vline(xintercept = 0,linetype='dashed',col='#999999')+
  #geom_abline(intercept=0,slope=1,linetype='dashed',col='#999999')+
  labs(x = 'Genetic association with the exposure', y = 'Genetic association with the outcome')+
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
scatterp
#save as: UVMR_scatterplot_urea

ggsave(paste0('UVMR_scatterplot_urea_new','.eps' ),
       plot = scatterp,
       path=paste0('/Users/haodongtian/Documents/TVMR(latex)/plots/'),
       height = 8, width = 8, units = "in",limitsize=TRUE)




###MPCMR fitting--------------------------------------------------------------------------------
###---------------------------------------------------------------------------------------------


dim(Dat)  #116378   460  # ID + 258SNPs + Xt(the string for longitudinal information)  + nuisance phenotypes
dim(DatY) #350223   461  # ID + 258SNPs +         urea  + timepoint                    + nuisance phenotypes



###use the partial time region------------------------------
get_exposure_and_time<-function(s){#string item
  v<-as.numeric(strsplit(  s, '_'   )[[1]])
  res<-list(Ly=v[ 1: (length(v)/2) ] , Lt=v[ -(1: (length(v)/2)) ]  )
  return(res)
}


#the maximal time region
leftrange<-0   ; rightrange<-1000
Ly_real<-list()
Lt_real<-list()
ID<-rep(NA, nrow(Dat))  #用来记录那些individual没有落入需要的time region里，方便排除
for(i in 1:nrow(Dat)){
  res<-get_exposure_and_time( Dat$Xt[i])
  ts<-res$Lt  ; ys<-res$Ly
  #肯定是存在至少一个SBP measurement for all individual
  which_in_sub<-(leftrange<=ts)&(ts<=rightrange)
  Ly_real[[i]]<-ys[which_in_sub]
  Lt_real[[i]]<-ts[which_in_sub]
  ID[i]<-  sum(which_in_sub)>0
}
CreateDesignPlot(Lt_real) #很完美的design plot
#design_plot_0_1000  #500 500


#a propriate sub time region such that FCPA can well work (a time gap = 25是最大的时间范围同时使得design plot很好)
leftrange<-50   ; rightrange<-75
Ly_real<-list()
Lt_real<-list()
ID<-rep(NA, nrow(Dat))  #用来记录那些individual没有落入需要的time region里，方便排除
for(i in 1:nrow(Dat)){
  res<-get_exposure_and_time( Dat$Xt[i])
  ts<-res$Lt  ; ys<-res$Ly
  #肯定是存在至少一个SBP measurement for all individual
  which_in_sub<-(leftrange<=ts)&(ts<=rightrange)
  Ly_real[[i]]<-ys[which_in_sub]
  Lt_real[[i]]<-ts[which_in_sub]
  ID[i]<-  sum(which_in_sub)>0
}
CreateDesignPlot(Lt_real) #很完美的design plot
#design_plot_50_75  #500 500


##先移除out-of-region的NA项，合成一个Dat_sub,再做fPCA

Dat_sub<-Dat[ID,]
dim(Dat_sub) #106498   460

Lt_real_sub<-Lt_real[ID] ; Ly_real_sub<-Ly_real[ID]
length( Lt_real_sub  ); length( Ly_real_sub  )  # 106498  #106498

###FPCA----------------------------------------------------------------------------------

#memory release: only maintain Dat, DatY, Dat_sub, Lt_real_sub, Ly_real_sub
all_variables <- ls()
variables_to_keep <- c("Dat", "DatY" , "Dat_sub" , "Lt_real_sub" , "Ly_real_sub")
variables_to_remove <- setdiff(all_variables, variables_to_keep)
rm(list = variables_to_remove)

###FPCA
res <- FPCA(Ly_real_sub, Lt_real_sub,list(dataType='Sparse', error=TRUE, verbose=TRUE))#默认的就是dataType='Sparse'
#time cost: ～2 mins



saveRDS(res,file=paste0("/Users/haodongtian/Documents/files/R new/MPCMR/sim_data_new2/",50,"_",75,".RData"))
#res<-readRDS(paste0("/Users/haodongtian/Documents/files/R new/MPCMR/sim_data_new2/50_75.RData"))

###review:
#res$selectK: K the number of PC remained
#res$lambda: the eigenvalues  (K = length(res$lambda) )
#res$cumFVE: the cumulative FVE ( res$cumFVE ~= res$lambda/sum(res$lambda )  )

#res$workGrid: length fixed to 51; used for present the values of the eigenfunction
#res$obsGrid: all observed timepoint from the original data (50.0 50.1 50.2 50.3 ... 74.8 74.9 75.0)

#res$xiEst: individual estimated PCs matrix ( res$xiVar: a list of length n, each of which indicates the variance of the estimated PCs)
#res$phi: the 51*K matrix: the eigenfunction values at workGrid

#quick see the values eg: sum(res$phi[,2]*res$phi[,3])
#t(res$phi)%*%(res$phi)  #非对角线元素都等于0  good!


##eigenfunctions
plotEifun(res) #易理解: 越往后的PC对应的eigenfunction越浮夸(波动多)
#the first 4 PCs just explain 95% FVE
#Fig9_new  700*500


##individual fitted curve vs measured values over the selected time region
PC<-res$xiEst
dim(res$xiEst)  #106498     6

plot_individual_exposure_curve<-function( iindex ){
  plot(  Lt_real_sub[[iindex]],Ly_real_sub[[iindex]] ,
         xlim=c(50,75) , ylim=c(80,160), col='red',main=iindex , xlab='Age' , ylab='SBP' )
  lines( res$workGrid,     (res$phi)%*%PC[iindex,]+res$mu   )
}
plot_individual_exposure_curve(1)  #1 = iidex 随便选的 358也行



###MPCMR fitting---------------------------------------------------------------

#review: 258SNPs
nrow(Dat_sub)==nrow(res$xiEst) #TRUE #确保(G X)data都做了fPCA
dim(Dat_sub) #106498 (for 50_75)   460  #ID + 258SNPs + Xt(the string for longitudinal information)  + nuisance phenotypes


#outcome time region selection (maybe no restriction at all to increae the sample size)
leftrangeY<-0;rightrangeY<-100
DatY_sub<-DatY[   (leftrangeY<DatY$age)&(DatY$age<rightrangeY)   ,   ]
dim(DatY_sub) #350223    461


###store some data, which can be directly used for MPCMR fitting - this is convenient for HPC
write.csv(Dat_sub  , paste0('/Users/haodongtian/Documents/files/R new/MPCMR/sim_data_new2/Dat_sub.csv'), row.names=F)
write.csv(DatY_sub  , paste0('/Users/haodongtian/Documents/files/R new/MPCMR/sim_data_new2/DatY_sub.csv'), row.names=F)

Dat_sub<-read.csv(paste0('/Users/haodongtian/Documents/files/R new/MPCMR/sim_data_new2/Dat_sub.csv'))
DatY_sub<-read.csv(paste0('/Users/haodongtian/Documents/files/R new/MPCMR/sim_data_new2/DatY_sub.csv'))

###IDmatch (for one-sample individual data fitting)
IDmatch_used<-match(Dat_sub$V1, DatY_sub$V1   )  #以Gmatrix为主导

sum(  !is.na(IDmatch_used)  )  #101279 ##the overlapping sample size



###MPCMR_GMM (one-sample individual-level data fitting)




#basisfunction = the first two eigenfunctions (equivalent to just use the first 2 PCs; ie nPC=2)
print(    Sys.time()  )
MPCMR_GMM_real<-MPCMR_GMM(    Gmatrix=Dat_sub[,2:259],
                              res=res,
                              Yvector=DatY_sub$Y,
                              Gymatrix=DatY_sub[,2:259],
                              IDmatch=IDmatch_used,
                              nPC=2,  ##how many nPC to be retained?
                              nL=NA,  ##how many polynomial basisfunctions used for fitting? only acailable for <= nPC; the default is = nPC
                              eigenfit=TRUE,
                              polyfit=FALSE,  #no basisfunction fitting to save time
                              LMCI=TRUE,
                              LMCI2=TRUE,
                              nLM=10, #with parallel
                              Parallel=TRUE, #默认=TRUE 且使用detectCores()-1 cores
                              cores_used=5,
                              XYmodel= NA
                         )
print(    Sys.time()  )

saveRDS(MPCMR_GMM_real,file=paste0("/Users/haodongtian/Documents/files/R new/MPCMR/sim_data_new2/MPCMR_GMM_real_nLM10.RData"))


#Running large datasets in parallel can quickly get you into trouble.
#If you run out of memory the system will either crash or run incredibly slow.

###running time summary:
#instrument strength analysis: ~3 mins
#one-time gmm_lm_onesample() fit: ~2 mins
#one-time vector_to_LM() fit without parallel: ~3 mins
#one-time vector_to_LM() fit for each core under the parallel: ~9 mins
#if use more cores in parallel (with the same fixed RAM), the laptop will become slower





###Note: if MPCMR fitting is too time-consuming, consider using HPC to fit
###store the fitted results as MPCMR_GMM_real1.RData

#MPCMR_GMM_real<-readRDS(paste0("/Users/haodongtian/Documents/files/R new/MPCMR/sim_data_new2/MPCMR_GMM_real1.RData"))




###result analysis

##weak instrument analysis
MPCMR_GMM_real$ISres
#             RR         F       cF   Qvalue  df       pvalue
#PC1 0.031308178 13.308724 9.455552 943.0482 257 0.000000e+00
#PC2 0.005159924  2.135769 1.517416 383.1561 257 5.358544e-07




MPCMR_GMM_real_linear<-MPCMR_GMM(    Gmatrix=Dat_sub[,2:259],
                                     res=res,
                                     Yvector=DatY_sub$Y,
                                     Gymatrix=DatY_sub[,2:259],
                                     IDmatch=IDmatch_used,
                                     nPC=4,  ##how many nPC to be retained?
                                     nL=2,  ##how many polynomial basisfunctions used for fitting? only acailable for <= nPC; the default is = nPC
                                     eigenfit=FALSE,
                                     polyfit=TRUE,  #no basisfunction fitting to save time
                                     LMCI=TRUE,
                                     LMCI2=TRUE,
                                     nLM=10, #with parallel
                                     Parallel=TRUE, #默认=TRUE 且使用detectCores()-1 cores
                                     cores_used=5,
                                     XYmodel= NA
)

saveRDS(MPCMR_GMM_real_linear,file=paste0("/Users/haodongtian/Documents/files/R new/MPCMR/sim_data_new2/MPCMR_GMM_real_linear_nLM10.RData"))


MPCMR_GMM_real_linear$ISres
#             RR         F       cF   Qvalue  df       pvalue
#PC1 0.031308178 13.308724 7.847500 727.6312 255 0.000000e+00
#PC2 0.005159924  2.135769 1.418863 389.3271 255 1.178888e-07
#PC3 0.002881737  1.190069 1.011032 283.6672 255 1.049225e-01
#PC4 0.002659176  1.097912 1.027578 281.5470 255 1.217432e-01




#plot
ggdata1<-MPCMR_GMM_real$ggdata1
ggdata2<-MPCMR_GMM_real_linear$ggdata2

ggdata1$basis<-'eigenfunction'
ggdata2$basis<-'polynomial'

ggdata<-rbind( ggdata1 , ggdata2  )

plotdif<- max(ggdata$effect_up)-min(ggdata$effect_low)
p<- ggplot(ggdata, aes(time, effect))+
  geom_hline(yintercept = 0,linewidth=0.5,linetype = 2,col='grey' )+
  geom_line(ggdata, mapping =aes(time, effect), alpha=1,linewidth=1  )+
  geom_line(ggdata, mapping =aes(time, effect_low), alpha=1,linewidth=1,linetype = 2  )+
  geom_line(ggdata, mapping =aes(time, effect_up), alpha=1,linewidth=1,linetype = 2  )+
  geom_line(ggdata, mapping =aes(time, LM_low), alpha=1,linewidth=1,linetype = 2 , col='#666666' )+
  geom_line(ggdata, mapping =aes(time, LM_up), alpha=1,linewidth=1,linetype = 2  , col='#666666' )+
  labs(x='Age',y='Time-varying effect')+
  facet_grid(cols=vars(basis)  )+
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  coord_cartesian( ylim = c(min(ggdata$effect_low)-0.5*plotdif  ,   max(ggdata$effect_up)+0.5*plotdif ) )
p

ggsave(paste0('real_fit_eigen_poly','.eps' ),
       plot = p,
       path=paste0('/Users/haodongtian/Documents/TVMR(latex)/plots/'),
       height = 6, width = 12, units = "in",limitsize=TRUE)





