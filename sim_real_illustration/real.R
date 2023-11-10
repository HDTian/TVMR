
###new real application


library(data.table)#as.data.table()
library(fdapace)#PACE
library(tidyverse) #一些data.table %>% 操作
library(corrplot)#correlated SNPs
library(parallel)

master<-read.csv('D:\\files\\R new\\MPCMR\\master.csv')

Dat<-read.csv(paste0('D:\\files\\R new\\MPCMR\\sim_data_new\\Dat.csv'))

DatY<-read.csv(paste0('D:\\files\\R new\\MPCMR\\sim_data_new\\DatY.csv'))

get_exposure_and_time<-function(s){#string item
  v<-as.numeric(strsplit(  s, '_'   )[[1]])
  res<-list(Ly=v[ 1: (length(v)/2) ] , Lt=v[ -(1: (length(v)/2)) ]  )
  return(res)
}


#leftrange<-25   ; rightrange<-50
leftrange<-50   ; rightrange<-75
#leftrange<-0   ; rightrange<-1000
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

##先移除out-of-region项，合成一个Dat_sub,再做fPCA

Dat_sub<-Dat[ID,]  #dim(Dat_sub) 40951   460

Lt_real_sub<-Lt_real[ID] ; Ly_real_sub<-Ly_real[ID]#length( Lt_real_sub  ); length( Ly_real_sub  )   40951  #58411

##fPCA---------------
#suggest using existent package
#res <- FPCA(Ly_real_sub, Lt_real_sub,list(dataType='Sparse', error=TRUE, verbose=TRUE))#默认的就是dataType='Sparse'
#saveRDS(res,file=paste0("D:\\files\\R new\\MPCMR\\sim_data_new\\",leftrange,"_",rightrange,".RData"))
res<-readRDS(paste0("D:\\files\\R new\\MPCMR\\sim_data_new\\",leftrange,"_",rightrange,".RData"))

plotEifun(res)#越往后的PC对应的eigenfunction越浮夸(波动多)- 易理解

##individual fitted curve vs measured values over the selected time region
PC<-res$xiEst #dim(res$xiEst)  #40951     4
iindex<-54#随便选的
plot(  Lt_real_sub[[iindex]],Ly_real_sub[[iindex]] ,
       xlim=c(leftrange,rightrange) , ylim=c(80,160), col='red',main=iindex  )
lines( res$workGrid,     (res$phi)%*%PC[iindex,]+res$mu   )



###MPCMR fitting------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------

#review: 258SNPs
nrow(Dat_sub)==nrow(res$xiEst) #TRUE
dim(Dat_sub) #58411   460 (for 50_75)


#outcome time region selection
leftrangeY<-70;rightrangeY<-100
DatY_sub<-DatY[   (leftrangeY<DatY$age)&(DatY$age<rightrangeY)   ,   ]

dim(DatY_sub)

#59695   461 (70_100)  (urea)

#IDmatch
IDmatch_used<-match(Dat_sub$V1, DatY_sub$V1   )  #以Gmatrix为主导

sum(  !is.na(IDmatch_used)  )  #19075 #the overlapping sample size



#MPCMR_GMM
MPCMR_GMM_real<-MPCMR_GMM(    Gmatrix=Dat_sub[,2:259],  #除了让MPCMR_GMM自己做整合外，也可以在之前就做好one-sample整合
                              res=res,
                              Yvector=DatY_sub$Y,
                              Gymatrix=DatY_sub[,2:259],
                              IDmatch=IDmatch_used,
                              #nPC=NA,
                              nL=1,
                              LMCI=FALSE, 
                              LMCI2=FALSE,
                              #nLM=20, #with parallel
                              #Parallel=TRUE, #默认使用detectCores()-1 cores
                              XYmodel= NA,
                            )
##MPCMR_GMM各部分耗时：
#IS: 5 mins   gmm_lm_onesample: 1 min   LMCI2 (L=1 case): 5 mins

#Warning messages: Removed 51 rows containing missing values (`geom_line()`). 不重要 NaN的操作而已 

###result analysis 

##weak instrument analysis
MPCMR_GMM_real$ISres
#             RR        F        cF   Qvalue  df    pvalue
#PC1 0.040215706 9.444240 6.5522023 634.5986 256 0.0000000
#PC2 0.007861667 1.786024 1.1033493 251.2169 256 0.5726799
#PC3 0.005094029 1.154049 0.9691316 229.5926 256 0.8810106


##fit plot
MPCMR_GMM_real$p1;MPCMR_GMM_real$p2

#can reproduce the ggplot via ggdata

saveRDS(MPCMR_GMM_real,file='C:\\Users\\Haodong\\Desktop\\MPCMR_GMM_real.RData') #479 MB
MPCMR_GMM_real<-readRDS("C:\\Users\\Haodong\\Desktop\\MPCMR_GMM_real.RData")



ggdata1<-MPCMR_GMM_real$ggdata1
ggdata2<-MPCMR_GMM_real$ggdata2
ggdata<-rbind(ggdata1, ggdata2  )
ggdata$Scenario<-rep( c('Eigenfunction', 'Polynomial'), each=nrow(ggdata1)  )

#plotdif<- max(ggdata$effect_up)-min(ggdata$effect_low)
scales_y <- list(
  `Eigenfunction` = scale_y_continuous(limits = c(-0.05, 0.04)),
  `Polynomial` = scale_y_continuous(limits = c(-0.001, 0.002))
)
library(facetscales)
p<- ggplot(ggdata, aes(time, effect))+
  geom_hline(yintercept = 0,linewidth=0.5,linetype = 2,col='grey' )+
  geom_line(ggdata, mapping =aes(time, true_shape), alpha=1,linewidth=1,col='blue'  )+
  geom_line(ggdata, mapping =aes(time, effect), alpha=1,linewidth=1  )+
  geom_line(ggdata, mapping =aes(time, effect_low), alpha=1,linewidth=1,linetype = 2  )+
  geom_line(ggdata, mapping =aes(time, effect_up), alpha=1,linewidth=1,linetype = 2  )+
  geom_line(ggdata, mapping =aes(time, LM_low), alpha=1,linewidth=1,linetype = 2 , col='#666666' )+
  geom_line(ggdata, mapping =aes(time, LM_up), alpha=1,linewidth=1,linetype = 2  , col='#666666' )+
  facet_grid_sc(rows=vars(Scenario) , scales = list(y = scales_y))+
  #facet_grid(rows=vars(Scenario),scales='free'  )+
  labs(x='Age',y='Time-varying effect')+
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())#+coord_cartesian( ylim = c(min(ggdata$effect_low)-0.5*plotdif  ,   max(ggdata$effect_up)+0.5*plotdif ) )

p
#Removed 51 rows containing missing values (`geom_line()`).  不要紧 (true shape的NaN values)
ggsave(paste0('real_fit','.eps' ),
       plot = p,
       path=paste0('C:\\Users\\Haodong\\Desktop\\TVMR\\plots\\'),
       height = 6, width = 6, units = "in",limitsize=TRUE)


##constant effect (time-invariant effect) estimation 
fitRES$MPCMRest_p  #constent estimate: 0.0007404146
as.numeric(sqrt(gmm_res$variance_matrix)) #s.e. #0.0002110514
fitRES$MPCMRest_p + c(-1,1)*1.96*as.numeric(sqrt(gmm_res$variance_matrix)) #CI:0.0003267539 0.0011540754

c(  fitRES$ggdata2$LM_low[1],   fitRES$ggdata2$LM_up[1]  )#LM CI: 0.0001628002 0.0013180290

#cumulative effect (i.e. *25)

fitRES$MPCMRest_p*25 #0.01851037
(fitRES$MPCMRest_p + c(-1,1)*1.96*as.numeric(sqrt(gmm_res$variance_matrix)) )*25 #0.008168846 0.028851884
c(  fitRES$ggdata2$LM_low[1],   fitRES$ggdata2$LM_up[1]  )*25 # 0.004070006 0.032950724











