###TVMR illustration

###this is for illustration; if you wish to reproduce the results in simulation; refer to pre_sim.R and formal_sim.R scripts

library(data.table)#as.data.table()
library(fdapace)#PACE
library(tidyverse) #一些data.table %>% 操作
library(corrplot)#correlated SNPs
library(parallel)

#f辅助functions: (mm_lm_onesample) (gmm_lm_twosample) (IS) (Qfunction)

#getDat -> plotEifun -> IVW: MPCMRfit     -> naiveMPCMRfit
#                    -> GMM: MPCMR_GMM    -> maiveMPCMR_GMM




###workflow/pipeline
# RES<-getX(J=30,ZXmodel='C')
# res <- FPCA(RES$Ly_sim, RES$Lt_sim,list(dataType='Sparse', error=TRUE, verbose=TRUE))#res: fPCAresults
# for(XYmodel_used in c('0','1','2','3','4','5')){
#   DAT<-getY(RES, XYmodel=XYmodel_used)#DAT: complete data
#   nPC<-sum(res$cumFVE<0.95)+1  #0.95 left to sensitivity analysis
#   fitres<-MPCMRfitting( Gmatrix=RES$DAT[,1:RES$details$J],
#                         PCmatrix=res$xiEst[,1:nPC],
#                         Yvector=DAT$Y,
#                         parfitting=TRUE,
#                         test=TRUE
#                         )
#   #MPCMRfitting: nonparametric fitting + polynomial parametric fitting + two-step test
# }






###先把res<-FPCA存储好;之后再分析
simFPCA<-function(seed){
  set.seed(seed)
  RES<-getX(J=30,ZXmodel='A')
  res <- FPCA(RES$Ly_sim, RES$Lt_sim,list(dataType='Sparse', error=TRUE, verbose=TRUE))
  saveRDS(res,file=paste0('D:\\files\\R new\\MPCMR\\fPCA_res\\sim_scenarioA_',seed,'.RData') )
  set.seed(seed)
  RES<-getX(J=30,ZXmodel='B')
  res <- FPCA(RES$Ly_sim, RES$Lt_sim,list(dataType='Sparse', error=TRUE, verbose=TRUE))
  saveRDS(res,file=paste0('D:\\files\\R new\\MPCMR\\fPCA_res\\sim_scenarioB_',seed,'.RData') )
  set.seed(seed)
  RES<-getX(J=30,ZXmodel='C')
  res <- FPCA(RES$Ly_sim, RES$Lt_sim,list(dataType='Sparse', error=TRUE, verbose=TRUE))
  saveRDS(res,file=paste0('D:\\files\\R new\\MPCMR\\fPCA_res\\sim_scenarioC_',seed,'.RData') )
  #res<-readRDS( paste0('D:\\files\\R new\\MPCMR\\fPCA_res\\sim_scenarioC_',seed,'.RData') )
  #write.csv(RES$DAT,paste0('D:\\files\\R new\\MPCMR\\fPCA_res\\sim_scenarioC_DAT_',seed,'.csv'), row.names=T)
  set.seed(seed)
  RES<-getX(J=30,ZXmodel='D')
  res <- FPCA(RES$Ly_sim, RES$Lt_sim,list(dataType='Sparse', error=TRUE, verbose=TRUE))
  saveRDS(res,file=paste0('D:\\files\\R new\\MPCMR\\fPCA_res\\sim_scenarioD_',seed,'.RData') )
}


simFPCA_onlyD<-function(seed){
  RES<-getX(J=30,ZXmodel='D')
  res <- FPCA(RES$Ly_sim, RES$Lt_sim,list(dataType='Sparse', error=TRUE, verbose=TRUE))
  saveRDS(res,file=paste0('D:\\files\\R new\\MPCMR\\fPCA_res\\sim_scenarioD_',seed,'.RData') )
}

library(parallel)
n_cores_used<-detectCores()-1
cl<-makeCluster(n_cores_used)
clusterEvalQ(cl=cl , expr=library(data.table))
clusterEvalQ(cl=cl , expr=library(fdapace))
clusterExport(  cl=cl ,  varlist=c('getX', 'simFPCA' )  )
Nb<-n_cores_used;Nb
parSapply(   cl ,  52:100, simFPCA )
parSapply(   cl ,  2:10, simFPCA_onlyD )
stopCluster(cl)




###Paper results:
#0. simple results (individual curves, fPCA fitting curve, eigenfunction curves, tec)
#1. MSE&Coverage table
#2. Sign plot (now abandon)
#3. Constant-effect test
#4. one-parameter (under the constant-effect assumption) fitting
#5. Weak instrument assessment
#6. IV validity assessment

#*only Constant-effect test not reject; then go for global causal null test
#* e.g. XY scenario3, naive G-Y association test may not powerful, but Constant-effect test can be of high power



###0.simple results (single fitted plot; Design plot, etc)------------------------------------
###-------------------------------------------------------------------------------------------

XYmodel_used<-'3' #用同一个固定的XYmodel作为演示

# 三种ZX model; 每种3种fitting methods (association vs MPCMR vs MPCMRp)
for(ZXmodel_used in c('A','B','C')){
  if(ZXmodel_used=='A'){seed<-1}
  if(ZXmodel_used=='B'){seed<-2}
  if(ZXmodel_used=='C'){seed<-3}
  res<-readRDS( paste0('D:\\files\\R new\\MPCMR\\fPCA_res\\sim_scenario',ZXmodel_used,'_',seed,'.RData') )
  set.seed(seed)
  RES<-getX(J=30,ZXmodel=ZXmodel_used)#和上方readRDS中文件名保持一致
  DAT<-getY(RES, XYmodel=XYmodel_used)#DAT: complete data
  ##naive MPC fitting----------------------------------------------
  MPCMRFiting_res<-naiveMPCMRfit( Gmatrix=RES$DAT[,1:RES$details$J],
                                  Yvector=DAT$Y,
                                  res=res,
                                  XYmodel=XYmodel_used )
  ggsave(paste0('association_',ZXmodel_used,'.eps' ),
         plot = MPCMRFiting_res$p ,
         path=paste0('C:\\Users\\Haodong Tian\\Desktop\\All_TVMR\\plots\\'),
         height = 4, width = 4.5, units = "in",limitsize=TRUE)


  ##nonparametrc  MPCMR fitting--------------------------------------
  MPCMRFiting_res<-MPCMRfit( Gmatrix=RES$DAT[,1:RES$details$J],
                             Yvector=DAT$Y,
                             res=res,
                             XYmodel=XYmodel_used )
  ggsave(paste0('MPCMR_',ZXmodel_used,'.eps' ),
         plot = MPCMRFiting_res$p1 ,
         path=paste0('C:\\Users\\Haodong Tian\\Desktop\\All_TVMR\\plots\\'),
         height = 4, width = 4.5, units = "in",limitsize=TRUE)

  #parametric fitting-----------------------------------------------
  ggsave(paste0('MPCMRp_',ZXmodel_used,'.eps' ),
         plot = MPCMRFiting_res$par_fit_res$p ,
         path=paste0('C:\\Users\\Haodong Tian\\Desktop\\All_TVMR\\plots\\'),
         height = 4, width = 4.5, units = "in",limitsize=TRUE)
}




###variable trajectory curve
set.seed(111);RES<-getX(J=30,ZXmodel='A') #和上方readRDS中文件名保持一致
indi_plot(RES,123)#图像生成模式是plot(...) 所以不方便储存
set.seed(222);RES<-getX(J=30,ZXmodel='B') #和上方readRDS中文件名保持一致
indi_plot(RES,123)
set.seed(333);RES<-getX(J=30,ZXmodel='C') #和上方readRDS中文件名保持一致
indi_plot(RES,123)
#ZXmodelA  600 400



###1. MSE&Coverge table---------------------------------------------------------------------
###-----------------------------------------------------------------------------------------
MSEvector<-c()
COVvector<-c()
ZXmodelvector<-rep( c('A','B','C') ,  each=6*3 )  #6 XYmodels #3 estimation types
XYmodelvector<-rep( c('0','1','2','3','4','5') , each=3  )  #3 estimation types
Type<-rep(  c('naive','nonpar','par') )

for(ii in c(1:50)){
  seed<-ii
  for(ZXmodel_used in c('A','B','C')){
    res<-readRDS( paste0('D:\\files\\R new\\MPCMR\\fPCA_res\\sim_scenario',ZXmodel_used,'_',seed,'.RData') )
    set.seed(seed)
    RES<-getX(J=30,ZXmodel=ZXmodel_used) #和上方readRDS中文件名保持一致
    for(XYmodel_used in c('0','1','2','3','4','5') ){
      #由于getY也具有随机性，最好前面也跟着set.seed(seed+1000)
      DAT<-getY(RES, XYmodel=XYmodel_used)#DAT: complete data
      ##naive MPC fitting----------------------------------------------
      MPCMRFiting_res<-naiveMPCMRfit( Gmatrix=RES$DAT[,1:RES$details$J],
                                      Yvector=DAT$Y,
                                      res=res,
                                      XYmodel=XYmodel_used )
      MSEvector<-c(MSEvector,MPCMRFiting_res$MSE)
      COVvector<-c(COVvector,MPCMRFiting_res$Coverage_rate)
      ##nonparametrc  MPCMR fitting--------------------------------------
      MPCMRFiting_res<-MPCMRfit( Gmatrix=RES$DAT[,1:RES$details$J],
                                 Yvector=DAT$Y,
                                 res=res,
                                 XYmodel=XYmodel_used )
      MSEvector<-c(MSEvector,MPCMRFiting_res$MSE)
      COVvector<-c(COVvector,MPCMRFiting_res$Coverage_rate)
      #parametric fitting-----------------------------------------------
      MSEvector<-c(MSEvector,MPCMRFiting_res$par_fit_res$MSE)
      COVvector<-c(COVvector,MPCMRFiting_res$par_fit_res$Coverage_rate)
    }
  }
}

resdat<-data.frame(   MSE=MSEvector, COV=COVvector,
                      ZXmodel= ZXmodelvector, XYmodel=XYmodelvector,Type=Type)
dim(resdat)
nrow(resdat)#  2700 = 50*6*4*3  (50 sim times)*(6 XYmodels)*(4 ZXmodels)*(3 estimation types)


TABLE<-c()
for(ZXmodel_used in c('A','B','C') ){
  resdat_sub<-resdat[resdat$ZXmodel==ZXmodel_used,]
  for(type_used in  c('naive','nonpar','par')){
    resdat_sub_sub<-resdat_sub[resdat_sub$Type==type_used,]
    MSEvec<-c()
    COVvec<-c()
    for(XYmodel_used in c('0','1','2','3','4','5') ){
      resdat_sub_sub_sub<-resdat_sub_sub[resdat_sub_sub$XYmodel==XYmodel_used,]
      #MSE
      MSEvec<-c(MSEvec, mean( resdat_sub_sub_sub$MSE ) )
      #COV
      COVvec<-c(COVvec , mean( resdat_sub_sub_sub$COV ) )
    }
    TABLE<-cbind(TABLE, MSEvec , COVvec)
  }
}
dim(TABLE) # 6 24=4*3*2
View(TABLE)

library(xtable)
xtable(TABLE*100,digits=3)
xtable(TABLE[,13:18]*100,digits=3)
resdat_sub<-resdat[ (resdat$ZXmodel=='A')&(resdat$XYmodel=='1')&(resdat$Type=='par'),]

hist(  resdat_sub$MSE ,n=50  )



###2. Sign plot-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------------------


#only focus the XYmodel='4' and '5'
for(ZXmodel_used in c('A','B','C')){

  for(XYmodel_used in c('4','5') ){

    Sig_points<-c()
    for(ii in c(1:50)){
      seed<-ii
      res<-readRDS( paste0('D:\\files\\R new\\MPCMR\\fPCA_res\\sim_scenario',ZXmodel_used,'_',seed,'.RData') )
      set.seed(seed)
      RES<-getX(J=30,ZXmodel=ZXmodel_used) #和上方readRDS中文件名保持一致
      DAT<-getY(RES, XYmodel=XYmodel_used)
      #MPCMR fit only
      MPCMRFiting_res<-MPCMRfit( Gmatrix=RES$DAT[,1:RES$details$J],
                                 Yvector=DAT$Y,
                                 res=res,
                                 XYmodel=XYmodel_used )
      Sig_points<-cbind(Sig_points,MPCMRFiting_res$sig_points)
    }
    #ggplot
    ggdata <- data.frame( Sim = rep(1:50, each=nrow(  Sig_points  ) ),
                          Age = MPCMRFiting_res$time_points,
                          Significance =  as.factor(  as.vector(Sig_points)    )   )
    p<-ggplot(ggdata, aes(Sim, Age)) +
      geom_tile(aes(fill = Significance))+
      theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
      coord_flip()+
      scale_fill_manual(values=c("-1"="blue","0"="grey", "1"="red"))

  }
}







###--------------------------------------------------------------------------------------------
###---------------------------------------------------------------------------------------------
#weak instrument analysis-----------------------------------------------------
cFF<-c()
QQ<-c()
for(ZXmodel_used in c('A','B','C')){
  cF<-c()#conditional F values
  Q<-c()#Q p-values
  for(ii in 1:50){
    seed<-ii;XYmodel_used<-'1'#XYmodel_used不重要；IS都会是一样的
    #get res
    res<-readRDS( paste0('D:\\files\\R new\\MPCMR\\fPCA_res\\sim_scenario',ZXmodel_used,'_',seed,'.RData') )
    #get RES (seed track)
    set.seed(seed)
    RES<-getX(J=30,ZXmodel=ZXmodel_used) #和上方readRDS中文件名保持一致
    DAT<-getY(RES, XYmodel=XYmodel_used)#DAT: complete data
    ##nonparametric MPCMR fitting--------------------------------------
    MPCMRfit_res<-MPCMRfit( Gmatrix=RES$DAT[,1:RES$details$J],
                            Yvector=DAT$Y,
                            res=res,
                            XYmodel=XYmodel_used )

    cF<-cbind(cF,   as.vector( MPCMRfit_res$ISres[,3]) )
    Q<-cbind(Q,   as.vector( MPCMRfit_res$ISres[,6]) )
  }
  cFF<-rbind(cFF,cF)
  QQ<-rbind(QQ,round(Q,3))
}
View(cFF);View(QQ)

#one-value weak instrument strength definition:   minimal cF value and maximal Q p-values

cFF_min<-c()
QQ_max<-c()
for(i in 1:3){ # 'A' 'B' 'C' 共3种ZX model scenarios
  cF_sub<-cFF[2*(i-1)+(1:2),]
  Q_sub<-QQ[2*(i-1)+(1:2),]
  cF_min<-apply(  cF_sub, 2, min )
  Q_max<-apply(  Q_sub, 2, max )
  cFF_min<-rbind(cFF_min,cF_min )
  QQ_max<-rbind(QQ_max,Q_max )
}
View(cFF_min);View(QQ_max)


ggdata<-data.frame(  Values=c( as.vector( t(cFF_min) ) ,  as.vector( t(QQ_max) )  ),
                     ZXmodel= rep(c('A','B','C') , each=50  ) ,
                     IStype=rep(c('Conditional F values' , 'Q p-values'  ) , each=3*50  ))

p<- ggplot(ggdata, aes(Values))+
  geom_histogram(aes(y = after_stat(ncount)),bins = 100) +
  facet_grid( rows=vars(ZXmodel), cols=vars(IStype) , scales ='free' )+
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  scale_x_continuous( n.breaks=10)
p

ggsave(paste0('IShist.eps' ),   #.eps
       plot = p ,  #非指定
       path=paste0('C:\\Users\\Haodong Tian\\Desktop\\All_TVMR\\plots\\'),
       height = 6, width = 8, units = "in",limitsize=TRUE)


#not-reject (i.e. believed weak IS) rate under ZXmodel='A' (the known weak instrument case)
sum(as.vector( t(QQ_max) )[1:50]>0.05)  #46



###3. Constant-effect test------------------------------------------------------------------
###-----------------------------------------------------------------------------------------

#array style
Parray<-array(NA, dim=c(15,3,50))  #50 就是fPCA有的数量

for(ii in c(1:50)){
  seed<-ii
  Pmatrix<-c()  #final dim(Pmatrix) #15 3
  for(ZXmodel_used in c('A','B','C')){
    Pvector<-c()
    res<-readRDS( paste0('D:\\files\\R new\\MPCMR\\fPCA_res\\sim_scenario',ZXmodel_used,'_',seed,'.RData') )
    set.seed(seed)
    RES<-getX(J=30,ZXmodel=ZXmodel_used) #和上方readRDS中文件名保持一致
    for(XYmodel_used in c('111','222','333') ){
      if(XYmodel_used=='111' ){ bcandidates<-c( 0,0.1, 0.5, 1, 10 ) }
      if(XYmodel_used=='222' ){ bcandidates<-c(0.001,0.01,0.05, 0.1 ,1 ) }
      if(XYmodel_used=='333' ){ bcandidates<-c(0.01,0.1,0.5,1, 10 ) }
      for(b_used in  bcandidates ){
        DAT<-getY(RES, XYmodel=XYmodel_used,b=b_used)#DAT: complete data
        ##parametric  MPCMR fitting with L=1--------------------------------------
        MPCMRFiting_res<-ParMPCMRfit( Gmatrix=RES$DAT[,1:RES$details$J],
                                      Yvector=DAT$Y,
                                      res=res,
                                      L=1,  #L=1 means time-varying effect is constant
                                      Fit=FALSE #fit 与否不重要
        )
        Pvector<-c(Pvector,MPCMRFiting_res$pvalue)
      }
    }
    Pmatrix<-cbind(Pmatrix, Pvector)
  }
  Parray[,,ii]<-Pmatrix
}

dim(Parray) # 15 3 50

Reject_array<-Parray<0.05

apply(Reject_array, c(1,2), mean   )

xtable( apply(Reject_array, c(1,2), mean   )    )



###one parameter fitting under the weak instrument (ZXmodel A)--------------------------
###-------------------------------------------------------------------------------------

ZXmodel_used<-'A'
MSEmatrix<-c();COVmatrix<-c()
CE_MSEmatrix<-c();CE_COVmatrix<-c()
for(ii in c(1:50)){
  seed<-ii
  res<-readRDS( paste0('D:\\files\\R new\\MPCMR\\fPCA_res\\sim_scenario',ZXmodel_used,'_',seed,'.RData') )
  set.seed(seed)
  RES<-getX(J=30,ZXmodel=ZXmodel_used) #和上方readRDS中文件名保持一致
  MSEvector<-c(); COVvector<-c(); CE_MSEvector<-c(); CE_COVvector<-c()
  for(XYmodel_used in c('0','1','2','3','4','5') ){#6种 XYmodels
    DAT<-getY(RES, XYmodel=XYmodel_used)#DAT: complete data
    ##One-parametric MPCMR fitting----------------------------------------------
    ParMPCMRfit_res<-ParMPCMRfit( Gmatrix=RES$DAT[,1:RES$details$J],
                                  Yvector=DAT$Y,
                                  res=res,
                                  L=1,  #人为选择L=1 due tO the weak instrument
                                  XYmodel=XYmodel_used,
                                  Fit=TRUE )
    ##global (curve) MSE and coverage rate--------------------------------------
    MSEvector<-c(MSEvector ,  ParMPCMRfit_res$MSE)  #这些MSE 是average over Gird timepoint; 这些都要再mean的
    COVvector<-c(COVvector ,  ParMPCMRfit_res$Coverage_rate)

    ##CE (cumulative effect) MSE and coverage rate-----------------------------------------
    est<-as.numeric(ParMPCMRfit_res$Est); se<-sqrt( as.numeric(ParMPCMRfit_res$VM)  )
    CEest<-50*est   ; CEse<-50*se
    if(XYmodel_used=='0'){CEtrue = 0  }
    if(XYmodel_used=='1'){CEtrue = 5  }
    if(XYmodel_used=='2'){CEtrue = 25  }
    if(XYmodel_used=='3'){CEtrue = 0  }
    if(XYmodel_used=='4'){CEtrue = 2  }
    if(XYmodel_used=='5'){CEtrue = 2  }
    CE_MSEvector<-c( CE_MSEvector ,  (CEest-CEtrue)^2 )
    CE_COVvector<-c( CE_COVvector ,  abs( CEest - CEtrue  )<1.96*CEse  )
  }
  MSEmatrix<-cbind(  MSEmatrix,MSEvector  )
  COVmatrix<-cbind(  COVmatrix,COVvector  )
  CE_MSEmatrix<-cbind(  CE_MSEmatrix,CE_MSEvector  )
  CE_COVmatrix<-cbind(  CE_COVmatrix,CE_COVvector  )
}
apply(  MSEmatrix, 1, mean )
apply(  COVmatrix, 1, mean )
apply(  CE_MSEmatrix, 1, mean )
apply(  CE_COVmatrix, 1, mean )

xtable(  data.frame(
  Scenario=paste0(  'A',0:5 )   , CE=c( 0,5,25,0,2,2) ,
  MSE=apply(  MSEmatrix, 1, mean )*100 ,Coverage=apply(  COVmatrix, 1, mean ),
  MSE2=apply(  CE_MSEmatrix, 1, mean ) ,Coverage2=apply(  CE_COVmatrix, 1, mean )    ), digits=3  )

###IV validity testing--------------------------------------------------------------------
###---------------------------------------------------------------------------------------

Rarray<-array(NA, dim=c(18,4,50))  #50 就是fPCA有的数量

for(ii in c(1:50)){
  seed<-ii
  Rmatrix<-c()  #Reject/Result matrix
  for(ZXmodel_used in c('A','B','C')){
    res<-readRDS( paste0('D:\\files\\R new\\MPCMR\\fPCA_res\\sim_scenario',ZXmodel_used,'_',seed,'.RData') )
    set.seed(seed)
    RES<-getX(J=30,ZXmodel=ZXmodel_used) #和上方readRDS中文件名保持一致
    for(XYmodel_used in c('0','1','2','3','4','5') ){
      Pvector<-c()#plriotrpic vector
      for(pleiN_used in c(0,5,10,30)){
        DAT<-getY(RES, XYmodel=XYmodel_used,plei=TRUE, pleiN=pleiN_used)#DAT: complete data

        ##nonparametrc  MPCMR fitting--------------------------------------
        MPCMRfit_res<-MPCMRfit( Gmatrix=RES$DAT[,1:RES$details$J],
                                Yvector=DAT$Y,
                                res=res,
                        XYmodel=XYmodel_used )
        Pvector<-c(Pvector,MPCMRfit_res$first_step_test$pvalue<0.05)
      }
      Rmatrix<-rbind(  Rmatrix, Pvector )
    }
  }
  Rarray[,,ii]<-Rmatrix
}

Results<-apply(Rarray, c(1,2), mean   )
xtable( Results, digits=3 )




