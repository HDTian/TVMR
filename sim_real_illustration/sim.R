
###new simulation script

###sim


library(data.table)#as.data.table()
library(fdapace)#PACE
library(tidyverse) #一些data.table %>% 操作
library(corrplot)#correlated SNPs
library(parallel)

###0: preparation (simulate FPCA and store)

###先把res<-FPCA存储好;之后再分析
###注意FPCA只和Z-X有关，和X-Y无关
simFPCA<-function(seed){#10 timepoints per individual: total takes 15 mins
  ###A
  set.seed(seed)
  RES<-getX(J=30,ZXmodel='A')
  res <- FPCA(RES$Ly_sim, RES$Lt_sim,list(dataType='Sparse', error=TRUE, verbose=TRUE))  #10 timepoints per individual: takes 5 mins
  saveRDS(res,file=paste0('D:\\files\\R new\\MPCMR\\fPCA_res\\sim_scenarioA_',seed,'.RData') )
  ###B
  set.seed(seed)
  RES<-getX(J=30,ZXmodel='B')
  res <- FPCA(RES$Ly_sim, RES$Lt_sim,list(dataType='Sparse', error=TRUE, verbose=TRUE))
  saveRDS(res,file=paste0('D:\\files\\R new\\MPCMR\\fPCA_res\\sim_scenarioB_',seed,'.RData') )
  ###C
  set.seed(seed)
  RES<-getX(J=30,ZXmodel='E')
  res <- FPCA(RES$Ly_sim, RES$Lt_sim,list(dataType='Sparse', error=TRUE, verbose=TRUE))
  saveRDS(res,file=paste0('D:\\files\\R new\\MPCMR\\fPCA_res\\sim_scenarioC_',seed,'.RData') )
  #res<-readRDS( paste0('D:\\files\\R new\\MPCMR\\fPCA_res\\sim_scenarioC_',seed,'.RData') )
  #write.csv(RES$DAT,paste0('D:\\files\\R new\\MPCMR\\fPCA_res\\sim_scenarioC_DAT_',seed,'.csv'), row.names=T)
}

simFPCA(1)

##or parallel version

n_cores_used<-detectCores()-1
cl<-makeCluster(n_cores_used)

clusterEvalQ(cl=cl , expr=library(data.table))
clusterEvalQ(cl=cl , expr=library(fdapace))

clusterExport(  cl=cl ,  varlist=c('getX', 'simFPCA' )  ) #HPC中记得在R script加上getX的表达式

parSapply(   cl ,  1:100, simFPCA ) #每个simFPCA耗时10*3=30mins
#不需要储存在一个对象里了，因为每个cores都有saveRDS()
#注意：跑的中途不会立刻显示RDS在文件夹中，但是强制结束还是会显示的
stopCluster(cl)

##or better to use HPC





###Paper results content:
#1. MPCMR fitting plot for different Z-X scenarios and fitting method; under the fixed X-Y scenario = '3'
#2. simple results (individual curves, fPCA fitting curve, eigenfunction curves, etc)

#3. MSE&Coverage table (main selling result)

#4. Weak instrument assessment
#5. IV validity assessment


#6. One-parameter fitting (under the constant-effect assumption; i.e. Scenario A)
#7. Constant-effect test (only valid when IV is valid; i.e. no pleiotropy effect)


###result1: MPCMR fitted curve: 9 subplots------------------------------------------------
###---------------------------------------------------------------------------------------
XYmodel_used<-'2'



GGDATA<-c()
for(ZXmodel_used in c('E','B','G')){
  seed=1
  #res<-readRDS( paste0('D:\\files\\R new\\MPCMR\\fPCA_res\\sim_scenario',ZXmodel_used,'_',seed,'.RData') )
  res<-readRDS( paste0('/Users/haodongtian/Documents/MPCMR/FPCAres/sim_scenario',ZXmodel_used,'_',seed,'.RData') )
  set.seed(seed)
  RES<-getX(J=30,ZXmodel=ZXmodel_used)#和上方readRDS中文件名保持一致
  #由于getY也具有随机性，最好前面也跟着set.seed(seed+1000)
  set.seed(seed+1000)
  DAT<-getY(RES, XYmodel=XYmodel_used)#DAT: complete data
  ##naive MPC fitting----------------------------------------------
  MPCMRFiting_res<-naiveMPCMRfit( Gmatrix=RES$DAT[,1:RES$details$J],
                                  Yvector=DAT$Y,
                                  res=res,
                                  XYmodel=XYmodel_used )
  GGdata<-MPCMRFiting_res$ggdata; GGdata$Scenario<-ZXmodel_used ; GGdata$Method<-'Association'
  GGDATA<-rbind( GGDATA , GGdata  )
  #or seperate plot
  # ggsave(paste0('association_',ZXmodel_used,'.eps' ),
  #        plot = MPCMRFiting_res$p ,
  #        path=paste0('C:\\Users\\Haodong\\Desktop\\TVMR\\plots\\'),
  #        height = 4, width = 4.5, units = "in",limitsize=TRUE)


  ##nonparametrc  MPCMR fitting (i.e. basiafunction = eigenfunctions)--------------------------------------
  MPCMRres<-MPCMR_GMM(  Gmatrix=RES$DAT[,1:RES$details$J],
                        res=res,
                        Yvector=DAT$Y,
                        #nL=1,
                        #nLM=20, #with parallel
                        XYmodel= XYmodel_used
                      )
  GGdata<-MPCMRres$ggdata1; GGdata$Scenario<-ZXmodel_used ; GGdata$Method<-'MPCMR'
  GGDATA<-rbind( GGDATA , GGdata  )
  #or seperate plot
  # ggsave(paste0('MPCMR_',ZXmodel_used,'.eps' ),
  #        plot = MPCMRres$p1 ,
  #        path=paste0('C:\\Users\\Haodong\\Desktop\\TVMR\\plots\\'),
  #        height = 4, width = 4.5, units = "in",limitsize=TRUE)

  #parametric fitting (basiafunctions = polynomials)-----------------------------------------------
  GGdata<-MPCMRres$ggdata2; GGdata$Scenario<-ZXmodel_used ; GGdata$Method<-'MPCMR(polynomial)'
  GGDATA<-rbind( GGDATA , GGdata  )
  #or seperate plot
  # ggsave(paste0('MPCMRp_',ZXmodel_used,'.eps' ),
  #        plot = MPCMRres$p2 ,
  #        path=paste0('C:\\Users\\Haodong\\Desktop\\TVMR\\plots\\'),
  #        height = 4, width = 4.5, units = "in",limitsize=TRUE)
}
dim(GGDATA) #51*9=459 9

#single ggplot (注意修改name)
ggdata<-GGDATA
ggdata$Scenario[   ggdata$Scenario=='E'  ]<-'A'  #'E' ->'A'
ggdata$Scenario[   ggdata$Scenario=='B'  ]<-'B'  #'B' ->'B'
ggdata$Scenario[   ggdata$Scenario=='G'  ]<-'C'  #'G' ->'C'


ggdata$Method[ ggdata$Method=='MPCMR'    ]<-'MPCMR(eigenfunction)'

pp<- ggplot(ggdata, aes(time, effect))+
      geom_hline(yintercept = 0,linewidth=0.5,linetype = 2,col='grey' )+
      geom_line(ggdata, mapping =aes(time, true_shape), alpha=1,linewidth=1,col='blue'  )+
      geom_line(ggdata, mapping =aes(time, effect), alpha=1,linewidth=1  )+
      geom_line(ggdata, mapping =aes(time, effect_low), alpha=1,linewidth=1,linetype = 2  )+
      geom_line(ggdata, mapping =aes(time, effect_up), alpha=1,linewidth=1,linetype = 2  )+
      geom_line(ggdata, mapping =aes(time, LM_low), alpha=1,linewidth=1,linetype = 2 , col='#666666' )+
      geom_line(ggdata, mapping =aes(time, LM_up), alpha=1,linewidth=1,linetype = 2  , col='#666666' )+
      labs(x='Age',y='Time-varying effect')+
      facet_grid(rows=vars(Scenario), cols=vars(Method)  )+
      theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())#+coord_cartesian( ylim = c(min(ggdata$effect_low)-0.5*plotdif  ,   max(ggdata$effect_up)+0.5*plotdif ) )
pp
ggsave(paste0('Fig4','.eps' ),
              plot = pp,
              path=paste0('C:\\Users\\Haodong\\Desktop\\TVMR\\plots\\'),
              height = 6, width = 6.5, units = "in",limitsize=TRUE)

ggsave(paste0('Fig4','.eps' ),
       plot = pp,
       path=paste0('/Users/haodongtian/documents/TVMR(latex)/plots/'),
       height = 6, width = 6.5, units = "in",limitsize=TRUE)


###result2: time-varying variable plot-------------------------------------------------------------------
###------------------------------------------------------------------------------------------------------
set.seed(111);RES<-getX(J=30,ZXmodel='A')
indi_plot(RES,123)#图像生成模式是plot(...) 所以不方便储存
set.seed(222);RES<-getX(J=30,ZXmodel='B')
indi_plot(RES,123)
set.seed(333);RES<-getX(J=30,ZXmodel='C')
indi_plot(RES,123)
#ZXmodelA  600 400



#or ggplot version

set.seed(111);RESA<-getX(J=30,ZXmodel='E')
set.seed(222);RESB<-getX(J=30,ZXmodel='B')
set.seed(333);RESE<-getX(J=30,ZXmodel='G')

get_variable_vector<-function(RES, #RES<-getX()
                              i=1 #which individual?
                              ){
  reslist<-list()
  G<-RES$details$G; times<-RES$details$times; X<-RES$details$X; MGX<-RES$details$MGX; UUU<-RES$details$UUU
  reslist$times<-times[-1]
  reslist$exposure<-X[i,]
  reslist$confounding<-UUU[i,]
  reslist$gene_score<-(G%*%t(MGX))[i,]

  return(reslist)
}

ggdata<-c()
for(ii in c('E','B','G')){
  if(ii=='E'){RES<-RESA};if(ii=='B'){RES<-RESB};if(ii=='G'){RES<-RESE}
  rres<-get_variable_vector(RES,1123)
  ggd<-cbind( rres$times,c(rres$exposure,rres$confounding,rres$gene_score), rep(c('Exposure','Confounding','Gene score'),each=1000  ),ii  )
  ggdata<-rbind(  ggdata,ggd )
}

#dim(ggdata) #9000 4=time variable typename ZXmodel

ggdata<-as.data.frame(ggdata)
names(ggdata)<-c('time','variable','variable_type','ZXmodel')
ggdata$time<-as.numeric(  ggdata$time )
ggdata$variable<-as.numeric(  ggdata$variable )

#rename the label
ggdata$ZXmodel[ggdata$ZXmodel=='E']<-'A'
ggdata$ZXmodel[ggdata$ZXmodel=='B']<-'B'
ggdata$ZXmodel[ggdata$ZXmodel=='G']<-'C'


p<- ggplot(ggdata, aes(time, variable))+
  geom_hline(yintercept = 0,linewidth=0.5,linetype = 2,col='grey' )+
  geom_line(ggdata, mapping =aes(time, variable,col=variable_type), alpha=1,linewidth=0.6  )+
  labs(x='Time',y='Variable level')+
  facet_grid(rows=vars(ZXmodel)  )+
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  guides(color=guide_legend(     title='Variable type' )  )
p
ggsave(paste0('variable_plot','.eps' ),
       plot = p,
       path=paste0('C:\\Users\\Haodong\\Desktop\\TVMR\\plots\\'),
       height = 4.5, width = 6, units = "in",limitsize=TRUE)



##R^2 under different scenarios at multiple timepoints
for(ZXmodel_used in  c('A','E','B','G')){
  print(ZXmodel_used )
  RRRresult<-c()
  for(seed_used in  1:1000){ #用1:50就足够稳定了
    set.seed(seed_used)
    RES<-getX(J=30,ZXmodel=ZXmodel_used)
    RRresults<-c()
    for(tt in c(1,10,20,30,40,50)  ){
      RRresults<- c(RRresults ,  summary(lm(  RES$DAT[,30+tt] ~ as.matrix( RES$DAT[,1:30]  )  ))$r.squared )
    }
    RRRresult<-rbind(RRRresult,RRresults )
  }
  RRRresult_mean<- apply( RRRresult, 2, mean )
  print(RRRresult_mean )
}

# "A"
# 0.04055127 0.03152270 0.02501107 0.02104389 0.01821491 0.01598982
# "E"
# 0.04157037 0.03814204 0.03548834 0.02836931 0.02286075 0.02089249
# "B"
# 0.04086078 0.03627661 0.03992328 0.04810993 0.05802020 0.06911997
# "G"
# 0.04132522 0.06023939 0.10852397 0.16422613 0.21832031 0.26899284


###results3: MSE and coverage rate (MSE + General GMM coverage +  robust LM coverage)-----------------------
###and the mean predict value scatterplot
###---------------------------------------------------------------------------------------------------------

MSEvector<-c()
COVvector<-c()
COVvectorLM<-c()
ZXmodelvector<-rep( c('A','B','E') ,  each=6*3 )  #6 XYmodels #3 estimation types
XYmodelvector<-rep( c('0','1','2','3','4','5') , each=3  )  #3 estimation types
Type<-rep(  c('naive','nonpar','par') )

for(ii in c(1:1000)){
  seed<-ii
  for(ZXmodel_used in c('A','B','E')){
    res<-readRDS( paste0('D:\\files\\R new\\MPCMR\\fPCA_res\\sim_scenario',ZXmodel_used,'_',seed,'.RData') )
    set.seed(seed)
    RES<-getX(J=30,ZXmodel=ZXmodel_used) #和上方readRDS中文件名保持一致
    for(XYmodel_used in c('0','1','2','3','4','5') ){
      #由于getY也具有随机性，最好前面也跟着set.seed(seed+1000)
      set.seed(seed+1000)
      DAT<-getY(RES, XYmodel=XYmodel_used)#DAT: complete data
      ##naive MPC fitting----------------------------------------------
      MPCMRFiting_res<-naiveMPCMRfit( Gmatrix=RES$DAT[,1:RES$details$J],
                                      Yvector=DAT$Y,
                                      res=res,
                                      XYmodel=XYmodel_used )
      MSEvector<-c(MSEvector,MPCMRFiting_res$MSE)
      COVvector<-c(COVvector,MPCMRFiting_res$Coverage_rate)
      COVvectorLM<-c(COVvectorLM,NA)  #association目前并不用LM
      ##nonparametrc  MPCMR fitting--------------------------------------
      MPCMRres<-MPCMR_GMM(  Gmatrix=RES$DAT[,1:RES$details$J],
                            res=res,
                            Yvector=DAT$Y,
                            #nL=1,
                            #nLM=20, #with parallel
                            XYmodel= XYmodel_used,
                            )
      MSEvector<-c(MSEvector,MPCMRres$MSE)
      COVvector<-c(COVvector,MPCMRres$Coverage_rate)
      COVvectorLM<-c(COVvectorLM,MPCMRres$Coverage_rate_LM)
      #parametric fitting-----------------------------------------------
      MSEvector<-c(MSEvector,MPCMRres$MSE_p)
      COVvector<-c(COVvector,MPCMRres$Coverage_rate_p)
      COVvectorLM<-c(COVvectorLM,MPCMRres$Coverage_rate_p_LM)
    }
  }
}

resdat<-data.frame(   MSE=MSEvector, COV=COVvector,
                      ZXmodel= ZXmodelvector, XYmodel=XYmodelvector,Type=Type)
dim(resdat)
nrow(resdat)#  2700 = 50*6*4*3  (50 sim times)*(6 XYmodels)*(4 ZXmodels)*(3 estimation types)


#or parallel in HPC
MSE_COV_<-c()
for(ZXmodel_used in c('A','B','E')){
  for(XYmodel_used in c('0','1','2','3') ){
    MSE_COV<-read.csv(paste0( 'C:\\Users\\Haodong\\Desktop\\MPCMRres\\FITRESmean_',ZXmodel_used,  XYmodel_used, '.csv'),header=T,na.strings ="?")
    MSE_COV_<-rbind(  MSE_COV_ ,   as.vector(MSE_COV)$x )
  }
}
MSE_COV_
rownames( MSE_COV_  )<-   c( paste0( 'A', 0:3 ) , paste0( 'B', 0:3 ) ,  paste0( 'E', 0:3 )  )
colnames( MSE_COV_  )<- c('MSE', 'COV' , 'MSE', 'COV' ,  'COV_LM' , 'MSE', 'COV' ,  'COV_LM' )
times100<-function(vector){   vector*100 }

MSE_COV__<-round( apply(  MSE_COV_ , 2, times100) ,3 )
xtable(MSE_COV__,digits = 3)






#使用4个时间点(t=10,20,30,40)上的MSE和coverage rate并使用逆转的表格形式 (in MacOS)
MSE_COV_<-c()
for(ZXmodel_used in c('E','B','G')){ #最后变成 A B C
  for(XYmodel_used in c('0','1','2','3','6','7') ){ #最后变成 1 2 3 4 5 6
    MSE_COV<-read.csv(paste0( '/Users/haodongtian/Documents/MPCMR/MPCMRres/FITRESmeannew_',ZXmodel_used,  XYmodel_used, '.csv'),header=T,na.strings ="?")
    originalvector<-as.vector(MSE_COV)$x  #length(originalvector) ==40
    originalmatrix<-rbind( originalvector[1:6],
                           originalvector[13:18] ,
                           originalvector[25:30]  ,
                           originalvector[7:12],
                           originalvector[19:24] ,
                           originalvector[31:36]     )
    MSE_COV_<-cbind(  MSE_COV_ , c(originalmatrix[,2], originalmatrix[,3], originalmatrix[,4], originalmatrix[,5] )    )
    #as.vector(MSE_COV)$x长度为40 = 36 + cF1 +cF2 +Qp1 +Qp2 where 36= (6 MSE + 6 COV)*3 estimation strategies
  }
}
dim(MSE_COV_)
rownames( MSE_COV_ )<- rep ( c( paste0('MSE', 1:3 ), paste0('COV', 1:3  )     ), length=24)
colnames(  MSE_COV_ )<-c( paste0( 'A', 1:6 ), paste0( 'B', 1:6 ), paste0( 'C', 1:6 )  )


# rownames( MSE_COV_  )<-   c( paste0( rep( c( 'A','B','E','G','D' ) ,each=6) , c(0:3,6,7) ) )
# colnames( MSE_COV_  )<- c(  paste0( rep(c('SE','cov'),each=6  ), 1:6 ),
#                             paste0( rep(c('SE','cov'),each=6  ), 1:6 ),
#                             paste0( rep(c('SE','cov'),each=6  ), 1:6 ),
#                             'cF1','cF2','Qp1','Qp2'
#                             )





View(MSE_COV_)

#change a unit
times100<-function(vector){   vector*100 }
MSE_COV__<-round( apply(  MSE_COV_ , 2, times100) ,3 )


xtable(MSE_COV__,digits = 3) #library(xtable)



#or parallel in HPC (mean predicts)

actual_time<-c( 0.05 ,10.04 ,20.03 ,30.02 ,40.01 ,50.00 )

Predict_<-c()
for(ZXmodel_used in c('E','B','G')){ #最后变成 A B C
  for(XYmodel_used in c('0','1','2','3','6','7') ){ #最后变成 1 2 3 4 5 6
    Predict<-read.csv(paste0( '/Users/haodongtian/Documents/MPCMR/MPCMRres/FITRESmeannew_predict_',ZXmodel_used,  XYmodel_used, '.csv'),header=T,na.strings ="?")
    originalvector<-as.vector(Predict)$x  #length(originalvector) ==18
    matrix_res<-cbind( rep( c( 0,10,20,30,40,50  ), length=24 )    ,
                   c(originalvector, get_true_shape_values(actual_time, XYmodel_used ) ),
                   rep(  c('Strategy 1', 'Strategy 2', 'Strategy 3','True effect'), each=6    ) ,
                   ZXmodel_used,
                   paste0(XYmodel_used,'_'))
    Predict_<-rbind( Predict_ ,  matrix_res )
  }
}

dim( Predict_ ) #432=24*18 5

ggdata<-as.data.frame(Predict_)
colnames(ggdata )<-c(  'Age', 'Estimated_effect','Strategy' , 'ZXmodel' , 'XYmodel'   )
ggdata$Age<-as.numeric(ggdata$Age  )
ggdata$Estimated_effect<-as.numeric(ggdata$Estimated_effect  )
#rename
ggdata$ZXmodel[ggdata$ZXmodel=='E']<-'A'
ggdata$ZXmodel[ggdata$ZXmodel=='G']<-'C'
ggdata$XYmodel[ggdata$XYmodel=='0_']<-'1'
ggdata$XYmodel[ggdata$XYmodel=='1_']<-'2'
ggdata$XYmodel[ggdata$XYmodel=='2_']<-'3'
ggdata$XYmodel[ggdata$XYmodel=='3_']<-'4'
ggdata$XYmodel[ggdata$XYmodel=='6_']<-'5'
ggdata$XYmodel[ggdata$XYmodel=='7_']<-'6'


#rename (as suggested by Steve)
ggdata$Strategy[ggdata$Strategy=='Strategy 1']<-'Association'
ggdata$Strategy[ggdata$Strategy=='Strategy 2']<-'MPCMR(eigenfunction)'
ggdata$Strategy[ggdata$Strategy=='Strategy 3']<-'MPCMR(polynomial)'

ggtrue<-data.frame(  Age=c( 0,10,20,30,40,50  ),
                     effect= as.numeric(  sapply( c(0:3,6,7), function(x){ get_true_shape_values(actual_time,x)  } )  ),
                     XYmodel=rep(  c('1','2','3','4','5','6') , each= 6    )
                       )

p<- ggplot(ggdata, aes(Age, Estimated_effect))+
  geom_hline(yintercept = 0,linewidth=0.5,linetype = 2,col='grey' )+
  geom_point(ggdata, mapping =aes(Age, Estimated_effect,shape=Strategy, col=Strategy) ,size=2.5,alpha=1.0 )+  #position='jitter'
  #geom_point(ggtrue, mapping =aes(Age, effect),size=2  )+
  labs(x='Age',y='Estimated effect')+
  facet_grid(rows=vars(XYmodel),cols=vars(ZXmodel)  )+
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  guides(color=guide_legend(     title=' ' ),shape=guide_legend(     title=' ' )  )
p



ggsave(paste0('predict_scatterplot_rename','.eps' ),
       plot = p,
       path=paste0('/Users/haodongtian/Documents/TVMR(latex)/plots/'),
       height = 9, width = 7, units = "in",limitsize=TRUE)


###result4. Weak instrument assessment-------------------------------------------------------------
###------------------------------------------------------------------------------------------------
cFF<-c()
QQ<-c()
n_sim_used<-1000
for(ZXmodel_used in c('A','B','E')){
  cF<-c()#conditional F values
  Q<-c()#Q p-values
  for(ii in 1:n_sim_used){
    seed<-ii;XYmodel_used<-'1'#XYmodel_used不重要；IS都会是一样的
    #get res
    res<-readRDS( paste0('D:\\files\\R new\\MPCMR\\fPCA_res\\sim_scenario',ZXmodel_used,'_',seed,'.RData') )
    #get RES (seed track)
    set.seed(seed)
    RES<-getX(J=30,ZXmodel=ZXmodel_used) #和上方readRDS中文件名保持一致
    #由于getY也具有随机性，最好前面也跟着set.seed(seed+1000)
    set.seed(seed+1000)
    DAT<-getY(RES, XYmodel=XYmodel_used)#DAT: complete data
    ##nonparametric MPCMR fitting--------------------------------------
    # MPCMRfit_res<-MPCMRfit( Gmatrix=RES$DAT[,1:RES$details$J],   #MPCMRfit uses naive Q statistic
    #                         Yvector=DAT$Y,
    #                         res=res,
    #                         XYmodel=XYmodel_used )
    ##MPCMR_GMM uses advanced Q statistic
    MPCMRfit_res<-MPCMR_GMM(  Gmatrix=RES$DAT[,1:RES$details$J],
                          res=res,
                          Yvector=DAT$Y,
                          LMCI=FALSE,LMCI2=FALSE, #save time
                          nLM=20, #with parallel  #by default =20; bit slower  nLM=10 quite quick
                          Parallel=TRUE,
                          XYmodel= XYmodel_used
    )

    cF<-cbind(cF,   as.vector( MPCMRfit_res$ISres[,3]) )
    Q<-cbind(Q,   as.vector( MPCMRfit_res$ISres[,6]) )
  }
  cFF<-rbind(cFF,cF)
  QQ<-rbind(QQ,round(Q,3))
  print(paste0( 'finished: ZX scenrio' ,ZXmodel_used  )  )
}
#View(cFF);View(QQ)
dim(cFF);dim(QQ) # 6  n_sim_used

#one-value weak instrument strength definition (for better visualization):   \
#minimal cF value  and  maximal Q p-values

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
                     ZXmodel= rep(c('A','B','C') , each=n_sim_used  ) ,
                     IStype=rep(c('Conditional F values' , 'Q p-values'  ) , each=3*n_sim_used  ))
write.csv(ggdata  , paste0('D:\\files\\R new\\MPCMR\\sim_data\\ggdata_weak_IV_results'), row.names=F)

facetlines <- data.frame(IStype = c('Conditional F values', 'Q p-values' ), position = c(10,0.05))
p<- ggplot(ggdata, aes(Values))+
  geom_histogram(aes(y = after_stat(ncount)),bins = 100) +
  facet_grid( rows=vars(ZXmodel), cols=vars(IStype) , scales ='free' )+
  geom_vline(aes(xintercept = position), facetlines,colour='grey', linewidth=0.75 ,linetype = 2)+
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  scale_x_continuous( n.breaks=10)
p

ggsave(paste0('IShist_new.eps' ),   #.eps
       plot = p ,  #非指定
       path=paste0('C:\\Users\\Haodong\\Desktop\\TVMR\\plots\\'),
       height = 5, width = 9, units = "in",limitsize=TRUE)


#not-reject (i.e. believed weak IS) rate under ZXmodel='A' (the known weak instrument case)
sum(as.vector( t(QQ_max) )[1:n_sim_used]>0.05)  #934


###result5: IV validity assessment-----------------------------------------------------------------------
###---------------------------------------------------------------------------------------------
#use ParMPCMRfit where L=nPC to achieve the standard MVMR testing (已证；是一样的)
n_sim_used<-1000

Rarray<-array(NA, dim=c(18,4,n_sim_used))   #18=3 ZXscenarios*6 XYscenarios  #4= 4 pleiotropy size case  #1000 就是fPCA有的数量

Rarray<-array(NA, dim=c(12,6,n_sim_used))

for(ii in c(1:n_sim_used)){
  seed<-ii
  Rmatrix<-c()  #Reject/Result matrix
  for(ZXmodel_used in c('A','B','D')){
    res<-readRDS( paste0('D:\\files\\R new\\MPCMR\\fPCA_res\\sim_scenario',ZXmodel_used,'_',seed,'.RData') )
    set.seed(seed)
    RES<-getX(J=30,ZXmodel=ZXmodel_used) #和上方readRDS中文件名保持一致
    for(XYmodel_used in c('0','1','2','3') ){
      Pvector<-c()#pleiotropic vector
      for(pleiN_used in c(0,1,3,5,10,30)){
        #由于getY也具有随机性，最好前面也跟着set.seed(seed+1000)
        set.seed(seed+1000)
        DAT<-getY(RES, XYmodel=XYmodel_used,plei=TRUE, pleiN=pleiN_used)#DAT: complete data

        ##ParMPCMR fitting--------------------------------------
        MPCMRfit_res<-ParMPCMRfit( Gmatrix=RES$DAT[,1:RES$details$J],
                                   Yvector=DAT$Y,
                                   res=res,
                                   #L=NA,  #默认为nonparametric case;对应的Q test就是IV validity test
                                   XYmodel=XYmodel_used,
                                   Fit=FALSE #没必要fit；直接做Q testing即可
                                   )
        Pvector<-c(Pvector,MPCMRfit_res$pvalue<0.05) #pvalue<0.05即拒绝，即表明存在pleiotropy, 即IV is invalid
      }
      Rmatrix<-rbind(  Rmatrix, Pvector )
    }
  }
  Rarray[,,ii]<-Rmatrix
}

Results<-apply(Rarray, c(1,2), mean   )
xtable( Results, digits=3 ) #library(xtable)




#or parallel in HPC (first download the FITRESmean files)
Result_vector<-c()
for(ZXmodel_used in c('A','B','E')){
  for(XYmodel_used in c('0','1','2','3') ){
    result_vector<-read.csv(paste0( 'C:\\Users\\Haodong\\Desktop\\Results\\Result5_',ZXmodel_used,  XYmodel_used, '.csv'),header=T,na.strings ="?")
    #length(result_vector) #6
    Result_vector<-rbind(  Result_vector ,   result_vector$x )
  }
}

rownames( Result_vector  )<-   c( paste0( 'A', 0:3 ) , paste0( 'B', 0:3 ) ,  paste0( 'E', 0:3 )  )
colnames( Result_vector  )<- c('0', '1' , '3', '5' ,  '10' , '30' )

xtable(Result_vector,digits = 3)







###if weak instrument exists: how we should do?
###result6. One-parameter fitting (under the constant-effect assumption; i.e. Scenario A)
###---------------------------------------------------------------------------------
n_sim_used<-1000


ZXmodel_used<-'A' ###here, with constant G-X effect we can identify the cumulative effect
MSEmatrix<-c();COVmatrix<-c();COVmatrixLM<-c()
CE_MSEmatrix<-c();CE_COVmatrix<-c();CE_COVmatrixLM<-c()
for(ii in c(1:n_sim_used)){
  seed<-ii
  res<-readRDS( paste0('D:\\files\\R new\\MPCMR\\fPCA_res\\sim_scenario',ZXmodel_used,'_',seed,'.RData') )
  set.seed(seed)
  RES<-getX(J=30,ZXmodel=ZXmodel_used) #和上方readRDS中文件名保持一致
  MSEvector<-c(); COVvector<-c();COVvectorLM<-c()
  CE_MSEvector<-c(); CE_COVvector<-c();CE_COVvectorLM<-c()
  for(XYmodel_used in c('0','1','2','3') ){#6种 XYmodels
    #由于getY也具有随机性，最好前面也跟着set.seed(seed+1000)
    set.seed(seed+1000)
    DAT<-getY(RES, XYmodel=XYmodel_used)#DAT: complete data
    ##One-parametric MPCMR fitting----------------------------------------------
    # ParMPCMRfit_res<-ParMPCMRfit( Gmatrix=RES$DAT[,1:RES$details$J],
    #                               Yvector=DAT$Y,
    #                               res=res,
    #                               L=1,  #人为选择L=1 due tO the weak instrument #L=1可以视为正确的time-varying form- 因为是等效的
    #                               XYmodel=XYmodel_used,
    #                               Fit=TRUE )

    ##One -parametric MPCMR GMM fitting-----------------------------------------
    MPCMRres<-MPCMR_GMM(  Gmatrix=RES$DAT[,1:RES$details$J],
                          res=res,
                          Yvector=DAT$Y,
                          #Gymatrix=NA,
                          #IDmatch=NA,
                          #nPC=NA,
                          nL=1,  #人为选择L=1 due tO the weak instrument #L=1可以视为正确的time-varying form- 因为是等效的
                          LMCI=FALSE,LMCI2=TRUE, #nL=1的Polynomial parametric fitting还需要用LM coverage
                          nLM=20, #with parallel  #by default =20; bit slower  nLM=10 quite quick
                          Parallel=TRUE,
                          XYmodel= XYmodel_used
                        )

    ##global (curve) MSE and coverage rate--------------------------------------
    MSEvector<-c(MSEvector ,  MPCMRres$MSE_p)  #这些MSE 是average over Gird timepoint; 这些都要再mean的
    COVvector<-c(COVvector ,  MPCMRres$Coverage_rate_p)
    COVvectorLM<-c(COVvectorLM ,  MPCMRres$Coverage_rate_p_LM)

    ##CE (cumulative effect) MSE and coverage rate-----------------------------------------

    #普通的Gaussian inference
    est<-as.numeric(MPCMRres$MPCMRest_p); se<-sqrt( as.numeric(MPCMRres$MPCMRvar_p)  )
    CEest<-50*est   ; CEse<-50*se
    if(XYmodel_used=='0'){CEtrue = 0  }
    if(XYmodel_used=='1'){CEtrue = 5  }
    if(XYmodel_used=='2'){CEtrue = 25  }
    if(XYmodel_used=='3'){CEtrue = 0  }
    if(XYmodel_used=='4'){CEtrue = 2  }
    if(XYmodel_used=='5'){CEtrue = 2  }
    CE_MSEvector<-c( CE_MSEvector ,  (CEest-CEtrue)^2 )
    CE_COVvector<-c( CE_COVvector ,  abs( CEest - CEtrue  )<1.96*CEse  )

    #LM的condidence set
    CE_COVvectorLM<-c( CE_COVvectorLM , ( ( CEtrue/50 -   MPCMRres$ggdata2$LM_low[1] )*( CEtrue/50 - MPCMRres$ggdata2$LM_up[1] ) )<0 )

  }
  MSEmatrix<-cbind(  MSEmatrix,MSEvector  )
  COVmatrix<-cbind(  COVmatrix,COVvector  )
  COVmatrixLM<-cbind(  COVmatrixLM,COVvectorLM  )
  CE_MSEmatrix<-cbind(  CE_MSEmatrix,CE_MSEvector  )
  CE_COVmatrix<-cbind(  CE_COVmatrix,CE_COVvector  )
  CE_COVmatrixLM<-cbind(  CE_COVmatrixLM,CE_COVvectorLM  )

}
apply(  MSEmatrix, 1, mean )
apply(  COVmatrix, 1, mean )
apply(  CE_MSEmatrix, 1, mean )
apply(  CE_COVmatrix, 1, mean )

xtable(  data.frame(
  Scenario=paste0(  'A',0:3 )   , CE=c( 0,5,25,0,2,2) ,
  MSE=apply(  MSEmatrix, 1, mean )*50^2 ,Coverage=apply(  COVmatrix, 1, mean ),CoverageLM=apply(  COVmatrixLM, 1, mean ),
  MSE2=apply(  CE_MSEmatrix, 1, mean ) ,Coverage2=apply(  CE_COVmatrix, 1, mean ), Coverage2LM=apply(  CE_COVmatrixLM, 1, mean )   ),
  digits=3  )
#注意curve MSE需要* [time lenght]^2才能位置和 CE MSE一个数量比较级






###result7: Constant-effect test (GMM IV validity test use nL=1)---------------------
###----------------------------------------------------------------------------------
#array style
n_sim_used<-1000
Parray<-array(NA, dim=c(15,3,n_sim_used))  #1000 就是模拟的总数量


for(ii in c(1:n_sim_used)){
  seed<-ii
  Pmatrix<-c()  #final dim(Pmatrix) #15 3
  for(ZXmodel_used in c('A','B','E')){#3 Z-X scenarios
    Pvector<-c()
    res<-readRDS( paste0('D:\\files\\R new\\MPCMR\\fPCA_res\\sim_scenario',ZXmodel_used,'_',seed,'.RData') )
    set.seed(seed)
    RES<-getX(J=30,ZXmodel=ZXmodel_used) #和上方readRDS中文件名保持一致
    for(XYmodel_used in c('111','222','333') ){
      if(XYmodel_used=='111' ){ bcandidates<-c( 0,0.1, 0.5, 1, 10 ) }
      if(XYmodel_used=='222' ){ bcandidates<-c(0.001,0.01,0.05, 0.1 ,1 ) }
      if(XYmodel_used=='333' ){ bcandidates<-c(0.01,0.1,0.5,1, 10 ) }
      for(b_used in  bcandidates ){
        #由于getY也具有随机性，最好前面也跟着set.seed(seed+1000)
        set.seed(seed+1000)
        DAT<-getY(RES, XYmodel=XYmodel_used,b=b_used)#DAT: complete data
        ##parametric  MPCMR fitting with L=1--------------------------------------
        # MPCMRres<-MPCMR_GMM(  Gmatrix=RES$DAT[,1:RES$details$J],
        #                       res=res,
        #                       Yvector=DAT$Y,
        #                       nL=1,
        #                       LMCI=FALSE #save time
        #                      )
        # Pvector<-c(Pvector,MPCMRres$IV_validity_and_basisfunction_test[3])

        ###alternative: IVW q test
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

dim(Parray) # 15 3 n_sim_used

Reject_array<-Parray<0.05

apply(Reject_array, c(1,2), mean   )

xtable( apply(Reject_array, c(1,2), mean   ),digit=3    ) #library(xtable )

#test
seed<-146
ZXmodel_used<-'C'; XYmodel_used<-'111' ;b_used<-10
res<-readRDS( paste0('D:\\files\\R new\\MPCMR\\fPCA_res\\sim_scenario',ZXmodel_used,'_',seed,'.RData') )
set.seed(seed)
RES<-getX(J=30,ZXmodel=ZXmodel_used) #和上方readRDS中文件名保持一致
DAT<-getY(RES, XYmodel=XYmodel_used,b=b_used)

#one-sample GMM
MPCMRres<-MPCMR_GMM(  Gmatrix=RES$DAT[,1:RES$details$J],
                      res=res,
                      Yvector=DAT$Y,
                      nL=1,
                      LMCI=FALSE #save time
)

MPCMRres$p1
MPCMRres$p2
MPCMRres$IV_validity_test
MPCMRres$IV_validity_and_basisfunction_test

#IVW Q
MPCMRFiting_res<-ParMPCMRfit( Gmatrix=RES$DAT[,1:RES$details$J],
                              Yvector=DAT$Y,
                              res=res,
                              L=1,  #L=1 means time-varying effect is constant
                              Fit=FALSE #fit 与否不重要
)

c(MPCMRFiting_res$Qvalue,  MPCMRFiting_res$Qdf,   MPCMRFiting_res$pvalue)

#conclusion:
#already check via simulation that the transformed exposure (via B) will lea to the same numric Q result as the previous IVW Q

#though IVW Q and GMM Q has the same asymptotic form/distribution
#IVW Q can better control the type I error


#or parallel in HPC (first download the FITRESmean files)
Result_vector<-c()
for(XYmodel_used in c('111','222','333')){
  if(XYmodel_used=='111' ){ bcandidates<-c( 0.00, 0.05, 0.1, 0.5, 1 ) }
  if(XYmodel_used=='222' ){ bcandidates<-c(0.001,0.005,0.01,0.05, 0.1  ) }
  if(XYmodel_used=='333' ){ bcandidates<-c(0.01,0.05,0.1, 0.5,1 )  }
  for(b_used in  bcandidates ){
    result_vector<-read.csv(paste0( 'C:\\Users\\Haodong\\Desktop\\Results\\Result7mean_', XYmodel_used,b_used ,'.csv'),header=T,na.strings ="?")
    #length(result_vector) #3
    Result_vector<-rbind(  Result_vector ,   c( b_used,  result_vector$x) )
  }
}

colnames( Result_vector  )<- c('betasize', 'A' ,  'B' , 'E' )

xtable(Result_vector,digits = 3)
