

###pre-real script: clean the data for formal real analysis (data stored in the folder Data in Git)

library(data.table)#as.data.table()
library(fdapace)#PACE
library(tidyverse) #一些data.table %>% 操作
library(corrplot)#correlated SNPs


###The 4 dataset:(各自的quality control和特殊处理)

###Dat1: ID reference and covariate dataset--------------------------------------------------
master<-read.csv('/Users/haodongtian/Documents/files/R new/MPCMR/master.csv')
dtp<-as.data.table(master)

###no any initial selection (to increase the sample size)
#dtp<-dtp[dtp$euro==1,]
#dtp<-dtp[dtp$race=='White',]
#dtp<-dtp[dtp$sex=='Female',]

Dat1<-dtp

dim(Dat1) #488366    204
write.csv(Dat1  , paste0('/Users/haodongtian/Documents/files/R new/MPCMR/sim_data_new2/Dat1.csv'), row.names=F)
#Dat1<-read.csv(paste0('/Users/haodongtian/Documents/files/R new/MPCMR/sim_data_new2/Dat1.csv'))

###Dat2: genetic information dataset---------------------------------------------------------
hannah_bp_euro<-read.csv('/Users/haodongtian/Documents/files/R new/MPCMR/hannah_bp_euro.csv')
dtg<-as.data.table(hannah_bp_euro)
#remove missing data + round SNP values + highly-correlated SNPs check
dtg<-na.omit(dtg)
dtg<-cbind(dtg[,1], round(    dtg[,   -1 ]      ))
table(  as.vector(as.matrix(dtg[,-1]) ))    #check SNP values 0,1,2  #checked
sort( as.numeric(cor( dtg[,   -1 ] ) ) , decreasing = TRUE  )[    (ncol(dtg[,   -1 ]) - 2): (ncol(dtg[,   -1 ]) - 2+10)    ]
#前3项必是1.000000 #第4项很关键：就是max cor     #max cor=0.3    so r^2~= 0.09 < 0.2 the common threshold
#其实全打入即使SNPs相关也不要紧；只要保证design matrix不要太奇异了即可
#但是如果使用Q statistics; 最好使用independent SNPs for better joint (profile) likelihood maximization
Dat2<-dtg

dim(Dat2) #367643    259= ID + 258 SNPs
write.csv(Dat2  , paste0('/Users/haodongtian/Documents/files/R new/MPCMR/sim_data_new2/Dat2.csv'), row.names=F)
#Dat2<-read.csv(paste0('/Users/haodongtian/Documents/files/R new/MPCMR/sim_data_new2/Dat2.csv'))

###Dat3: longitudinal exposure dataset---------------------------------------------------------
load("/Users/haodongtian/Documents/files/R new/MPCMR/ukb_gp_qrisk_SB.rdata")
Dat3_sbp<-ukb_gp_qrisk[ ukb_gp_qrisk$exposure_type =='sbp', ] #only focus on sbp as the exposure

#coarsened + exposure mean; final returned as a string
getString<-function( exposures,
                     times ,
                     coarsend_level=1   #level 使用1好了（而不是0），即保留一个小数点
){
  ctimes<-round(  times , coarsend_level )  #coarsened times
  dat<-data.frame( X= exposures, time= ctimes     )
  dat_mean<- dat %>% group_by( time) %>%  summarise(mean(X)) #group下的分类求mean/sort等操作
  single_string<-paste(c(dat_mean$`mean(X)` ,   dat_mean$time),collapse='_')
  return(single_string  )
}

Dat3<- Dat3_sbp %>% group_by( idno) %>%  summarise(getString(exposure_value,exposure_age)) #group下的分类求mean/sort等操作

dim(Dat3)#158310      2

#checked: Dat3没问题
write.csv(Dat3  , paste0('/Users/haodongtian/Documents/files/R new/MPCMR/sim_data_new2/Dat3.csv'), row.names=F)
#Dat3<-read.csv(paste0('/Users/haodongtian/Documents/files/R new/MPCMR/sim_data_new2/Dat3.csv'))

###Dat4: outcome dataset-------------------------------------------------------------------
#关于标号所对应的信息，参考https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=30670
#修改后面的id=*信息就好, 至于30670.0.0  30670.1.0 这种，其实是Instance的标号,往往用0就行，因为这个instance样本数量最多
all_biomarker<-read.csv('/Users/haodongtian/Documents/files/R new/MPCMR/all_biomarker_variables.csv')

##creatinine = all_biomarker$"30700-0.0"   #measurement date = all_biomarker$"30701-0.0"
#outcome_dat<-cbind( all_biomarker$eid, all_biomarker$X30700.0.0, all_biomarker$X30701.0.0 )

##urea = all_biomarker$"30670-0.0"  #measurement date = all_biomarker$"30671-0.0"
outcome_dat<-cbind( all_biomarker$eid, all_biomarker$X30670.0.0, all_biomarker$X30671.0.0 )

#match the dob (date of birth) values according to master reference dataset
dpt_index<-match( outcome_dat[,1], master$Adiposity_sample_ID    )
outcome_dat<-cbind(outcome_dat, master$dob[dpt_index]  )

outcome_dat<-na.omit( outcome_dat)

get_ages<-function( dom , dob   ){ #dom: date of measurement ;   dob: date of birth
  DOM<-matrix( as.numeric(unlist(strsplit( dom , '-'   )  )) , length( dom )   , 3,byrow = TRUE ) #3= 年 - 月 -日
  DOB<-matrix( as.numeric(unlist(strsplit( dob , '-'   )  )) , length( dob )   , 3,byrow = TRUE )
  age<- DOM%*%c( 1, 1/12 , 1/365  ) - DOB%*%c( 1, 1/12 , 1/365  )
  return(age)
}
outcome_dat<-cbind( outcome_dat , get_ages(outcome_dat[,3],outcome_dat[,4])      )

Dat4<-data.frame( eid=outcome_dat[,1] , Y=outcome_dat[,2] , age=outcome_dat[,5] )

dim(Dat4)  #465031      3 #ID + urea + timepoint

write.csv(Dat4  , paste0('/Users/haodongtian/Documents/files/R new/MPCMR/sim_data_new2/Dat4.csv'), row.names=F)
#Dat4<-read.csv(paste0('/Users/haodongtian/Documents/files/R new/MPCMR/sim_data_new2/Dat4.csv'))














###data combination




###------------------------------------------------------------------------------------------
###One-sample individual data: (instruments exposure) dataset, combining Dat1+Dat2+Dat3->Dat123 ----------------------------

#review:     match(c(3,1,5,6) , c(4,5,3,9) ) # 3 NA  2 NA
#Dat1 交 Dat2 -> Dat12------------------------------------------------
ID<-match(  Dat1$UKB_sample_ID, Dat2$X )  #按照前者的顺序

Dat2_<-Dat2[ID,]
Dat12<-cbind(  Dat1$Adiposity_sample_ID ,  Dat2_[,-1]  , Dat1[,-(1:4)]  )  #ID + 258 SNPs + nuisance phenotype variables
names(Dat12)[1]<-'V1'
ncol( Dat2_[,-1] )  # 258 SNPs

###Dat12存储一下， 主要是为了后面的naive MR analysis (naive是指不考虑time-varying)
Dat_n<-cbind(Dat12[,1:259],Dat12$sbp)  #dim(Dat_n) #197128    260
#removes (ID,Gs,Xt) NA items
ID<-  apply(  is.na(Dat_n[,1:260]  )  ,  1, sum  )==0
Dat_n<-Dat_n[ID,]
dim(Dat_n) #367283        260  #ID + 258 SNPs + SBP(measured at study baseline)
names(Dat_n)[260]<-'V2'
write.csv(Dat_n  , paste0('/Users/haodongtian/Documents/files/R new/MPCMR/sim_data_new2/Dat_n.csv'), row.names=F)
#Dat_n<-read.csv(paste0('/Users/haodongtian/Documents/files/R new/MPCMR/sim_data_new2/Dat_n.csv'))



#Dat12 交 Dat3 -> Dat123 ----------------------------------------------
ID<-match(  Dat12$V1, as.numeric(Dat3$idno) )

Dat3_<-Dat3[ID,]
Dat123<-cbind(  Dat12[,1:259]  ,  Dat3_[,2]  ,   Dat12[,-c(1:259)])

#Dat123 remove missing data -> storage as Dat---------------------------
dim(Dat123)  #488366    460
#names(Dat123)[260]
names(Dat123)[260]<-'Xt'

#removes (ID,Gs,Xt) NA items
ID<-  apply(  is.na(Dat123[,1:260]  )  ,  1, sum  )==0
Dat<-Dat123[ID,]


dim(Dat)  #116378    460

write.csv(Dat  , paste0('/Users/haodongtian/Documents/files/R new/MPCMR/sim_data_new2/Dat.csv'), row.names=F)
#Dat<-read.csv(paste0('/Users/haodongtian/Documents/files/R new/MPCMR/sim_data_new2/Dat.csv'))

nrow( na.omit(Dat[,1:260]))==nrow(Dat)#TRUE




###------------------------------------------------------------------------------------------
###One-sample outcome dataset, combining Dat1+Dat2+Dat4->Dat124 ----------------------------
#Dat1 交 Dat2 -> Dat12------------------------------------------------
ID<-match(  Dat1$UKB_sample_ID, Dat2$X )  #按照前者的顺序

Dat2_<-Dat2[ID,]
Dat12<-cbind(  Dat1$Adiposity_sample_ID ,  Dat2_[,-1]  , Dat1[,-(1:4)]  )
names(Dat12)[1]<-'V1'
ncol( Dat2_[,-1] )  # 258 SNPs

#Dat12 交 Dat4 -> Dat124 ----------------------------------------------
ID<-match(  Dat12$V1, as.numeric(Dat4$eid) )

Dat4_<-Dat4[ID,]
Dat124<-cbind(  Dat12[,1:259]  ,  Dat4_[,2:3]  ,   Dat12[,-c(1:259)])

#Dat124 remove missing data -> storage as DatY---------------------------
dim(Dat124)  #488366    461 #ID + 258SNPs + urea + timepoint + nuisance phenotypes
#names(Dat124)

#removes (ID,Gs,Y, age) NA items
ID<-  apply(  is.na(Dat124[,1:261]  )  ,  1, sum  )==0
DatY<-Dat124[ID,]
dim(DatY)  #350223    461

DatY$Y<-as.numeric(DatY$Y);DatY$age<-as.numeric(DatY$age)
write.csv(DatY  , paste0('/Users/haodongtian/Documents/files/R new/MPCMR/sim_data_new2/DatY.csv'), row.names=F)
#DatY<-read.csv(paste0('/Users/haodongtian/Documents/files/R new/MPCMR/sim_data_new2/DatY.csv'))

nrow( na.omit(DatY[,1:261]))==nrow(DatY)#TRUE




