###plotEigenfunction


plotEifun<-function(res #res<-FPCA(...)
                    ){
  K<-ncol(res$phi)
  lambdas<-round(res$lambda,3)
  VE<-lambdas/sum(lambdas)*100  #fraction variance explained
  VE<-round(VE,1  )
  ggdata<-data.frame( time=res$workGrid ,
                      values=as.vector(res$phi) ,
                      PC=rep(   paste0( 'PC' , 1:K, '(',VE  ,'%)'  ),each=nrow(res$phi) )) #PC: principal components

  p<-ggplot(ggdata, aes(time, values,PC))+
      geom_hline(yintercept = 0,linewidth=1,linetype = 2,col='grey' )+
      geom_line(ggdata, mapping =aes(x=time, y=values,color=PC), alpha=1,linewidth=1  )+
      labs(x='Age',y='Eigenfunction')+
      theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())#+coord_cartesian( ylim = c(min(ggdata$values)-0.5  ,   max(ggdata$values)+0.5 ) )
  return(p)
}
