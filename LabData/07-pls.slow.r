pls.slow<-function(x,y,iter=999){
  data.all<-cbind(x,y)
  data.all<-scale(data.all,scale=FALSE)
  XY.vcv<-cov(data.all)
  S12<-XY.vcv[1:dim(x)[2],(dim(x)[2]+1):(dim(x)[2]+dim(y)[2])]
  pls<-svd(S12)
  U<-pls$u; V<-pls$v
  XScores<-x%*%U[,1]; YScores<-y%*%V[,1]
  PLS.obs<-cor(XScores,YScores)
  P.val<-1
  PLS.rand<-NULL
  for(i in 1:iter){
    y.r<-y[sample(nrow(y)),]	
    data.r<-cbind(x,y.r)
    data.r<-scale(data.r,scale=FALSE)
    XY.vcv.r<-cov(data.r)
    S12.r<-XY.vcv.r[1:dim(x)[2],(dim(x)[2]+1):(dim(x)[2]+dim(y.r)[2])]
    pls.r<-svd(S12.r)
    U.r<-pls.r$u; V.r<-pls.r$v
    XScores.r<-x%*%U.r[,1]; YScores.r<-y.r%*%V.r[,1]
    corr.rand<-cor(XScores.r,YScores.r)
    PLS.rand<-rbind(PLS.rand,corr.rand)
    P.val<-ifelse(corr.rand>=PLS.obs, P.val+1,P.val) 
  }  
  P.val<-P.val/(iter+1)
  return(list(PLS.corr=PLS.obs,pvalue=P.val))
}
