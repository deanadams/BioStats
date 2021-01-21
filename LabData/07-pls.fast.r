#Speed it up using recommendations

pls.fast<-function(x,y,iter=999){
  n<-nrow(x)
  x<-scale(x,scale=FALSE); y<-scale(y,scale=FALSE)
  #permutation. 1st is observed
  ind <- c(list(1:n), (Map(function(x) sample.int(n, n), 1:iter)))
  y.rand <- lapply(1:(iter+1), function(i) y[ind[[i]], ])
  S12.r<-lapply(1:(iter+1), function(i) crossprod(x, y.rand[[i]])/(n - 1))
  pls.r<-lapply(1:(iter+1), function(i) La.svd(S12.r[[i]], 1, 1))
  r.rand<-sapply(1:(iter+1), function(i) cor(x %*% pls.r[[i]]$u, y.rand[[i]] %*% t(pls.r[[i]]$vt)))
  p.val<- 1-(rank(r.rand)[1])/(iter+1)
  return(list(PLS.corr=r.rand[[1]],pvalue=p.val))
}

