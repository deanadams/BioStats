## a quick comparison of computational time: how implementation matters
   #Packages microbenchmark

library(microbenchmark)

x<-matrix(rnorm(10000),ncol=2)
xf<-cbind(1,x)
y<-matrix(rnorm(nrow(x)))

lm(y~x)  #Common method
solve(t(xf)%*%xf)%*%t(xf)%*%y
crossprod(solve(crossprod(xf)),crossprod(xf,y))
lm.fit(xf,y)$coefficients
.lm.fit(xf,y)$coefficients  ### NOTE: a very low-level function (cannot use in packages submitted to CRAN)
qr.coef(qr(xf),y)

microbenchmark(
  lm(y~x),
  solve(t(xf)%*%xf)%*%t(xf)%*%y,
  crossprod(solve(crossprod(xf)),crossprod(xf,y)),
  lm.fit(xf,y),.lm.fit(xf,y),
  qr.coef(qr(xf),y)
)

###NOTE that the best implementation can change with the size of the data matrix
#Large X univ. Y
x<-matrix(rnorm(10000),ncol=50)
xf<-cbind(1,x)
y<-matrix(rnorm(nrow(x)))
microbenchmark(
  lm(y~x),
  solve(t(xf)%*%xf)%*%t(xf)%*%y,
  crossprod(solve(crossprod(xf)),crossprod(xf,y)),
  lm.fit(xf,y),.lm.fit(xf,y),
  qr.coef(qr(xf),y)
)

##Large Y univ. X
y<-matrix(rnorm(10000),ncol=100)
x<-matrix(rnorm(nrow(y)))
xf<-cbind(1,x)
microbenchmark(
  lm(y~x),
  solve(t(xf)%*%xf)%*%t(xf)%*%y,
  crossprod(solve(crossprod(xf)),crossprod(xf,y)),
  lm.fit(xf,y),.lm.fit(xf,y),
  qr.coef(qr(xf),y)
)

#large Y and X
y<-matrix(rnorm(20000),ncol=100)
x<-matrix(rnorm(10000),ncol=50)
xf<-cbind(1,x)
microbenchmark(
  lm(y~x),
  solve(t(xf)%*%xf)%*%t(xf)%*%y,
  crossprod(solve(crossprod(xf)),crossprod(xf,y)),
  lm.fit(xf,y),.lm.fit(xf,y),
  qr.coef(qr(xf),y)
)

