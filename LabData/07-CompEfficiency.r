#Thoughts on speeding up (and cleaning up) R code

############ let us begin with a comparison of implementations
mymean<-function(x){
  n<-length(x)
  tmp<-0
  for (i in 1:n){
    tmp<-tmp+x[i]
  }
  mn<-tmp/n
  return(mn)
}

x<-matrix(rnorm(1000))
n<-length(x)

library(microbenchmark)
library(ggplot2)
microbenchmark(mean(x),apply(x,2,mean), sum(x)/length(x),mymean(x),colSums(x)/length(x))
  #Note that the function 'mean' is NOT the fastest: sum(x)/length(x) is

############ Evaluating choke-points in code
#Example
library(aprof)
source("07-pls.slow.r")
tmp<-tempfile() #create tmp file for saving profiler output
Rprof(tmp,line.profiling=TRUE)  #profile the function
x<-matrix(rnorm(1000),ncol=10)
y<-matrix(rnorm(1000),ncol=10)
pls.slow(x,y)
Rprof(append=FALSE)
fooaprof<-aprof("07-pls.slow.r",tmp) #Create aprof object
plot(fooaprof)

###### Some General suggestions for speeding up code
#1: pre-allocate memory  

iter<-99
SS<-array(NA,iter) #pre-allocate
newSS<-function(iter){ #'on the fly'
  SS<-NULL
  for (i in 1:iter){SS<-rbind(SS,NA)}
  return(SS)
}

microbenchmark(SS<-array(NA,99), x<-newSS(99))
microbenchmark(SS<-array(NA,9999), x<-newSS(9999),times=10)


#2: pre-calculate objects 

x<-cbind(1,matrix(rnorm(1000),ncol=10))
y<-matrix(rnorm(100))

all.calc<-function(x,y){
  coef.r<-array(NA,dim=c(999,ncol(x)))
  for (i in 1:999){
    y.r<-y[sample(nrow(y)),]	
    coef.r[i,]<-solve(t(x)%*%x)%*%t(x)%*%y.r
  }
}

hat.calc<-function(x,y){
  hat<-solve(t(x)%*%x)%*%t(x)
  coef.r<-array(NA,dim=c(999,ncol(x)))
  for (i in 1:999){
    y.r<-y[sample(nrow(y)),]	
    coef.r[i,]<-hat%*%y.r
  }  
}

microbenchmark(all.calc(x,y),hat.calc(x,y),times=10)


#3:  Use lower-level functions 
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

#4: Vectorize when possible. Don't speak R with a 'C accent'
fn1<-function(x){
  means<-array(0,ncol(x))
  for(i in 1:ncol(x)){
    for(j in 1:nrow(x)){
      means[i]<-means[i]+x[j,i]
    }
  }
  means<-means/nrow(x)
  return(means)
}

x<-matrix(rnorm(1000*1000),ncol=1000)
microbenchmark(fn1(x),colMeans(x),apply(x,2,mean),times=10)


###
x <- matrix(rnorm(1000*10000), ncol=1000)
fn1<-function(x){
  mx <- rep(NA, nrow(x))
  for(i in 1:nrow(x)){ mx[i] <- max(x[i,])  } 
  return(mx)
}

microbenchmark(fn1(x),apply(x,1,max),times=10)  #loop is faster here


#####PLS: compare new and old versions
source('07-pls.fast.r')
source('07-pls.slow.r')
x<-matrix(rnorm(10000),ncol=10)
y<-matrix(rnorm(20000),ncol=20)
microbenchmark(pls.slow(x,y),pls.fast(x,y),times=5)

