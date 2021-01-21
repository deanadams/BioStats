# Performing Matrix Computations in R & Basic GLM
       #packages: MASS 

bumpus<-read.csv("Data/bumpus.csv",header=T)
bumpus.data<-as.matrix(bumpus[,(5:13)])  #note: 'as.matrix' forces data type
  cor(bumpus.data)    #correlation  matrix
  vcv.bumpus<-var(bumpus.data)    #covariance matrix
  var(bumpus.data)

#__________________________________________________________________________#
# Elementary matrix algebra

a<-matrix(c(1,0,4,2,-1,1),nrow=3)
b<-matrix(c(1,-1,2,1,1,0),nrow=2)
a
b

c<-t(a)	#matrix transpose
a
c

d<-matrix(c(2,3,1,4),nrow=2)
det(d)	#matrix determinant
exp(determinant(d)$modulus)  #another way

2*a	#scalar multiplication

#Matrix addition and subtraction
b+c
b-c
a+b		##NOTE: non-conformable matrices (check rxc of your matrices!)

a#elementwise multiplication (Hadamard product)
c
b
c*b

diag(5)   #identity matrix with 5 rows/cols

# matrix multiplication
a%*%b		## %*% is symbol for matrix multiplication
b%*%a		## matrix order matters

#__________________________________________________________________________#
# matrix inversion
d
solve(d)
solve(vcv.bumpus)

### example of redundant variables causing singularities
a <- rnorm(10)
b <-rnorm(10)
c <- a+b
cov.dat <- cov(cbind(a,b,c))
solve(cov.dat)


# Generalized Inverse
library(MASS)
ginv(vcv.bumpus)		#useful for singular covariance matrices (when P>N, or reduncancies in data)

#__________________________________________________________________________#
### Matrix decomposition
#eigen-analysis
eigen(vcv.bumpus)			#decomposition of square-symmetric matrices

#singular-value decomposition
svd(vcv.bumpus)			#Same results as 'eigen'
svd(bumpus.data)			#Can also decompose rectangular matrices

#__________________________________________________________________________#
#GLM in matrix form
rm(list=ls())
x<-matrix(1,10)
x2<-matrix(c(1,1,1,1,1,0,0,0,0,0))
xnew<-cbind(x,x2)
y<-matrix(c(5,4,4,4,3,7,5,6,6,6))
yreg<-matrix(c(1,3,4,6,8,7,9,11,10,13))
xreg<-matrix(c(1,2,3,4,5,6,7,8,9,10))
xnewreg<-cbind(x,xreg)

#__________________________________________________________________________#
#ANOVA
model1<-lm(y~x2)
  model.matrix(model1)
b<-solve(t(xnew)%*%xnew)%*%t(xnew)%*%y
summary(model1)
b

predict(model1)
yhat<-xnew%*%b
yhat

SSFull<-t(y-yhat)%*%(y-yhat)
bred<-solve(t(x)%*%x)%*%t(x)%*%y
yhatred<-x%*%bred
yhatred
SSRed<-t(y-yhatred)%*%(y-yhatred) #Is SSTot
SSModel<-SSRed-SSFull
anova(model1)
SSModel
SSFull
SSRed

#__________________________________________________________________________#
#Regression
model2<-lm(yreg~xreg)
breg<-solve(t(xnewreg)%*%xnewreg)%*%t(xnewreg)%*%yreg
summary(model2)
breg

predict(model2)
yhatreg<-xnewreg%*%breg
yhatreg

SSFullreg<-t(yreg-yhatreg)%*%(yreg-yhatreg)
bredreg<-solve(t(x)%*%x)%*%t(x)%*%yreg
yhatredreg<-x%*%bredreg
yhatredreg
SSRedreg<-t(yreg-yhatredreg)%*%(yreg-yhatredreg) #Is SSTot
SSModelreg<-SSRedreg-SSFullreg
anova(model2)
SSModelreg
SSFullreg
SSRedreg

