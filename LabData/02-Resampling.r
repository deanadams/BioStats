# Various resampling procedures
#      NOTE: Some routines are not the most efficient means of accomplishing the
     #   Task at hand. Advanced R-programming may be used to enhance computational performance.

        #Packages: mvtnorm, boot

#__________________________________________________________________________#
#Generating Random Numbers
	#Of paramount importance when running Monte Carlo simulations 
rm(list=ls())
#normal distribution
x <- rnorm(5000)
hist(x,20,freq=T)
x2<-runif(5000)
hist(x2,20,freq=T)

z1<-rnorm(500, mean = 1, sd = 0.3)
hist(z1,20,freq=T)

z2<-rnorm(500, mean = 3, sd = 0.3)
plot(z1,z2)

z3<-runif(500, 0, 1)
z4<-runif(500, 0,1)
plot(z3,z4)

##Setting the 'seed' for random number generation
	#IMPORTANT for checking simulations prior to full run
set.seed(2)  #sets the seed
rnorm(10)
set.seed(2)  #sets the seed
rnorm(10)
rnorm(10)		#note that the seed is not held indefinitely in memory)

## Simulate CORRELATED Data Using MVTNORM
	##basic idea: simulate rnorm vectors, & multiply by decomposition of a covariance matrix
		## covariance matrix specifies correlation structure
library(mvtnorm) 
corr.val<-0.7
a <- rmvnorm(n=500,mean=c(0,0),sigma=matrix(c(1,corr.val,corr.val,1),2,2)) 
cor(a)
plot(a)

#__________________________________________________________________________#
###Shuffle Data: i.e. resampling without replacement
rm(list=ls())
x <- 1:10
#dim(x)<-c(10,1) #force dimensions to make column vector
x
sample(x)	# Randomize the order of locations
sample(x,replace=FALSE)	#more explicit
sample(x,(length(x)-2),replace=FALSE)	#sub-sample

#__________________________________________________________________________#
# WHEN SAMPLE CAN GET US IN TROUBLE!!!!!!
x <- 1:10
y<-cbind(x,x,x)
y

sample(y,replace=FALSE)    #We did not tell it to preserve rows!!
y[sample(nrow(y)),]   	   # ROWS Preserved (needed for resampling multivariate data)

#__________________________________________________________________________#
#Simple Randomization Test
rm(list=ls())
bumpus<-read.csv("Data/bumpus.csv",header=T)
bumpus.data<-log(bumpus[,(5:13)]) # set of log-linear measurements
sex<-as.factor(bumpus[,2])
TL<-bumpus.data$TL

#Observed data
t.test(formula=bumpus.data$TL~sex)	
plot(sex,bumpus.data$TL,ylab="Total Length")
t.obs<-t.test(formula=bumpus.data$TL~sex)$statistic[[1]]  #grab t-statistic from frame	
t.obs

#Randomization Test
permute<-999
P.1tailed<-P.2tailed<-1
t.rand.vec<-array(NA,(permute+1))
t.rand.vec[permute+1]<-t.obs
for(i in 1:permute){
  ###Shuffle Data
	y.rand<-sample(bumpus.data$TL)  #NOTE: notation works ONLY for single variable 
  ###Run analysis on random data
	t.rand.vec[i]<-t.test(formula=y.rand~sex)$statistic[[1]]
}  #end permute

##Significance assessment
P.1tailed<-length(which(t.rand.vec<=t.rand.vec[permute+1])) / (permute+1)  #because observed is negative
P.2tailed<-length(which(abs(t.rand.vec)>=abs(t.rand.vec[permute+1]))) / (permute+1)  
P.1tailed
P.2tailed
####Plot
hist(t.rand.vec,20,freq=T,col="gray")
segments(t.obs, 0, t.obs, 50)  ##Plot Observed value

#####NOTE: Distribution centered on 0.0.  This is the expected value under the null hypothesis:
    #thus, we have correctly shuffled the valid exchangeable units under Ho.



#__________________________________________________________________________#
# Bootstrapping Data 
x <- 1:10
sample(x,replace=TRUE)

#Observed mean and CI
mean.TL<-mean(TL)	
mean.TL
int<-1.96*sqrt(var(TL)/length(TL))
CIlow<-mean.TL-int
CIhi<-mean.TL+int
CIlow
CIhi
#Bootstrap data
boot.mean<-numeric(1000)
for (i in 1:1000){
  boot.mean[i]<-mean(sample(TL,replace=T))
}

#__________________________________________________________________________#
### PERCENTILE BOOTSTRAP CI
quantile(boot.mean,c(0.025,0.975))

#__________________________________________________________________________#
###Standard Boostrap CI
int.boot<-1.96*sd(boot.mean)
CIlow.boot<-mean.TL-int.boot
CIhi.boot<-mean.TL+int.boot
CIlow.boot
CIhi.boot

#__________________________________________________________________________#
###### ANOTHER approach with library 'boot' [From Crawley: The R Book]
#install.packages("boot")
library(boot)
mymean<-function(TL,i)mean(TL[i])
myboot<-boot(TL,mymean,R=1000)
myboot

# Percentile Bootstrap CI
quantile(myboot$t,c(0.025,0.975))
#Various Bootstrap CI
boot.ci(myboot)

rm(list=ls())

