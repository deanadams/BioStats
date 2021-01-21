# Various regression models
    #packages: RRPP, lmodel2, scatterplot3d

bumpus<-read.csv("Data/bumpus.csv",header=T)
bumpus.data<-log(as.matrix(bumpus[,(5:13)])) # matrix of log-linear measurements
sex<-as.factor(bumpus[,2])
surv<-as.factor(bumpus[,4])
SexBySurv<-as.factor(paste(sex,surv))
Y<-as.matrix(bumpus.data[,1])
X1<-bumpus.data[,2]
X2<-bumpus.data[,3]

library(RRPP)
mydat <- rrpp.data.frame(Y = Y, X1 = X1, X2 = X2, sex = sex, surv = surv, SexBySurv = SexBySurv)

#__________________________________________________________________________#
#simple linear regression  (Model I Regression)
lm(Y~X1)

#more information
model1<-lm(Y~X1)
summary(model1)	#provides regression coefficients	
anova(model1)	#provides model term tests
plot(X1,Y,pch=21, bg="black", cex=2)
abline(model1,lwd=2,col="red")

  #Lots of components in model, such as:
model1$fitted.values
model1$residuals

#Additional Diagnostic plots
plot(model1)
	
#plot 1: fitted vs. residuals (structure in plot BAD)
	#plot 2: residual quantiles (near line GOOD)
	#plot 3: fitted vs. sqrt resid (in a triangle BAD)
	#plot 4: cooks distance (identify points with high leverage)

#Regression evaluated via residual randomization (RRPP)
model2 <- lm.rrpp(Y~X1, print.progress = FALSE, data = mydat)
anova(model2)
anova(model1)	#Identical to parametric results 
coef(model2)
coef(model1)
#__________________________________________________________________________#
#### Model II Regression 
library(lmodel2)
lmodel2(Y~X1,nperm=999)
  RMA<-lmodel2(Y~X1)
plot(RMA, pch=21,cex=2, bg="black")
abline(model1,lwd=2,col="blue")

#__________________________________________________________________________#
#multiple regression
summary(lm(Y~X1+X2))
anova(lm(Y~X1+X2))

#via RRPP
anova(lm.rrpp(Y~X1+X2,print.progress = FALSE, data=mydat))
  cor(X1,X2)  #hmm, there is multicollinearity in the X-variables. Perhaps use type II SS.
anova(lm.rrpp(Y~X1+X2,print.progress = FALSE, data=mydat,SS.type = "II"))

# Plot for multiple regression
#install.packages("scatterplot3d")
library(scatterplot3d)
plot<-scatterplot3d(X1,X2,Y)
plot$plane3d(lm(Y~X1+X2))

#__________________________________________________________________________#
#polynomial regression
fit  <- lm(Y~X1) #first degree
fit2 <- lm(Y~poly(X1,2,raw=TRUE))#second degree
fit3 <- lm(Y~poly(X1,3,raw=TRUE))#third degree
fit4 <- lm(Y~poly(X1,4,raw=TRUE))#fourth degree

#evaluate models
anova(fit)
anova(fit,fit2)  #In this case, adding polynomial terms NOT an improvement
anova(fit2,fit3)
anova(fit3,fit4)

plot(X1,Y)
abline(model1,col="red")
xx <-seq(min(X1),max(X1),length=100)
lines(xx, predict(fit2, data.frame(X1=xx)), col="blue")
lines(xx, predict(fit3, data.frame(X1=xx)), col="green")
lines(xx, predict(fit4, data.frame(X1=xx)), col="purple")

#__________________________________________________________________________#
###ANCOVA
anova(lm(Y~X2*SexBySurv))

#Implies slopes for M and F statistically equivalent, so drop and re-run
anova(lm(Y~X2+SexBySurv))
model.ancova<-lm(Y~X2+SexBySurv)

colo<-rep("blue",nrow(bumpus.data)); colo[which(SexBySurv == 'f TRUE')]<-"red";
	colo[which(SexBySurv == 'f FALSE')]<-"pink"; colo[which(SexBySurv == 'm FALSE')]<-"lightblue";

plot(X2,Y,pch=21,cex=2,bg=colo)
abline(model.ancova$coefficients[1],model.ancova$coefficients[2])
abline((model.ancova$coefficients[1]+model.ancova$coefficients[3]),model.ancova$coefficients[2])
abline((model.ancova$coefficients[1]+model.ancova$coefficients[4]),model.ancova$coefficients[2])
abline((model.ancova$coefficients[1]+model.ancova$coefficients[5]),model.ancova$coefficients[2])
	#note: parameters are intercept, slope, deviations of int. for groups

#via RRPP
model.anc2 <- lm.rrpp(Y~X2*SexBySurv, print.progress = FALSE, data = mydat)
anova(model.anc2)

model.anc3 <- lm.rrpp(Y~X2+SexBySurv, print.progress = FALSE, data = mydat)
#NOTE: `plot` works on lm.rrpp objects, but is most useful for  #multivariate data

plot(X2,Y,pch=21,cex=2,bg=colo)
abline(model.ancova$coefficients[1],model.ancova$coefficients[2])
abline((model.ancova$coefficients[1]+model.ancova$coefficients[3]),model.ancova$coefficients[2])
abline((model.ancova$coefficients[1]+model.ancova$coefficients[4]),model.ancova$coefficients[2])
abline((model.ancova$coefficients[1]+model.ancova$coefficients[5]),model.ancova$coefficients[2])

### Pairwise comparisons: test of means
pairwise.t.test(model.ancova$fitted.values, SexBySurv, p.adj = "none")  #standard approach

#Pairwise comparisons of means via RRPP
summary(pairwise(model.anc3, groups = SexBySurv),test.type = "dist") 

#Pairwise comparisons of SLOPES (NOTE: makes sense only when interaction term significant)
PW <- pairwise(model.anc3, groups = SexBySurv)
summary(PW, test.type = "VC", angle.type = "deg")  #comparison of angles between vectors


#__________________________________________________________________________#
## SIMPLE RESAMPLING EXAMPLE
F.obs<-anova(lm(Y~X1))$F[[1]]  #Find Test value and save
permute<-1999
F.rand.vec<-array(NA,(permute+1))
F.rand.vec[permute+1]<-F.obs
for(i in 1:permute){
  ###Shuffle Data
	y.rand<-sample(Y)	#Resample vector 
	F.rand.vec[i]<-anova(lm(y.rand~X1))$F[[1]]  
}  
F.obs
P.Ftest<-rank(F.rand.vec[permute+1])/(permute+1)
P.Ftest
####Plot
hist(F.rand.vec,40,freq=T,col="gray")
segments(F.obs, 0, F.obs, 50)  ##Plot Observed value
