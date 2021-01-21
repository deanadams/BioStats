# Various ANOVA models
   #Packages:   RRPP, lme4, lattice

bumpus<-read.csv("Data/bumpus.csv",header=T)
bumpus.data<-log(as.matrix(bumpus[,(5:13)])) # matrix of log-linear measurements
sex<-as.factor(bumpus[,2])
surv<-as.factor(bumpus[,4])
SexBySurv<-as.factor(paste(sex,surv))
 TotalLength<-bumpus.data[,1]

 ## Multi-group histogram
 sexdata<-split(TotalLength,sex)
 hist(sexdata$m,col="blue",freq=TRUE)
 hist(sexdata$f,col="red",freq=TRUE,add=TRUE)

  #ANOTHER WAY, SHOWING OVERLAP
 hist(sexdata$m, col=rgb(0,0,0,1),freq=TRUE)
 hist(sexdata$f, col=rgb(0.8,0.8,0.8,0.5), freq=TRUE,add=T)

#__________________________________________________________________________#
#Single-Factor ANOVA
 #  Some components of the linear model
 model1<-lm(TotalLength~sex)
 anova(model1)  	#Statistical summary
 coef(model1)   #for ANOVA, these are means gp1 and dev for mean gp 2
 predict(model1)		#Predicted values
 resid(model1)		#Residuals
 
 anova(lm(TotalLength~sex))
plot(sex,TotalLength,ylab="Total Length")
  summary(aov(TotalLength~sex))    #another way to do univariate ANOVA

##### ANOVA via residual randomization
library(RRPP) 
mydat <- rrpp.data.frame(TotalLength = TotalLength, sex = sex, surv=surv, SexBySurv = SexBySurv)

res <- lm.rrpp(TotalLength~sex, print.progress = FALSE, data = mydat)
anova(res)     #note: all summary and parameter values, as well as permutation components, are saved: see res$ANOVA$...
coef(res)

### Single-Factor ANOVA (4 groups)
anova(lm(TotalLength~SexBySurv))

#via RRPP
model2 <- lm.rrpp(TotalLength~SexBySurv, print.progress = FALSE, data = mydat)
anova(model2)

#post-hoc tests
pairwise.t.test(TotalLength, SexBySurv, p.adj = "none")  #standard

#pairwise comparisons via RRPP
posthoc <- pairwise(fit = model2,groups = SexBySurv, print.progress = FALSE)
summary(posthoc)

#__________________________________________________________________________#
###Factorial ANOVA
anova(lm(TotalLength~sex+surv+sex:surv))

###SHORTHAND for FULL Factorial
anova(lm(TotalLength~sex*surv))

#RRPP
model3 <- lm.rrpp(TotalLength~sex*surv, print.progress = FALSE, data = mydat)
anova(model3)

##Full randomization (FRPP) vs. residual randomization (RRPP)
anova(lm.rrpp(TotalLength~sex*surv, print.progress = FALSE, RRPP = TRUE, data = mydat))
anova(lm.rrpp(TotalLength~sex*surv, print.progress = FALSE,RRPP = FALSE, data = mydat))

model.r <- lm.rrpp(TotalLength~sex+surv, print.progress = FALSE, data = mydat)
model.n <- lm.rrpp(TotalLength~1, print.progress = FALSE, data = mydat)

###NOTE: pairwise comparison interpretations can change when using different reduced models
summary(pairwise(model3, groups = SexBySurv, print.progress = FALSE))  #  Default: against H_r: Y=sex+surv
summary(pairwise(model3, fit.null = model.r,groups = SexBySurv, print.progress = FALSE))  #against H_r: Y=sex+surv
summary(pairwise(model3, fit.null = model.n,groups = SexBySurv, print.progress = FALSE))  #against H_r: Y~1

#__________________________________________________________________________#
###Nested ANOVA 
anova(lm(TotalLength~sex/surv)) #survival nested within sex  
  #NOTE: F-tests not correct: all are based on MSE; but sex should be tested against nested effect!

#Using RRPP
model.nested <- lm.rrpp(TotalLength~sex/surv, effect.type = "F", print.progress = FALSE, data = mydat)
anova(model.nested)  #unadjusted
anova(model.nested, error = c("sex:surv", "Residuals"))  #  Correct MSE specification
   #inspect F for sex relative to MS-sex and MS-sex:surv

#__________________________________________________________________________#
## Mixed Models
library(lme4)
library(lattice)

#plot the data: see the random slopes/intercepts issue
xyplot(Reaction ~ Days | Subject, sleepstudy, type = c("g","p","r"),
       index = function(x,y) coef(lm(y ~ x))[1],
       xlab = "Days of sleep deprivation",
       ylab = "Average reaction time (ms)", aspect = "xy")

# Random slopes & intercepts model
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
anova(fm1)
summary(fm1)

#Random intercepts model
fm2 <- lmer(Reaction ~ Days + (1|Subject), data = sleepstudy) 
summary(fm2)
anova(fm2)
anova(lm(Reaction~Days*Subject, data= sleepstudy))
anova(lm(Reaction~Days, data= sleepstudy))

#__________________________________________________________________________#
## Type I vs II vs III SS:  

### easy to do in RRPP (harder in 'base' R)
anova(lm.rrpp(TotalLength~sex*surv, SS.type="I", print.progress = FALSE, data = mydat))  #sequential SS
anova(lm.rrpp(TotalLength~sex*surv, SS.type="II", print.progress = FALSE, data = mydat)) #conditional SS
anova(lm.rrpp(TotalLength~sex*surv, SS.type="III", print.progress = FALSE, data = mydat)) #marginal SS

#__________________________________________________________________________#
## SIMPLE RESAMPLING EXAMPLE
F.obs<-anova(lm(TotalLength~sex))[1,4]  #Find Test value and save
permute<-999
F.rand.vec<-array(NA,(permute+1))
F.rand.vec[permute+1]<-F.obs
for(i in 1:permute){
  ###Shuffle Data
	y.rand<-sample(TotalLength)  #SHUFFLE FIRST COL
	F.rand.vec[i]<-anova(lm(y.rand~sex))[1,4]  
}  
P.Ftest<-rank(F.rand.vec[permute+1])/(permute+1) 
F.obs
P.Ftest

hist(F.rand.vec,40,freq=T,col="gray")
segments(F.obs, 0, F.obs, 50)  ##Plot Observed value

