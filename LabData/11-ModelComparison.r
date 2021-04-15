# Model Selection AIC Examples 

# Read data from Smith and Collyer 2008
snake<-read.csv("data/Lab-11-snake.csv",header=T)
site<-as.factor(snake[,1]);region<-as.factor(snake[,2]);sex<-as.factor(snake[,3])
svl<-snake[,4]; hs<-snake[,5] # hs= head size, calculated as centroid size from lm
Y<-as.matrix(snake[,-(1:5)])

##Setup: Univariate and multivariate models
    #Univariate
hs.svl<-lm(hs~svl, model=T,x=T)
hs.sex<-lm(hs~sex, model=T,x=T)
hs.reg<-lm(hs~region, model=T,x=T)
hs.svl.reg<-lm(hs~svl + region, model=T,x=T)
hs.svl.by.reg<-lm(hs~svl*region, model=T,x=T)
hs.svl.sex<-lm(hs~svl + sex, model=T,x=T)
hs.svl.by.sex<-lm(hs~svl*sex, model=T,x=T)
hs.ancova<-lm(hs~svl + sex*region, model=T,x=T)
hs.full<-lm(hs~svl*sex*region, model=T,x=T)
    #multivariate 
Y.svl<-lm(Y~svl, model=T,x=T)
Y.sex<-lm(Y~sex, model=T,x=T)
Y.reg<-lm(Y~region, model=T,x=T)
Y.svl.reg<-lm(Y~svl+region, model=T,x=T)
Y.svl.by.reg<-lm(Y~svl*region, model=T,x=T)
Y.svl.sex<-lm(Y~svl + sex, model=T,x=T)
Y.svl.by.sex<-lm(Y~svl*sex, model=T,x=T)
Y.mancova<-lm(Y~svl + sex*region, model=T,x=T)
Y.full<-lm(Y~svl*sex*region, model=T,x=T)

########### Model Comparisons

#1: Nested Model comparisons:  LRT
anova(hs.full)  #with type I (sequential SS), each term addition is an LRT
anova(hs.svl.reg)
anova(hs.svl,hs.svl.reg)

#Multivariate is the same
anova(Y.full)
anova(Y.svl.reg,Y.svl)
anova(Y.svl.reg,Y.reg)
anova(Y.mancova,Y.svl.reg)  #NOTE: sex not an improvement of model

#2: AIC comparisons
aic.summary<-AIC(hs.svl,hs.sex,hs.reg,hs.svl.reg,hs.svl.by.reg,hs.svl.sex,hs.svl.by.sex,
hs.ancova,hs.full)

aic.summary  #smallest AIC is preferred model

#Model Averaging
c1<-exp(-.5*0)
c2<-exp(-.5*(AIC(hs.svl.by.reg)-AIC(hs.svl.reg)))
w1<-c1/(c1+c2)  #AIC weights
w2<-c2/(c1+c2)
  w1
  w2
beta1<-w1*(hs.svl.reg$coef)
beta2<-w2*(hs.svl.by.reg$coef)

beta.avg<-c((beta1[1]+beta2[1]),(beta1[2]+beta2[2]),(beta1[3]+beta2[3]),
            (beta1[4]+beta2[4]),beta2[5],beta2[6])
beta.avg #model averaged coefficients

#FIND AIC for averaged model (need likelihood and its parameters)
Terms<-terms(formula(hs.svl.by.reg))  #set up X matrix
X<-model.matrix(Terms)
resid<-hs-X%*%beta.avg
E<-t(resid)%*%resid  #Note: univariate, so this IS the det(E)

n<-nrow(resid); p<-1  #NOTE: p<-ncol(resid) for multivariate
k<-hs.svl.by.reg$rank # same as df in "AIC" except "AIC" adds one for variance
LLik<- (-2*((-.5*n*(log(E/n^p)+p))-(0.5*(n*log(2*pi)))))
pen<-2*(p*k+0.5*p*(p+1))
AIC.mdl.avg<-LLik+pen

AIC(hs.svl.by.reg)
AIC(hs.svl.reg)
AIC.mdl.avg  #NOTE: AIC of average model is WORSE than the 2 input models

#NOTE, these are nested models, so standard LRT is used
anova(hs.svl.reg)
anova(hs.svl.by.reg)
anova(hs.svl.reg,hs.svl.by.reg) # Adding interaction NOT an improvement

#2b: Multivariate model comparison
AIC(Y.svl)  #Need multivariate equivalent (and with correct parameter penalty)
  # Multivariate AIC function (from M. Collyer)
maic<-function(x,...){ # each x is a model
    y<-resid(x)
    E<-t(y)%*%y
    d<-det(E)
    if(length(E)>1) n<-nrow(y) else n<-length(y)
    if(length(E)>1) p<-ncol(y) else p<-1
    k<-x$rank # same as df in "AIC" except "AIC" adds one for variance
    lh<-n*(log(d/n^p)+p)
    LLik<- (-2*((-.5*n*(log(d/n^p)+p))-(0.5*(n*log(2*pi)))))
    pen<-2*(p*k+0.5*p*(p+1))
    m<-LLik+pen
    maic<-c(k,m)
}  

maic.svl<-maic(Y.svl)
maic.sex<-maic(Y.sex)
maic.reg<-maic(Y.reg)
maic.svl.reg<-maic(Y.svl.reg)
maic.svl.by.reg<-maic(Y.svl.by.reg)
maic.svl.sex<-maic(Y.svl.sex)
maic.svl.by.sex<-maic(Y.svl.by.sex)
maic.mancova<-maic(Y.mancova)
maic.full<-maic(Y.full)

rbind(maic.svl,maic.sex,maic.reg,maic.svl.reg,maic.svl.by.reg,maic.svl.sex,
	maic.svl.by.sex,maic.mancova,maic.full)

# Note, maic simplifies to univariate
aic.svl<-maic(hs.svl)
aic.reg<-maic(hs.reg)
aic.svl.reg<-maic(hs.svl.reg)

aic.svl
AIC(hs.svl)

# compare maic to AIC
aic.reg[2]-aic.svl.reg[2] # using maic
AIC(hs.reg)-AIC(hs.svl.reg) # canned R AIC function

##3: Cross validation
x<-seq(1:10)
y<-c(2,2,3,4,4,7,7,7,8,8)
plot(x,y)
anova(lm(y~x))

#subsample, fit model, residuals for cross-validated data, 
  #variance across iterations
iter=1000
rep<-lapply(1:iter, function(j) sample(seq(1:10),5))
diff<-lapply(1:iter, function(j) setdiff(seq(1:10),rep[[j]]))
  #full model
model.x<- lapply(1:iter, function(j) lm(y[rep[[j]]]~x[rep[[j]]])) 
resid.y<-lapply(1:iter, function(j)  resid(model.x[[j]],newdata=data.frame(c(x[diff[[j]]]))))
ss.x<-unlist(lapply(1:iter, function(j) crossprod(resid.y[[j]])))
  #reduced model
model.1<- lapply(1:iter, function(j) lm(y[rep[[j]]]~1)) 
resid.1<-lapply(1:iter, function(j)  resid(model.1[[j]],newdata=data.frame(c(x[diff[[j]]]))))
ss.1<-unlist(lapply(1:iter, function(j) crossprod(resid.1[[j]])))

c(mean(ss.1),var(ss.1))
c(mean(ss.x),var(ss.x))
