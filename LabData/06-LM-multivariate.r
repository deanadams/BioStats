# Multivariate GLM: MANOVA, Multivariate Regression, MANCOVA
         #Packages: RRPP, mvabund

library(RRPP) 
library(vegan)
library(mvabund)

bumpus <- read.csv("Data/bumpus.csv",header=T)
bumpus.data <- log(as.matrix(bumpus[,(5:13)])) # matrix of linear measurements
sex <- as.factor(bumpus[,2])
surv <- as.factor(bumpus[,4])
SexBySurv <- paste(sex,surv)
TotalLength <- bumpus.data[,1]
Y <- bumpus.data[,-1]
mydat <- rrpp.data.frame(Y=Y,sex=sex,surv=surv,TotalLength=TotalLength,SexBySurv = SexBySurv)

#__________________________________________________________________________#
#Describe your data
(100*apply(bumpus.data,2,sd)) /apply(bumpus.data,2,mean)  # CV: Coefficient of Variation

cor(bumpus.data)
pairs(bumpus.data)
var(bumpus.data)
var(scale(bumpus.data))
as.matrix(dist(bumpus.data))

#__________________________________________________________________________#
#single factor MANOVA
model1 <- lm(bumpus.data~sex)
summary(model1)	#yields a set of univariate analyses

summary(manova(model1))	#does multivariate test (using Pillai's)
summary(manova(model1),test="Wilks")	#does multivariate test (using Wilks)

##### MANOVA via RRPP
model.rrpp <- lm.rrpp(Y~sex,data = mydat, print.progress = FALSE)
anova(model.rrpp)
plot(model.rrpp, type = "PC", pch=21, bg = sex)  #PC PLOT!
legend("topright", levels(sex), pch = 21, pt.bg = 1:4)

##When parametric methods break down
Ynew <- matrix(rnorm(100000), nrow=100) #100 x 100 matrix: N=p
gp <- gl(2,50)
summary(manova(lm(Ynew~gp)))  #parametric algebra cannot be completed
anova(lm.rrpp(Ynew~gp,print.progress = FALSE))  #NO problem with RRPP!

#__________________________________________________________________________#
#Factorial MANOVA
model2<-lm(bumpus.data~sex*surv)
summary(manova(model2))

#Factorial MANOVA via RRPP
model2.rrpp <- lm.rrpp(Y~sex*surv,data = mydat, print.progress = FALSE)
anova(model2.rrpp)
groups <- interaction(mydat$sex, mydat$surv)
plot(model2.rrpp, type = "PC", pch=21, bg = groups)
legend("topright", levels(groups), pch = 21, pt.bg = 1:4)

#Also in vegan (FRPP only)
adonis(formula = Y~sex*surv, method = "euclidean" )

#Pairwise comparisons 
PW <- pairwise(model2.rrpp,groups = SexBySurv, print.progress = FALSE)
summary(PW,test.type = "dist")  #No pairwise groups significant AFTER accounting for main effects!

model.null <- lm.rrpp(Y~1,data=mydat, print.progress = FALSE)
summary(pairwise(model2.rrpp, model.null, groups = SexBySurv, 
                 print.progress = FALSE),test.type = "dist")
   #Male vs. female seem to be biggest differences

#__________________________________________________________________________#
### Multivariate Regression
summary(manova(lm(Y~TotalLength)))

model.reg <- lm.rrpp(Y~TotalLength, data = mydat, print.progress = FALSE)
anova(model.reg)

### Visualizing multivariate regression 
plot(model.reg, type = "regression", reg.type = "RegScore", 
     predictor = mydat$TotalLength, pch=19)

#__________________________________________________________________________#
#MANCOVA	
summary(manova(lm(Y~TotalLength*sex*surv)))
summary(manova(lm(Y~TotalLength+sex*surv))) # FIT COMMON SLOPE

#MANCOVA via RRPP
model.mancova <- lm.rrpp(Y~TotalLength*sex*surv, data =mydat, print.progress = FALSE)
anova(model.mancova)

#Pairwise comparisons of slopes: NOTE: Used here for illustrative purposes only! 
   # Interaction should show evidence justifying this analysis
PW.mancova <- pairwise(model.mancova, groups = SexBySurv, covariate = TotalLength, print.progress = FALSE)
summary(PW.mancova, test.type = "VC", angle.type = "deg")

### Visualizing MANCOVA
plot(model.mancova, type = "regression", reg.type = "RegScore", 
     predictor = mydat$TotalLength, pch=19, col = as.numeric(groups))

plot(model.mancova, type = "regression", reg.type = "PredLine", 
     predictor = mydat$TotalLength, pch=19,
     col = as.numeric(groups))

#________________________________________________________________________#
# A more complex example
data(Pupfish)
Pupfish$Group <- interaction(Pupfish$Sex, Pupfish$Pop)

#MANOVA via RRPP
fit.m <- lm.rrpp(coords ~ Pop * Sex, data = Pupfish, print.progress = FALSE) 
anova(fit.m)$table
summary(pairwise(fit.m,groups = Pupfish$Group),test.type = "dist")

plot(fit.m, type = "PC", pch=21, bg = Pupfish$Group, cex=2)
legend("topright", levels(Pupfish$Group), pch = 21, pt.bg = 1:4)

#MANCOVA via RRPP
Pupfish$logSize <- log(Pupfish$CS)
fit.slopes <- lm.rrpp(coords ~ logSize * Pop * Sex, data = Pupfish, print.progress = FALSE)
anova(fit.slopes)$table
PWS <- pairwise(fit.slopes, groups = Pupfish$Group, covariate = Pupfish$logSize)
summary(PWS, test.type = "VC", angle.type = "deg")

plot(fit.slopes, type = "regression", reg.type = "RegScore", pch=21, bg = Pupfish$Group, predictor = Pupfish$logSize, cex=2)
legend("topleft", levels(Pupfish$Group), pch = 21, pt.bg = 1:4)

plot(fit.slopes, type = "regression", reg.type = "PredLine", pch=21, bg = Pupfish$Group, predictor = Pupfish$logSize, cex=2)
legend("topright", levels(Pupfish$Group), pch = 21, pt.bg = 1:4)

#________________________________________________________________________#
#MANCOVA with Y as Distance matrix 
mydat$Ydist <- dist(Y)
mancova.dist <- lm.rrpp (Ydist~TotalLength*sex*surv, data =mydat, print.progress = FALSE)
anova(mancova.dist)

#________________________________________________________________________#
#Methods for species abundance matrices:
  #1: MANOVA via RRPP: performed on log(abundance)
  #2: MANOVA via RRPP: performed on Bray-Curtis distances
  #3: mvAbund: generalized linear models on each species separately, then summed

data("spider")  #species abundance data across localities
dim(spider$abund)

# transformation of raw data (log+1)
spid <- data.frame(spider$x)
spid$abund <- as.matrix(spider$abund)
spid$logY <- as.matrix(log(spider$abund + 1)) # log-transform

### 1: MANOVA via RRPP on log(Y) 
fitLM <- lm.rrpp(logY ~ soil.dry + bare.sand + fallen.leaves + moss
                 + herb.layer + reflection, 
                 data = spid, iter = 999, SS.type = "II",
                 print.progress = FALSE)
anova(fitLM)
pca <- plot(fitLM, type = "PC")
text(pca$PC.points, rownames(pca$PC.points), pos = 1, cex = 0.5)

### 2: MANOVA via RRPP on Bray-Curtis distance (ala Anderson)
Db <- vegdist(spid$abund, method="bray")

fitD <- lm.rrpp(Db ~ soil.dry + bare.sand + fallen.leaves + moss
                + herb.layer + reflection, 
                data = spid, iter = 999, SS.type = "II",
                print.progress = FALSE)
anova(fitD)  #pretty similar

### 3: Generalized linear models: mvAbund (ala Wharton)
spiddat <- mvabund(spider$abund)
X <- spider$x

#To fit a log-linear model assuming counts are poisson:
  #NOTE: mvabund fits generlized LM to each Y-variable separately, then sums F-values
glm.spid <- manyglm(spiddat~X, family="poisson")
summary(glm.spid)   #not the same as above
anova(glm.spid)
