## GLS Simulations comparing FRPP, RRPP, RRPP-alt, rrpp-restricted
  #### NOTE: find way to simulated BM with uniform, log-normal, and exp^3
     #modify geomorph fast.sim.BM ??

library(RRPP)
library(phytools)
library(caper)
source('rrpp-alt.R')
source('rrpp-rest.r')

nsims <- 100
iter <- 999
n <- 100
groups <- 2

result <- t(sapply(1:nsims, function(j){
  b1 <- 0.0
  b2 <- 0.0  
  tree <- pbtree (n = n)
  Cov <- vcv.phylo(tree)
  e <- fastBM(tree)
  x <- rep(seq(0,(groups-1)), n/groups)
  x2 <- rep(seq(0,(groups-1)), each = n/groups)
  names(x) <- names(x2) <- names(e)
  y <- (x * b1) + (x2 * b2) + (x * b1)*(x2 * b2) + e 
  x <- factor(x); x2 <- factor(x2)
  df <- data.frame(x = x, x2 = x2, y = y, species = names(x))
  df2 <- comparative.data(data = df, tree, 'species')
  
  fit.param <- pgls(y~x*x2, data = df2)
  fit.RRPP <- lm.rrpp(y ~ x*x2, data = df, Cov = Cov,
                    print.progress = FALSE, RRPP = TRUE, iter= iter)
  fit.FRPP <- lm.rrpp(y ~ x*x2, data = df, Cov = Cov, 
                       print.progress = FALSE, RRPP = FALSE, iter= iter)
  fit.RRPP.alt <- rrpp.alt(y ~ x*x2, data = df, Cov = Cov, 
                           print.progress = FALSE, RRPP = TRUE, iter= iter)
  fit.rest <- rrpp.rest(y ~ x+x2, data = df, SS.type = "II", Cov = Cov, 
                        print.progress = FALSE, RRPP = TRUE, iter= iter)

  aov.param <- anova(fit.param)
  aov.RRPP <- anova(fit.RRPP)
  aov.FRPP <- anova(fit.FRPP)

  cat("iter", j, " ")
  c(b.par = coef(fit.param)[-1],
    b.RRPP = fit.RRPP$LM$gls.coefficients[-1],
    b.FRPP = fit.FRPP$LM$gls.coefficients[-1],
    F.par = aov.param[1:3,4], 
    F.RRPP = aov.RRPP$table$F[1:3], 
    F.FRPP = aov.FRPP$table$F[1:3],
    F.Rest =fit.rest[,1],
    F.RRPP.tBr = fit.RRPP.alt[,1],
    p.par = aov.param[1:3,5], 
    p.RRPP = aov.RRPP$table$`Pr(>F)`[1:3], 
    p.FRPP = aov.FRPP$table$`Pr(>F)`[1:3],
    p.Rest =  apply(fit.rest,1,RRPP:::pval),
    p.RRPP.tBr =  apply(fit.RRPP.alt,1,RRPP:::pval)
  )
}))

#save(result, file = "Results/Results.GLSFact.Rdata")
ResultB <- data.frame(result[,c(1:9)])
ResultF <- result[,c(10:23)]
ResultP <- result[,c(24:37)]

#####

png(file="Figs/GLS.Fact-f1.png")
pairs(ResultF[,c(1,4,7,10,12)], pch=19, 
      labels = c("Parametric","RRPP", "FRPP","Restricted","RRPP-Alt"), main = "F-X1")
dev.off()

png(file="Figs/GLS.Fact-f2.png")
pairs(ResultF[,c(2,5,8,11,13)], pch=19, 
      labels = c("Parametric","RRPP", "FRPP","Restricted","RRPP-Alt"), main = "F-X2")
dev.off()

png(file="Figs/GLS.Fact-f3.png")
pairs(ResultF[,c(3,6,9,14)], pch=19, 
      labels = c("Parametric","RRPP", "FRPP","RRPP-Alt"), main = "F-X1:X2")
dev.off()

png(file="Figs/GLS.Fact-p1.png")
pairs(ResultP[,c(1,4,7,10,12)], pch=19, 
      labels = c("Parametric","RRPP", "FRPP","Restricted","RRPP-Alt"), main = "P-X1")
dev.off()

png(file="Figs/GLS.Fact-p2.png")
pairs(ResultP[,c(2,5,8,11,13)], pch=19, 
      labels = c("Parametric","RRPP", "FRPP","Restricted","RRPP-Alt"), main = "P-X2")
dev.off()

png(file="Figs/GLS.Fact-p3.png")
pairs(ResultP[,c(3,6,9,14)], pch=19, 
      labels = c("Parametric","RRPP", "FRPP","RRPP-Alt"), main = "P-X1:X2")
dev.off()

##
colSums(ifelse(ResultP <= 0.05, 1, 0))/nsims # type I error
cor(ResultP)

#####################  From 2-group
#png(file="Figs/GLS.Fact-b1.png")
#pairs(ResultB[,c(1,4,7)], pch=19, 
#      labels = c("Parametric","RRPP", "FRPP"), main = "Beta-X1")
#dev.off()

#png(file="Figs/GLS.Fact-b2.png")
#pairs(ResultB[,c(2,5,8)], pch=19, 
#      labels = c("Parametric","RRPP", "FRPP"), main = "Beta-X2")
#dev.off()

#png(file="Figs/GLS.Fact-b3.png")
#pairs(ResultB[,c(3,6,9)], pch=19, 
#      labels = c("Parametric","RRPP", "FRPP"), main = "Beta-X1:X2")
#dev.off()