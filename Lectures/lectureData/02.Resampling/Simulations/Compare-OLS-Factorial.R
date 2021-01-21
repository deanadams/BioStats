## OLS Simulations comparing FRPP, RRPP, RRPP-alt, rrpp-restricted
  #rerun with different error distributions ala Anderson terBraak 2003

library(RRPP)
source('rrpp-alt.R')
source('rrpp-rest.r')

nsims <- 100
iter <- 999
n <- 100
betas <- seq(0, 1, 0.2)
groups <- 2

result <- t(sapply(1:nsims, function(j){
  b1 <- 0.0
  b2 <- 0.0 
  e <- rnorm(n) # e <- runif(n) # e <- rlnorm(n) # e <- rexp(n,rate = 3)
  x <- rep(seq(0,(groups-1)), n/groups)
  x2 <- rep(seq(0,(groups-1)), each = n/groups)
  y <- (x * b1) + (x2 * b2) + (x * b1)*(x2 * b2) + e 
  x <- factor(x); x2 <- factor(x2)
  df <- data.frame(x = x, x2 = x2, y = y)

  fit.param <- lm(y ~ x*x2, data = df)
  fit.RRPP <- lm.rrpp(y ~ x*x2, data = df, 
                    print.progress = FALSE, RRPP = TRUE, iter= iter)
  fit.FRPP <- lm.rrpp(y ~ x*x2, data = df, 
                       print.progress = FALSE, RRPP = FALSE, iter= iter)
  fit.RRPP.alt <- rrpp.alt(y ~ x*x2, data = df, 
                           print.progress = FALSE, RRPP = TRUE, iter= iter)
  fit.rest <- rrpp.rest(y ~ x+x2, data = df, SS.type = "II",
                        print.progress = FALSE, RRPP = TRUE, iter= iter)

  aov.param <- anova(fit.param)
  aov.RRPP <- anova(fit.RRPP)
  aov.FRPP <- anova(fit.FRPP)

  cat("iter", j, " ")
  c(b.par = fit.param$coefficients[-1],
    b.RRPP = fit.RRPP$LM$coefficients[-1],
    b.FRPP = fit.FRPP$LM$coefficients[-1],
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

##save(result, file = "Results/Results.OLSFact.Rdata")

ResultB <- data.frame(result[,c(1:9)])
ResultF <- result[,c(10:23)]
ResultP <- result[,c(24:37)]

#####

png(file="Figs/OLS.Fact-f1.png")
pairs(ResultF[,c(1,4,7,10,12)], pch=19, 
      labels = c("Parametric","RRPP", "FRPP","Restricted","RRPP-Alt"), main = "F-X1")
dev.off()

png(file="Figs/OLS.Fact-f2.png")
pairs(ResultF[,c(2,5,8,11,13)], pch=19, 
      labels = c("Parametric","RRPP", "FRPP","Restricted","RRPP-Alt"), main = "F-X2")
dev.off()

png(file="Figs/OLS.Fact-f3.png")
pairs(ResultF[,c(3,6,9,14)], pch=19, 
      labels = c("Parametric","RRPP", "FRPP","RRPP-Alt"), main = "F-X1:X2")
dev.off()

png(file="Figs/OLS.Fact-p1.png")
pairs(ResultP[,c(1,4,7,10,12)], pch=19, 
      labels = c("Parametric","RRPP", "FRPP","Restricted","RRPP-Alt"), main = "P-X1")
dev.off()

png(file="Figs/OLS.Fact-p2.png")
pairs(ResultP[,c(2,5,8,11,13)], pch=19, 
      labels = c("Parametric","RRPP", "FRPP","Restricted","RRPP-Alt"), main = "P-X2")
dev.off()

png(file="Figs/OLS.Fact-p3.png")
pairs(ResultP[,c(3,6,9,14)], pch=19, 
      labels = c("Parametric","RRPP", "FRPP","RRPP-Alt"), main = "P-X1:X2")
dev.off()

##
colSums(ifelse(ResultP <= 0.05, 1, 0))/nsims # type I error
cor(ResultP)

## ### FOR 2 Group case only!
#png(file="Figs/OLS.Fact-b1.png")
#pairs(ResultB[,c(1,4,7)], pch=19, 
#      labels = c("Parametric","RRPP", "FRPP"), main = "Beta-X1")
#dev.off()

#png(file="Figs/OLS.Fact-b2.png")
#pairs(ResultB[,c(2,5,8)], pch=19, 
#      labels = c("Parametric","RRPP", "FRPP"), main = "Beta-X2")
#dev.off()

#png(file="Figs/OLS.Fact-b3.png")
#pairs(ResultB[,c(3,6,9)], pch=19, 
#      labels = c("Parametric","RRPP", "FRPP"), main = "Beta-X1:X2")
#dev.off()