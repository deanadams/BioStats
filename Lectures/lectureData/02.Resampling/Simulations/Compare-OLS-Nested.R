## OLS Nested model (with additional factor)
   #rerun with different error distributions ala Anderson terBraak 2003

library(RRPP)
source('rrpp-alt-nest.R')

nsims <- 100
iter <- 999
n <- 100
g1 <- 5
g2 <- 10
g3 <- 2

result <- t(sapply(1:nsims, function(j){
  b1 <- 0.0
  b2 <- 0.0 
  b3 <- 0.0
#  e <- rnorm(n) # e <- runif(n) # e <- rlnorm(n) # e <- rexp(n,rate = 3)
  x <- sample(rep(seq(1,g1),n/g1))
  x2 <- c(rep(seq(0,(n/4-1)),  each = 2),rep(seq(0,(n/4-1)),  each = 2))      
  x3 <- rep(seq(0,(g3-1)), each = n/g3) 
  y <- (x * b1) + (x2 * b2) + (x3 * b3) + e 
  x <- factor(x); x2 <- factor(x2); x3 <- factor(x3)
  df <- data.frame(x = x, x2 = x2, x3 = x3, y = y)

  fit.param <- lm(y ~ x + x2/x3, data = df)
  fit.RRPP <- lm.rrpp(y ~ x + x2/x3, data = df, print.progress = FALSE, SS.type = "III")
  fit.FRPP <- lm.rrpp(y ~ x + x2/x3, data = df, 
                      RRPP = FALSE, print.progress = FALSE, SS.type = "III")

  fit.RRPP.alt <- rrpp.alt(y ~ x + x2/x3, data = df, SS.type = "III",
                           print.progress = FALSE, RRPP = TRUE, iter= iter)

  aov.RRPP <-  anova(fit.RRPP, error = c("Residuals", 
                                         "x2:x3", "Residuals"))
  aov.FRPP <- anova(fit.FRPP, error = c("Residuals", 
                                        "x2:x3", "Residuals"))
  aov.param <- NULL
  aov.param[1] <- 1-pf(aov.RRPP$table[1,5],aov.RRPP$table[1,1], aov.RRPP$table[4,1])
  aov.param[2] <- 1-pf(aov.RRPP$table[2,5],aov.RRPP$table[2,1], aov.RRPP$table[3,1])
  aov.param[3] <- 1-pf(aov.RRPP$table[3,5],aov.RRPP$table[3,1], aov.RRPP$table[4,1])

  cat("iter", j, " ")
  c(F.par = aov.RRPP$table$F[1:3],  # confirmed identical to lm as per above 
    F.RRPP = aov.RRPP$table$F[1:3], 
    F.FRPP = aov.FRPP$table$F[1:3],
    F.RRPP.tBr = fit.RRPP.alt[,1],
    p.par = aov.param, 
    p.RRPP = aov.RRPP$table$`Pr(>F)`[1:3], 
    p.FRPP = aov.FRPP$table$`Pr(>F)`[1:3],
    p.RRPP.tBr =  apply(fit.RRPP.alt,1,RRPP:::pval)
  )
}))

#save(result, file = "Results/Results.OLSNest.Rdata")

ResultF <- result[,c(1:12)]
ResultP <- result[,c(13:24)]

#####
png(file="Figs/OLS.Nest-f1.png")
pairs(ResultF[,c(1,4,7,10)], pch=19, 
      labels = c("Parametric","RRPP", "FRPP","RRPP-Alt"), main = "F-X1")
dev.off()

png(file="Figs/OLS.Nest-f2.png")
pairs(ResultF[,c(2,5,8,11)], pch=19, 
      labels = c("Parametric","RRPP", "FRPP","RRPP-Alt"), main = "F-X2")
dev.off()

png(file="Figs/OLS.Nest-f3.png")
pairs(ResultF[,c(3,6,9,12)], pch=19, 
      labels = c("Parametric","RRPP", "FRPP","RRPP-Alt"), main = "F-X2:X3")
dev.off()

png(file="Figs/OLS.Nest-p1.png")
pairs(ResultP[,c(1,4,7,10)], pch=19, 
      labels = c("Parametric","RRPP", "FRPP","RRPP-Alt"), main = "P-X1")
dev.off()

png(file="Figs/OLS.Nest-p2.png")
pairs(ResultP[,c(2,5,8,11)], pch=19, 
      labels = c("Parametric","RRPP", "FRPP","RRPP-Alt"), main = "P-X2")
dev.off()

png(file="Figs/OLS.Nest-p3.png")
pairs(ResultP[,c(3,6,9,12)], pch=19, 
      labels = c("Parametric","RRPP", "FRPP","RRPP-Alt"), main = "P-X2:X3")
dev.off()

##
colSums(ifelse(ResultP <= 0.05, 1, 0))/nsims # type I error
cor(ResultP)
