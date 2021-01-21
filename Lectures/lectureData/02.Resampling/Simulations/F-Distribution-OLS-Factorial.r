##Simulation code for F-distributions for ANOVA
library(RRPP)
source('rrpp-alt.R')
source('rrpp-rest-4Gp.r')

nsim <- 100
iter <- 999
n <- 100
groups <- 4
df2 <- n-(2*(groups-1)+(groups-1)^2)-1 #full factorial

####
b1 <- 0.0
b2 <- 0.0
e <- lapply(1:nsim,function(j) rnorm(n))
x <- rep(seq(0,(groups-1)), n/groups)
x2 <- rep(seq(0,(groups-1)), each = n/groups)
y <- lapply(1:nsim, function(j) (x * b1) + (x2 * b2) + (x * b1)*(x2 * b2) + e[[j]]) 
x <- factor(x); x2 <- factor(x2)
df <- lapply(1:nsim, function(j) data.frame(x = x, x2 = x2, y = y[[j]]))

fit.RRPP <- lapply(1:nsim, function(j) lm.rrpp(y ~ x*x2, data = df[[j]], 
                     print.progress = FALSE, RRPP = TRUE, iter= iter))

fit.FRPP <- lapply(1:nsim, function(j) lm.rrpp(y ~ x*x2, data = df[[j]], 
                     print.progress = FALSE, RRPP = FALSE, iter= iter))

Fs.rrpp <- lapply(1:nsim, function(j) fit.RRPP[[j]]$ANOVA$Fs)
Fs.frpp <- lapply(1:nsim, function(j) fit.FRPP[[j]]$ANOVA$Fs)

Fs.rrpp.alt <- lapply(1:nsim, function(j) rrpp.alt(y ~ x*x2, data = df[[j]], 
                         print.progress = FALSE, RRPP = TRUE, iter= iter) )
Fs.rrpp.rest <- lapply(1:nsim, function(j) rrpp.rest(y ~ x+x2, data = df[[j]], 
                         print.progress = FALSE, RRPP = TRUE, iter= iter) )

sim.results <- list(Fs.rrpp = Fs.rrpp, Fs.rrpp.alt = Fs.rrpp.alt,
                    Fs.frpp = Fs.frpp, Fs.rrpp.rest = Fs.rrpp.rest)
#save(sim.results, file = "Results/F-dist-OLS-Fact.results.Rdata")


################# PLOTS
png(file="Figs/OLS.FDist-1.png")
#X1
df.theory <- df(density(Fs.rrpp[[1]][1,])$x, groups-1, df2)
df.rrpp.alt <- lapply(1:nsim, function(j) density(Fs.rrpp.alt[[j]][1,]))
df.rrpp.rest <- lapply(1:nsim, function(j) density(Fs.rrpp.rest[[j]][1,]))
df.rrpp <- lapply(1:nsim, function(j) density(Fs.rrpp[[j]][1,]))
df.frpp <- lapply(1:nsim, function(j) density(Fs.frpp[[j]][1,]))

plot(df.rrpp[[1]]$x, df.theory, xlab = "F", ylab = "Density", 
     main = "X1", type = "n",ylim = c(0,2), xlim = c(0, 20))
for (j in 1:nsim){ points(df.rrpp.alt[[j]]$x,df.rrpp.alt[[j]]$y,
                          type ="l", lwd=3,lty=1, col= 4 )}
for (j in 1:nsim){ points(df.rrpp.rest[[j]]$x,df.rrpp.rest[[j]]$y,
                          type ="l", lwd=3,lty=1, col= "orange" )}
for (j in 1:nsim){ points(df.frpp[[j]]$x,df.frpp[[j]]$y,
                                  type ="l", lwd=3,lty=1, col="gray")}
for (j in 1:nsim){ points(df.rrpp[[j]]$x, df.rrpp[[j]]$y, 
                                  type = "l", lwd=1,col = 2)}
points(df.rrpp[[1]]$x,df.theory,type ="l", lwd=3,lty=2)
legend("topright", c("Theory", "RRPP", "FRPP", "RRPP-alt", "Restricted"), lwd=2, lty = 1, 
       col = c(1,2,"gray",4,"orange"))
dev.off()

png(file="Figs/OLS.FDist-2.png")
#X2
df.theory <- df(density(Fs.rrpp[[1]][2,])$x, groups-1, df2)
df.rrpp.alt <- lapply(1:nsim, function(j) density(Fs.rrpp.alt[[j]][2,]))
df.rrpp.rest <- lapply(1:nsim, function(j) density(Fs.rrpp.rest[[j]][2,]))
df.rrpp <- lapply(1:nsim, function(j) density(Fs.rrpp[[j]][2,]))
df.frpp <- lapply(1:nsim, function(j) density(Fs.frpp[[j]][2,]))

plot(df.rrpp[[1]]$x, df.theory, xlab = "F", ylab = "Density", 
     main = "X2", type = "n",ylim = c(0,2), xlim = c(0, 20))
for (j in 1:nsim){ points(df.rrpp.alt[[j]]$x,df.rrpp.alt[[j]]$y,
                          type ="l", lwd=3,lty=1, col= 4 )}
for (j in 1:nsim){ points(df.rrpp.rest[[j]]$x,df.rrpp.rest[[j]]$y,
                          type ="l", lwd=3,lty=1, col= "orange" )}
for (j in 1:nsim){ points(df.frpp[[j]]$x,df.frpp[[j]]$y,
                                  type ="l", lwd=3,lty=1, col="gray")}
for (j in 1:nsim){ points(df.rrpp[[j]]$x, df.rrpp[[j]]$y, 
                                  type = "l", lwd=1,col = 2)}
points(df.rrpp[[1]]$x,df.theory,type ="l", lwd=3,lty=2)
legend("topright", c("Theory", "RRPP", "FRPP", "RRPP-alt", "Restricted"), lwd=2, lty = 1, 
       col = c(1,2,"gray",4,"orange"))
dev.off()

png(file="Figs/OLS.FDist-3.png")
#X1:X2
df.theory <- df(density(Fs.rrpp[[1]][3,])$x, 2*groups-1, df2)
df.rrpp.alt <- lapply(1:nsim, function(j) density(Fs.rrpp.alt[[j]][3,]))
df.rrpp <- lapply(1:nsim, function(j) density(Fs.rrpp[[j]][3,]))
df.frpp <- lapply(1:nsim, function(j) density(Fs.frpp[[j]][3,]))

plot(df.rrpp[[1]]$x, df.theory, xlab = "F", ylab = "Density", 
     main = "X1:X2", type = "n",ylim = c(0,2), xlim = c(0, 20))
for (j in 1:nsim){ points(df.rrpp.alt[[j]]$x,df.rrpp.alt[[j]]$y,
                          type ="l", lwd=3,lty=1, col= 4 )}
for (j in 1:nsim){ points(df.frpp[[j]]$x,df.frpp[[j]]$y,
                                  type ="l", lwd=3,lty=1, col="gray")}
for (j in 1:nsim){ points(df.rrpp[[j]]$x, df.rrpp[[j]]$y, 
                                  type = "l", lwd=1,col = 2)}
points(df.rrpp[[1]]$x,df.theory,type ="l", lwd=3,lty=2)
legend("topright", c("Theory", "RRPP", "FRPP", "RRPP-alt"), lwd=2, lty = 1, 
       col = c(1,2,"gray",4))
dev.off()
