##Simulation code for F-distributions for Nested ANOVA
library(phytools)
library(RRPP)
source('anova.nested.source.r')
source('rrpp-alt-nest.R')

nsim <- 100
iter <- 999
n <- 100
g1 <- 5
g2 <- 10
g3 <- 2

#### 
b1 <- 0.0
b2 <- 0.0 
b3 <- 0.0
tree <- pbtree (n = n,nsim = nsim)
Cov <- lapply(1:nsim, function(j) vcv.phylo(tree[[j]]))
e <- lapply(1:nsim,function(j) fastBM(tree[[j]]))

x <- lapply(1:nsim, function(j) sample(rep(seq(1,g1),n/g1)) )
x2 <- lapply(1:nsim, function(j) c(rep(seq(0,(n/4-1)),  each = 2),rep(seq(0,(n/4-1)),  each = 2)) )
x3 <- lapply(1:nsim, function(j) rep(seq(0,(g3-1)), each = n/g3))
for (i in 1:nsim){ names(x[[i]]) <- names(x2[[i]]) <- names(x3[[i]]) <- names(e[[i]]) }

y <- lapply(1:nsim, function(j) (x[[j]] * b1) + (x2[[j]] * b2) + 
              (x[[j]] * b1)*(x2[[j]] * b2) + (x3[[j]] * b3) + e[[j]] )
for (i in 1:nsim){
  x[[i]] <- factor(x[[i]]); x2[[i]] <- factor(x2[[i]]); x3[[i]] <- factor(x3[[i]])
}
df <- lapply(1:nsim, function(j) data.frame(x = x[[j]], x2 = x2[[j]], x3 = x3[[j]],
             y = y[[j]]))

fit.RRPP <- lapply(1:nsim, function(j) lm.rrpp(y ~ x + x2/x3, data = df[[j]], 
            Cov = Cov[[j]], print.progress = FALSE, RRPP = TRUE, SS.type = "III"))

fit.FRPP <- lapply(1:nsim, function(j) lm.rrpp(y ~ x + x2/x3, data = df[[j]], 
            Cov = Cov[[j]], print.progress = FALSE, RRPP = FALSE, SS.type = "III"))

anova(fit.FRPP[[1]])
#Grab Nested F-values
Fs.rrpp <- lapply(1:nsim, function(j) aov.single.model.F(fit.RRPP[[j]],
              error = c("Residuals", "x2:x3", "Residuals")))

Fs.frpp <- lapply(1:nsim, function(j) aov.single.model.F(fit.RRPP[[j]],
           error = c("Residuals", "x2:x3", "Residuals")))

Fs.rrpp.alt <- lapply(1:nsim, function(j) rrpp.alt(y ~ x + x2/x3, data = df[[j]], SS.type = "III",
                Cov = Cov[[j]], print.progress = FALSE, RRPP = TRUE, iter= iter))

sim.results <- list(Fs.rrpp = Fs.rrpp, Fs.rrpp.alt = Fs.rrpp.alt,
                    Fs.frpp = Fs.frpp)
#save(sim.results, file = "Results/F-dist-GLS-Nest.results.Rdata")

################# PLOTS
png(file="Figs/GLS-Nest.FDist-1.png")
#X1
df.theory <- df(density(Fs.rrpp[[1]][1,])$x, g1, 46)
df.rrpp <- lapply(1:nsim, function(j) density(Fs.rrpp[[j]][1,]))
df.rrpp.alt <- lapply(1:nsim, function(j) density(Fs.rrpp.alt[[j]][1,]))
df.frpp <- lapply(1:nsim, function(j) density(Fs.frpp[[j]][1,]))

plot(df.rrpp[[1]]$x, df.theory, xlab = "F", ylab = "Density", 
     main = "X1", type = "n",ylim = c(0,2), xlim = c(0, 20))
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

png(file="Figs/GLS-Nest.FDist-2.png")
#X2
df.theory <- df(density(Fs.rrpp[[1]][2,])$x, 24, 25)
df.rrpp <- lapply(1:nsim, function(j) density(Fs.rrpp[[j]][2,]))
df.rrpp.alt <- lapply(1:nsim, function(j) density(Fs.rrpp.alt[[j]][2,]))
df.frpp <- lapply(1:nsim, function(j) density(Fs.frpp[[j]][2,]))

plot(df.rrpp[[1]]$x, df.theory, xlab = "F", ylab = "Density", 
     main = "X2", type = "n",ylim = c(0,2), xlim = c(0, 20))
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

png(file="Figs/GLS-Nest.FDist-3.png")
#X2:X3
df.theory <- df(density(Fs.rrpp[[1]][3,])$x, 25, 46)
df.rrpp.alt <- lapply(1:nsim, function(j) density(Fs.rrpp.alt[[j]][3,]))
df.rrpp <- lapply(1:nsim, function(j) density(Fs.rrpp[[j]][3,]))
df.frpp <- lapply(1:nsim, function(j) density(Fs.frpp[[j]][3,]))

plot(df.rrpp[[1]]$x, df.theory, xlab = "F", ylab = "Density", 
     main = "X2:X3", type = "n",ylim = c(0,3), xlim = c(0, 20))
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