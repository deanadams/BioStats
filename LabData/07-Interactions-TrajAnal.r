#####  Interaction Term exploration: Trajectory Analysis	
  #Packages: RRPP
library(RRPP)

##Interaction terms and sub-components

#Pairwise comparisons from MANOVA
data(Pupfish)
Pupfish$Group <- interaction(Pupfish$Sex, Pupfish$Pop)
fit.m <- lm.rrpp(coords ~ Pop * Sex, data = Pupfish, print.progress = FALSE) 
anova(fit.m)$table
summary(pairwise(fit.m,groups = Pupfish$Group))
plot(fit.m, type = "PC", pch=21, bg = Pupfish$Group, cex=2)
legend("topright", levels(Pupfish$Group), pch = 21, pt.bg = 1:4)

#Pairwise group comparisons while accounting for a covariate
Pupfish$logSize <- log(Pupfish$CS)
fit.slopes <- lm.rrpp(coords ~ logSize * Pop * Sex, data = Pupfish, print.progress = FALSE)
anova(fit.slopes)$table
PWS.gp <- pairwise(fit.slopes, groups = Pupfish$Group)
summary(PWS.gp)
plot(fit.slopes, type = "PC", pch=21, bg = Pupfish$Group, cex=2)
legend("topright", levels(Pupfish$Group), pch = 21, pt.bg = 1:4)


#Comparison of slopes and lengths in MANCOVA
PWS <- pairwise(fit.slopes, groups = Pupfish$Group, covariate = Pupfish$logSize)
summary(PWS)
summary(PWS, test.type = "VC", angle.type = "deg")

plot(fit.slopes, type = "regression", reg.type = "RegScore", pch=21, bg = Pupfish$Group, predictor = Pupfish$logSize, cex=2)
legend("topleft", levels(Pupfish$Group), pch = 21, pt.bg = 1:4)

plot(fit.slopes, type = "regression", reg.type = "PredLine", pch=21, bg = Pupfish$Group, predictor = Pupfish$logSize, cex=2)
legend("topright", levels(Pupfish$Group), pch = 21, pt.bg = 1:4)


##### Trajectory Analysis
#1: Estimate trajectories from LS means in factorial mode
fit <- lm.rrpp(coords ~ Pop * Sex, data = Pupfish, iter = 999)
TA <- trajectory.analysis(fit, groups = Pupfish$Pop, 
                          traj.pts = Pupfish$Sex, print.progress = FALSE)

# Magnitude difference (absolute difference between path distances)
summary(TA, attribute = "MD") 

# Correlations (angles) between trajectories
summary(TA, attribute = "TC", angle.type = "deg") 

# Plot results
TP <- plot(TA, pch = as.numeric(Pupfish$Pop) + 20, bg = as.numeric(Pupfish$Sex),
           cex = 1, col = "black")
add.trajectories(TP, traj.pch = c(21, 22), start.bg = 1, end.bg = 2)
legend("topright", levels(Pupfish$Pop), pch =  c(21, 22), pt.bg = 1)


#2: Compare groups of pre-existing trajectories
data(motionpaths)
fit <- lm.rrpp(trajectories ~ groups, data = motionpaths, iter = 999)
TA <- trajectory.analysis(fit, groups = motionpaths$groups, traj.pts = 5)

# Magnitude difference (absolute difference between path distances)
summary(TA, attribute = "MD") 

# Correlations (angles) between trajectories
summary(TA, attribute = "TC", angle.type = "deg") 

# Shape differences between trajectories 
summary(TA, attribute = "SD") 

TP <- plot(TA, pch = 21, bg = as.numeric(motionpaths$groups),
           cex = 0.7, col = "gray")
add.trajectories(TP, traj.pch = 21, traj.bg = 1:4)


