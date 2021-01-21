#####  Interaction Term exploration: Trajectory Analysis	
  #Packages: geomorph, RRPP
library(geomorph)
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

#Pairwise comparisons while accounting for a covariate
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
class(Pupfish) <- "geomorph.data.frame" # Help R and geomorph work with the data

TA <- trajectory.analysis(f1 = coords ~ Pop * Sex, data = Pupfish, print.progress = TRUE)
summary(TA, angle.type = "deg")
plot(TA)

#2: Compare groups of pre-existing trajectories
data(motionpaths)
gdf <- geomorph.data.frame(trajectories = motionpaths$trajectories,
                           groups = motionpaths$groups)
TA <- trajectory.analysis(f1 = trajectories ~ groups, 
                          traj.pts = 5, data=gdf,print.progress=FALSE)
summary(TA)
plot(TA)


