###Some Ordination Methods
      #Packages: geomorph [for data], vegan, MASS

########################
#	PCA: Principal Components Analysis
library(RRPP)
library(vegan)
bumpus<-read.csv("data/bumpus.csv",header=T)
Y<-bumpus[,5:12]
Y <- scale(Y, scale = FALSE) #center data
gp.bumpus <- as.factor(bumpus$sex)
pca.bumpus<-prcomp(Y) 
summary(pca.bumpus)
PC.scores<-pca.bumpus$x
pca.bumpus$rotation[,1]  #1st PC axis only

plot(PC.scores,xlab="PC I", ylab="PC II",asp=1,pch=21,bg=gp.bumpus,cex = 1.5)
legend("topright", levels(gp.bumpus), pch = 21,pt.bg=1:2)

#Eigenvalue analysis via broken stick model
screeplot(pca.bumpus,bstick = TRUE)  #implies 2 PCAs sufficient graphical representation

##Plot of actual vs. PC1-2 distances
plot(dist(Y),dist(PC.scores[,1:2]))
cor(dist(Y),dist(PC.scores[,1:2]))

### PCA via svd
svd.res<-svd(Y)  
svd.res$d^2/sum(svd.res$d^2)   #same % variation per PC axis
  pc.scores.svd<-svd.res$u%*%diag(svd.res$d)  #PCA scores
plot(pc.scores.svd,asp=1, pch=21,bg=gp.bumpus,cex = 1.5)
legend("topright", levels(gp.bumpus), pch = 21,pt.bg=1:2)

#### PCA "by hand" via eigen-analysis
vcv.bumpus<-var(Y)	#Calculate PC axes
pc.bumpus<-eigen(vcv.bumpus) 
pc.bumpus$values/sum(pc.bumpus$values)   #same % variation per PC axis
pc.scores<-Y%*%pc.bumpus$vectors	#Projection
plot(pc.scores,xlab="PC I", ylab="PC II",asp=1,pch=21,bg=gp.bumpus,cex = 1.5)
legend("topright", levels(gp.bumpus), pch = 21,pt.bg=1:2)

    #note, axes can be reflected and still retain correct information  

#############Biplot with PCA
biplot(pca.bumpus)

###########################
# 	PCoA
bumpus.dist<-dist(Y)
PCoA<-cmdscale(bumpus.dist)   #from vegan
plot(-1*PCoA[,1], PCoA[,2],pch=21,bg=gp.bumpus,cex=1.5,asp=1)
legend("topright", levels(gp.bumpus), pch = 21,pt.bg=1:2)
#############################
#	NMDS
bumpus.nmds <- metaMDS(bumpus.dist, autotransform=FALSE, k=2)
    #A nice function. Runs 20 times with different starting points; tries to avoid local optima
plot(bumpus.nmds$points, asp=1,pch=21, bg=gp.bumpus, cex=1.5,xlab="NMDS1", ylab="NMDS2")
legend("topright", levels(gp.bumpus), pch = 21,pt.bg=1:2)

  #Plot of actual to plot distances (often curved)
plot(bumpus.dist, dist(scores(bumpus.nmds, display='sites'), method='eucl'), xlab='D.obs', ylab='D.plot')

data(dune)  #from vegan
dune
dune.dist<-vegdist(dune)  #default = Bray-Curtis distance
dune.nmds <- metaMDS(dune.dist, autotransform=FALSE, k=2)
#A nice function. Runs 20 times with different starting points; tries to avoid local optima
plot(dune.nmds)
plot(dune.nmds,type='t')

#Plot of actual to plot distances (often curved)
plot(dune.dist, dist(scores(dune.nmds, display='sites'), method='eucl'), xlab='D.obs', ylab='D.plot')


#################
#Correspondence Analysis
dune.cca<-cca(dune)  #from vegan
plot(dune.cca)

#Detrended Correspondence Analysis
dune.dca<-decorana(dune)
plot(dune.dca)

