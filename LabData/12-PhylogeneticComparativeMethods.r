# Some Phylogenetic Comparative Methods: PCMs
   ## packages: nlme, ape, geiger, phytools, geomorph

library(nlme)	# Contains GLS 
library(ape)	# Many phylogenetics-related analyses (see Paradis, Analysis of Phylogenetics and Evolution with R)
library(geiger) # Many evolutionary analyses
library(phytools) 
library(geomorph) #contains some multivariate PCMs
library(RRPP)

## Read data, phylogeny, and match the two
tree.best<-read.nexus("Data/PlethodontidTree.nex") #Maximum Credible Tree
plethdata<-read.csv("Data/meandata-CinGlutOnly.csv",row.names=1, header=TRUE )
size<-plethdata[,2];names(size)<-row.names(plethdata) 
gp<-as.factor(plethdata[,1]); names(gp)<-row.names(plethdata)
# Prune tree to match taxa
plethtree<-treedata(phy = tree.best,data = plethdata, warnings = FALSE)$phy
plot(tree.best)
plot(plethtree)
axisPhylo(1) 
#Relative body proportions: trait/size
Y<-plethdata[,c(3:8)]
Y<-apply(Y,2,as.numeric);row.names(Y)<-row.names(plethdata)
TL <- Y[,1]
Y.sz<-Y/size
gdf <- rrpp.data.frame(Y=Y, Y.sz=Y.sz, TL = TL, size=size,gp=gp) 
gdf$Cov <- vcv.phylo(plethtree)

###################
### Analyses
   #Regression: non-phylogenetic
anova(lm(TL~size))  
summary(lm(TL~size))
plot(TL~size)
abline(lm(TL~size))
summary(gls(TL~size, data=data.frame(TL=TL,size=size)))

# 1: PhylogeneticRegression (PGLS)
   #using GLS
df <- data.frame(TL=TL,size=size, species = names(TL), gp=gp)
bm.gls<-gls(TL~size, correlation=corBrownian(phy=plethtree, form = ~species), data= df)
summary(bm.gls)  #Here the correlation structure of the phylogeny is used
anova(bm.gls)

#using RRPP: same
pgls.res<-lm.rrpp(TL~size,data=gdf, Cov = gdf$Cov,  print.progress = FALSE) 
anova(pgls.res)
pgls.res$LM$gls.coefficients

#using Independent Contrasts: same
picTL<-pic(TL, plethtree)
picSz<-pic(size, plethtree)
cor.test(picTL, picSz)
summary(lm(picTL~picSz + 0)) #Contrasts run through origin: see Garland et al. 1992
plot(picTL~picSz)
abline(lm(picTL~picSz + 0))

#Phylogenetic anova
anova(lm.rrpp(TL ~ gp, Cov = gdf$Cov, data = gdf, iter = 999, print.progress = FALSE) )
anova(gls(TL~gp, correlation=corBrownian(phy=plethtree, form = ~species), data= df))  #same


  #multivariate phy-anova/regression (even when p>N)
PGLS.reg<-lm.rrpp(Y ~ size, Cov = gdf$Cov, data = gdf, iter = 999, print.progress = FALSE) 
anova(PGLS.reg)
plot(PGLS.reg, type = "regression", predictor=size,reg.type = "RegScore")

PGLS.aov<-lm.rrpp(Y ~ gp, Cov = gdf$Cov, data = gdf, iter = 999, print.progress = FALSE) 
anova(PGLS.aov)  

# 2: Phylogenetic PLS: multivariate
PLS.Y <- phylo.integration(A = Y[,1:3], A2 = Y[,4:6], phy= plethtree, print.progress = FALSE)
summary(PLS.Y)
plot(PLS.Y)

# 3: Phylogenetic ordination
  #Note: copied tree w/ simplified names for ease of visualization
   #one should retain the original names in one's dataset!
rownames(Y.sz) <- 1:nrow(Y.sz)
plethtree2 <- plethtree; plethtree2$tip.label <- rownames(Y.sz)
PCA.w.phylo <- gm.prcomp(Y.sz, phy = plethtree2) #: phylomorphospace
phylo.PCA <- gm.prcomp(Y.sz, phy = plethtree2, GLS = TRUE) #phylo-PCA
PaCA.gls <- gm.prcomp(Y.sz, phy = plethtree2, 
                      align.to.phy = TRUE, GLS = TRUE) # PACA

par(mfrow=c(2,2))
plot(PCA.w.phylo, phylo = TRUE, main = "PCA.w.phylo", 
     phylo.par = list(node.labels = FALSE))
plot(phylo.PCA, phylo = TRUE, main = "phylo PCA",
     phylo.par = list(node.labels = FALSE))
plot(PaCA.gls, phylo = TRUE, main = "PaCA",
     phylo.par = list(node.labels = FALSE))
par(mfrow=c(1,1))


############
# 4: Phylogenetic signal
phylosig(plethtree, TL, method="K", test=T, nsim=1000)  #phytools
physignal(gdf$TL,plethtree, print.progress = FALSE)   #geomorph

#Multivariate
PS.shape <- physignal(A=Y,phy=plethtree,iter=999, print.progress = FALSE)
summary(PS.shape)
plot(PS.shape)

#5a: Compare Evolutionary Rates Among Clades
ER<-compare.evol.rates(A=Y, phy=plethtree,gp=gp,iter=999, print.progress = FALSE)
summary(ER)   #significantly higher rate of morphological evolution 'large' Plethodon
plot(ER)

#5b: Compare Evolutionary Rates Among Traits (multivariate)
var.gp <- c("B","A","A","B","B","B")  #head (A) vs. body (B)
EMR<-compare.multi.evol.rates(A=Y,gp=var.gp, Subset=TRUE, phy= plethtree,iter=999, print.progress = FALSE)
summary(EMR) #Limb traits evolve faster than head traits
plot(EMR)

# 6: Comparing evolutionary models

## Compare BM vs. OU models using GEIGER
geo=get(data(geospiza))      #a smaller dataset (for time reasons)
tmp=treedata(geo$phy, geo$dat)
phy=tmp$phy; dat=tmp$data[,"tarsusL"] 
plot(phy)

bm.M<-fitContinuous(phy, dat, model="BM")
ou.M<-fitContinuous(phy, dat, bounds = list(alpha = c(min = exp(-500), max = exp(2))),model="OU")
  bm.M
  ou.M

# Compare models: LRT & AIC
LRT.M<-(2*(ou.M$opt$lnL-bm.M$opt$lnL))
prob.M<-pchisq(LRT.M, 1, lower.tail=FALSE)
LRT.M
prob.M

bm.M$opt$aic
ou.M$opt$aic #OU does not provide a better fit: use simpler model (BM)

# NOTE: The package OUCH, MVSLOUCH, OUwie, and otherscan compare more complicated models: BM1, BMM, OU1, OUM, etc.).
