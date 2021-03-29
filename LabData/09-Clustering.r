###Some Ordination Methods
   #Packages: geomorph (for data)

rm(list=ls())
#######################
#  Some Clustering Methods
#	NOTE: phylogenetic clustering approaches found in package APE, Phangorn & others
#	  see: 'Analysis of phylogenetic and evolution in R. 2012. E. Paradis
##
#######################
rm(list=ls())
mole.data<-read.csv("Data/Lab-09.Moles.csv",header=T,row.names=1)
mole.dist<-as.dist(mole.data)

### PCoA of mole data
PCoA<-cmdscale(mole.dist)
plot(PCoA,pch=21,bg='black',cex=1.5,asp=1)
  text(PCoA[,1]+1.5,PCoA[,2],row.names(mole.data))

#Some clustering methods
mole.single<-hclust(mole.dist,method="single")       #Single-link
mole.complete<-hclust(mole.dist,method="complete")   #Complete-link
mole.upgma<-hclust(mole.dist,method="average")       #UPGMA = average-link
mole.upgmc<-hclust(mole.dist,method="centroid")      #UPGMC
mole.wpgma<-hclust(mole.dist,method="mcquitty")      #WPGMA
mole.wpgmc<-hclust(mole.dist,method="median")        #WPGMC
mole.wards<-hclust(mole.dist,method="ward.D")          #Ward's

##PLOTS
plot(mole.single, hang=-1,lwd=2, main="Single Linkage")
plot(as.dendrogram(mole.single),horiz=TRUE,lwd=4,xlim=c(16,-1), main="Single Linkage")  #single-link
plot(as.dendrogram(mole.complete),horiz=TRUE,lwd=4,xlim=c(16,-1), main="Complete Linkage")  #complete-link
plot(as.dendrogram(mole.upgma),horiz=TRUE,lwd=4,xlim=c(16,-1), main="UPGMA")  #UPGMA
plot(as.dendrogram(mole.wpgma),horiz=TRUE,lwd=4,xlim=c(16,-1), main="WPGMA")  #WPGMA
plot(as.dendrogram(mole.upgmc),horiz=TRUE,lwd=4,xlim=c(16,-1), main="UPGMC")  #UPGMC
plot(as.dendrogram(mole.wpgmc),horiz=TRUE,lwd=4,xlim=c(16,-1), main="WPGMC")  #WPGMC
plot(as.dendrogram(mole.wards),horiz=TRUE,lwd=4,xlim=c(16,-1), main="Ward's")  #Ward's

cophenetic(mole.upgma)

plot(mole.dist,cophenetic(mole.upgma))  #NOTE that small distances better preserved
plot(mole.dist,dist(PCoA[,1:2]))  #PCoA

#some data (head shape in Plethodon salamanders) 
library(geomorph)
data(plethodon)
shape <- two.d.array(gpagen(plethodon$land)$coords) #data prep
  PC.scores<-prcomp(shape)$x
plot(PC.scores,pch=21,bg=as.factor(paste(plethodon$species,plethodon$site)),asp=1)

##UPGMA
pleth.dist<-dist(shape)  #choose your distance
pleth.upgma<-hclust(pleth.dist,method="average") 
plot(as.dendrogram(pleth.upgma),horiz=TRUE,lwd=4, main="UPGMA")  #UPGMA

#PLOT of actual vs. UPGMA distances
plot(pleth.dist,cophenetic(pleth.upgma))

# SAME from PC
plot(pleth.dist,dist(PC.scores[,1:2]))

#K-means
kclusters4<-kmeans(PC.scores,4)
plot(PC.scores[,1:2],col=kclusters4$cluster,, main="K=4")
points(kclusters4$centers, col = 1:4, pch = 8, cex=2)

kclusters3<-kmeans(PC.scores,3)
plot(PC.scores[,1:2],col=kclusters3$cluster, main="K=3")
points(kclusters3$centers, col = 1:3, pch = 8, cex=2)

kclusters2<-kmeans(PC.scores,2)
plot(PC.scores[,1:2],col=kclusters2$cluster, main="K=2")
points(kclusters2$centers, col = 1:2, pch = 8, cex=2)

#NOTE: repeating k-means at a given level can lead to differing results

#compare TESS
TESS<-array(NA,40)
for (i in 1:40){
  TESS[i]<-kmeans(PC.scores,i)$tot.withinss
}
plot( TESS)  #seems to bottom out at 3 groups

plot(PC.scores[,1:2],col=kclusters3$cluster, main="K=3")
points(kclusters3$centers, col = 1:3, pch = 8, cex=2)
