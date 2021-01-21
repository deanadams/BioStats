### Multivariate Association and Canonical Ordination Methods
    #Packages: geomorph, vegan, MASS

#Multivariate Association
library(geomorph)
library(vegan)
library(RRPP)

data(pupfish)
p1<-c(4,10:17, 39:56) #variables [landmarks] in Y1
Y1<-two.d.array(pupfish$coords[p1,,]) # head data as 2D matrix
Y2<-two.d.array(pupfish$coords[-p1,,]) # body data as 2D matrix
Group<-as.factor(paste(pupfish$Pop,pupfish$Sex)) 
cols <- rep(1, 56)
cols[p1] <- 2
plotAllSpecimens(pupfish$coords, mean = FALSE, plot.param=list(pt.bg = cols)) #to show the data


## 1: Mantel
mantel(dist(Y1),dist(Y2),permutations = 9999) 
 plot(dist(Y1),dist(Y2))
mantel.partial(dist(Y1),dist(Y2),dist(pupfish$CS))  #3-way Mantel

#2: RV coefficient
RV<-function(x,y,iter=999){
  n<-nrow(x)
  x<-scale(x,scale=FALSE); y<-scale(y,scale=FALSE)
  ind <- c(list(1:n), (Map(function(x) sample.int(n, n), 1:iter)))
  y.rand <- lapply(1:(iter+1), function(i) y[ind[[i]], ])
  S11<-crossprod(x, x)/(n - 1)
  S22<-crossprod(y, y)/(n - 1)
  S12.r<-lapply(1:(iter+1), function(i) crossprod(x, y.rand[[i]])/(n - 1))
  RV.r<-unlist(lapply(1:(iter+1), function(i) sum(diag(S12.r[[i]]%*%t(S12.r[[i]]))) /
                 sqrt(sum(diag(S11%*%S11)) * sum(diag(S22%*%S22)))))
  p.val<- 1-(rank(RV.r)[1])/(iter+1)
  p.val<-ifelse(p.val==0,1/(iter+1),p.val)
  return(list(RV=RV.r[[1]],pvalue=p.val))
}

RV(Y1,Y2)
sqrt(RV(Y1,Y2)$RV)  #closer to a correlation

#3: PLS
pls.res<-two.b.pls(Y1,Y2)
summary(pls.res)
plot(pls.res)


######################
#Canonical Ordination

# Read Data
mydata<-read.csv("data/Lab-10.Pupfish.csv",header=T)
species<-as.factor(mydata[,1]); sex<-as.factor(mydata[,2])
SL<-(mydata[,3]); Y<-as.matrix(mydata[,-(1:3)])	
Group<-as.factor(paste(species,sex))  


Y<-prcomp(Y)$x	
rdf <- rrpp.data.frame(Y=Y, SL=SL, sex=sex, species=species)
 col.gp<-rep("green",nrow(Y));   col.gp[which(species== 'FW')]<-"red"
 shape.gp<-rep(21,nrow(Y));   shape.gp[which(sex== 'M')]<-22
plot(Y,pch=shape.gp,bg=col.gp,asp=1,cex=1.5)

## CVA
library(MASS)
lda.pupfish<-lda(Y,Group)
cva.pupfish<-predict(lda.pupfish,Y)
cv.scores<-cva.pupfish$x

plot(cv.scores,pch=shape.gp,bg=col.gp,asp=1,cex=1.5)

#### CVA ISSUES:

#1: equidistant groups are not represented as such
library(mvtnorm) 
corr.val<-0.7
groups<-gl(3,50)
a<- rmvnorm(n=50,mean=c(-3.4,0),sigma=matrix(c(1,corr.val,corr.val,1),2,2))
b<- rmvnorm(n=50,mean=c(3.4,0),sigma=matrix(c(1,corr.val,corr.val,1),2,2))
c<- rmvnorm(n=50,mean=c(0,6),sigma=matrix(c(1,corr.val,corr.val,1),2,2))
orig.data<-rbind(a,b,c)
col.gp.r<-rep("black",nrow(orig.data)); col.gp.r[which(groups== '2')]<-"red"; col.gp.r[which(groups== '3')]<-"green"

plot(orig.data,pch=21,bg=col.gp.r,asp=1,cex=1.5)
ordiellipse(orig.data, groups,conf=0.95)
plot(predict(lda(orig.data,groups))$x,pch=21,bg=col.gp.r,asp=1,cex=1.5)
ordiellipse(predict(lda(orig.data,groups))$x, groups,conf=0.95)


#2: CVA can generate differences even when there are none.
data.rand<-matrix(rnorm(150*150),ncol=150)
plot(data.rand,pch=21,bg=col.gp.r,asp=1,cex=1.5)
plot(predict(lda(data.rand[,1:150],groups))$x,pch=21,bg=col.gp.r,asp=1,cex=1.5)

## RDA
plot(Y,pch=shape.gp,bg=col.gp,asp=1,cex=1.5)
pupfish.rda<-rda(Y~SL+species+sex+species:sex)
rda.scores<-predict(pupfish.rda)
plot(rda.scores,pch=shape.gp,bg=col.gp,asp=1,cex=1.5,xlab="RDA 1", ylab="RDA 2")

anova(lm.rrpp(Y~SL*species*sex, data=rdf, print.progress = FALSE))$table
#BUT, species:SL interaction significant (see PC plot). 
  #There are different slopes in data.

##change X matrix, obtain distinct plots
plot(Y,pch=shape.gp,bg=col.gp,asp=1,cex=1.5)
plot(predict(rda(Y~SL*species*sex)),pch=shape.gp,
     bg=col.gp,asp=1,cex=1.5,xlab="RDA 1", ylab="RDA 2")

plot(predict(rda(Y~species*sex)),pch=shape.gp,
     bg=col.gp,asp=1,cex=1.5,xlab="RDA 1", ylab="RDA 2")

