####Some Spatial Statistics
   #Packages: spatstat, vegan, spdep, nlme, RRPP

library(spatstat)	
library(RRPP)
library(vegan)
library(spdep)
library(gstat)

# Point patterns

#1: Ripley's K
data(cells)
plot(cells)
plot(Kest(cells))   #shows K for various models: isotropic, poisson, etc.
E<-envelope(cells,Kest,nsim=100,rank=2)
plot(E) # plot shows that observed are underdispersed at small spatial scales


# Generate Spatially Autocorrelated Data  (one could also read some in!)
set.seed(2345)
lat<-runif(50,0,5); long<-runif(50,0,5)
g<-cbind(lat,long)  #create an XY spatial grid
y <- sqrt(diag(g%*%t(g))) + rnorm(nrow(g))  #a spatially-autocorrelated variable
        # Value is associated with distance from origin
t <- sample(rep(c(1,2),length(y)/2))  #2 ecological 'groups' constrained by spatial contingency
sc <- 0.0   #spatial contingency: 0->1
for(i in 1:length(y)){
  z = scale(y,scale=sd(y))
  crit = 1-sc
  crit =c(-crit/2,crit/2)*3
  if(z[i]<=min(crit)) t[i]=1
  if(z[i]>=max(crit)) t[i]=2 
}

plot(g, pch=21, bg=t, cex=y, asp=1, main="Species diversity proprotional to circle size 
Color designates ecological type")

#2: Covariation of Geography and Data
mantel(dist(g), dist(y), permutations = 9999)  # Mantel association
mantel.partial(dist(t), dist(y), dist(g), permutations = 9999)  #3-way Mantel holding group constant

#3: Spatial Autocorrelation
W<-tri2nb(g) #weights with Delauney tesselation
moran.test(y,nb2listw(W))   #positive autocorrelation

#3b: Semivarigram
df <- data.frame(g,y)
coordinates(df) = ~g
res <- variogram(y~1,df)
plot(res, type = "b", main = "Variogram: y") 

#Plot vs. Gaussian model
VarMdl <- vgm(psill=2, model="Gau", nugget=0.1, range=1)
plot(res, model=VarMdl) 

#4: Account for Spatial Non-Independence
ols.fit <- lm(y~t, x=T)  #spatial proximity not considered
summary(ols.fit)
anova(ols.fit)

# GLS models with different spatial autocorrelation structure
  # Warning!  False convergences possible
t <- as.factor(t)
geo=data.frame(g,t,y)

gls.fit.exp = gls(y~t, data=geo, correlation=corExp(form=~lat+long)) #exponential
gls.fit.gaus = gls(y~t, data=geo,  correlation=corGaus(form=~lat+long)) #gaussian
gls.fit.spher = gls(y~t, data=geo,  correlation=corSpher(form=~lat+long)) #spherical
gls.fit.lin = gls(y~t, data=geo,  correlation=corLin(form=~lat+long))

# look at coeffcients: VERY DIFFERENT when spatial non-independence considered
ols.fit
gls.fit.exp 
gls.fit.gaus 
gls.fit.spher 
gls.fit.lin

# model comparisons
AIC(gls.fit.exp, gls.fit.gaus, gls.fit.lin) #Careful when y is multivariate, see Model Selection lecture
   ## Exponential decay of spatial dependence is best model
# Anova
anova(gls.fit.exp)
anova(gls.fit.gaus)
anova(gls.fit.spher)
anova(gls.fit.lin)

#  One could estimate spatial covariance matrix from distances via Gower's approach 
spatCov <- function(x){
  x <- as.matrix(x)
  if(ncol(x) != 2) stop("Need two columns for lat long!")
  d <- dist(x)
  P <- cmdscale(d, 2)
  tcrossprod(P)/(ncol(P) - 1)
}

spat.cov<-spatCov(g)
rdf <- rrpp.data.frame(g=g, y=y,t=t, spat.cov = spat.cov)
res <- lm.rrpp(y~t,Cov = spat.cov, data = rdf)
anova(res)

# NOTE: A more sophisticated spatial covariance matrices is probably preferred, that makes
  #use of the decay of spatial dependent signal.  
# Here, one can extract spatial covariance from nlme and use in RRPP

#exponential matrix
#my.exp <- corExp(1, form=~lat+long)
#my.exp <- Initialize(my.exp,data = geo)
#mycor <- corMatrix(my.exp)

mycor <- corMatrix(Initialize(corExp(1,form = ~lat+long),data=geo))
anova(lm.rrpp(y~t,Cov = mycor, data = geo))

