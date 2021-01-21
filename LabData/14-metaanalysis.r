##############Meta-Analysis
    #Packages: metafor 

library(metafor)
library(RRPP)
mydata<-read.csv("Data/Lab-14-gur_hedge.csv",header=T)
habitat<-mydata[,(1)]
effectd<-as.matrix(mydata[,(2)])                                                             
var<-as.matrix(mydata[,(4)])

#1: Fixed Effects m-a
#no structure
ma.no<-rma.uni(yi=effectd,vi=var,data=mydata,method="FE")
summary(ma.no)
forest(ma.no)  #plot of study effect sizes

#2: Fixed effects categorical m-a
ma.cat<-rma.uni(yi=effectd,vi=var,mods= ~habitat, data=mydata,method="FE")
summary(ma.cat)

#Group means
ma.cat$b[1] #lentic
ma.cat$b[1]+ma.cat$b[2] #marine
ma.cat$b[1]+ma.cat$b[3] #terrestrial

#4: Random Effects m-a
#no structure
rma.uni(yi=effectd,vi=var,data=mydata)

#2: Random effects categorical m-a
ma.catr<-rma.uni(yi=effectd,vi=var,mods= ~habitat, data=mydata)
summary(ma.catr)

#Group means
ma.catr$b[1] #lentic
ma.catr$b[1]+ma.catr$b[2] #marine
ma.catr$b[1]+ma.catr$b[3] #terrestrial

## Some other analyses
#Funnel plot for outliers
funnel(ma.no)

# fail safe number
fsn(yi=effectd,vi=var,data=mydata)
fsn(yi=effectd,vi=var,data=mydata,type="Rosenberg")
#essentially, more than 1500 unpublished non-significant studies need to be published to alter the result

## meta-analysis as GLS
rdf <- rrpp.data.frame(effectd = effectd, habitat=habitat, var.e = diag(as.vector(var)))

anova(lm.rrpp(effectd~1,Cov = rdf$var.e, data = rdf, print.progress = FALSE))
summary(ma.no)  #same SS!

res<- lm.rrpp(effectd~habitat,Cov = rdf$var.e, data = rdf, print.progress = FALSE)
anova(res)
summary(ma.cat) #same Q tests!
res$LM$gls.coefficients #same group effect sizes





