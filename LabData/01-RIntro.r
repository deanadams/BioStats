### An Introduction to R

#R is a statistical programming environment with many built-in mathematical functions
  #and many others that are found in packages that can be installed. Geomorph is one
  #such package.  
#Analyses in R are performed using a series of commands which are written in script files
  #and passed to the main console. The general workflow is:
#1: OPEN R
#2: OPEN R-script
#3: change to working directory 
#4: run analyses

# Some basics in R
3->a  #Assign value to variable
a

b<-c(3,4,5)   # 'c' combines values into vector or list
b
b[2]	#access items in vectors by calling their position

a <- rnorm(50) #generate random normal vector
b <- rnorm(a)
plot(a, b)      # a simple plot

c<-cbind(a,b)    #binds columns together (rbind does same by rows)
c
c[1]   #first element
c[1,]	 #first row
c[,1]  #first column

rbind(a,b)

ls() # See which R objects are now in the R workspace.

### Some base functions
sum(a)
mean(a)
min(a)
max(a)
var(a)
a^2		#square values
sqrt(a)	#NaN for negative values
abs(a)
cor(a,b)

?cor   #call man page for information

rm(list=ls())   #remove items in memory

# Matrix operations
a<-matrix(c(1,0,4,2,-1,1),nrow=3)
b<-matrix(c(1,-1,2,1,1,0),nrow=2)
a
b

c<-t(a)	#matrix transpose
a
c

2*a	#scalar multiplication

#Matrix addition and subtraction
b+c
b-c
a+b		##NOTE: non-conformable matrices (check rxc of your matrices!)

#elementwise multiplication (hadamard product)
c
b
c*b

# matrix multiplication
a%*%b		## %*% is symbol for matrix multiplication
b%*%a		## matrix order matters
	
rm(list=ls())


gl(2,10)	#Generate levels of a factor

## Read data
mydata<-read.csv(file="Data/Lab-01-RIntroData.csv",header=T)
mydata
Y<-as.matrix(mydata[,(2:3)])
FactorA<-as.factor(mydata[,4])
Y
FactorA

# The 'apply' functions (apply, sapply, tapply, etc.) loop over data and do things
apply(Y,2,sd)    #here, we obtain the std for each column of a matrix

tapply(Y[,1],FactorA,mean)	#Obtain means for first column for levels of FactorA
tapply(Y[,2],FactorA,mean)	#Obtain means for first column for levels of FactorA
tapply(Y,FactorA,mean)		#Try entire matrix: doesn't work

rowsum(Y, FactorA)/as.vector(table(FactorA))    #This obtains means.  Could also use a loop.

#### functions
  #In R one can write their own functions

mymean<-function(x){
  n<-length(x)
  tmp<-0
  for (i in 1:n){
    tmp<-tmp+x[i]
  }
  mn<-tmp/n
  return(mn)
}

x<-rnorm(10)
mean(x)  
mymean(x)  #works!

#Some basic statistics
model1<-lm(mydata$y~mydata$x)  #run regression
summary(model1)
anova(model1)  #generates anova table of results

#A plot with regression line
plot(mydata$x,mydata$y)
abline(coef(model1))

model2<-lm(mydata$y~mydata$groups)  #run anova
summary(model2)
anova(model2)

##### Attaching Data
attach(mydata)  #link data table to call columns by name

anova(lm(y~groups))     #ANOVA
anova(lm(y~x))		#regression
anova(lm(y~x*groups))   #ANCOVA

###MULTIVARIATE
ymatrix<-cbind(y,y2)  #generate multivariate Y-matrix
ymatrix
summary(manova(lm(ymatrix~x*groups))) ####NOTE: notation differs slightly

detach(mydata)	#be sure to detach later

rm(list=ls()) 

