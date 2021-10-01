# Copyright 2019-2021 by Robert M. Kirkpatrick
# Licensed under CC BY 4.0 <http://creativecommons.org/licenses/by/4.0/>

# Demonstration of fitting an independent-pathway model.
# Parameterization: loadings onto the IPs are free, IP's variances fixed to 1.0.

set.seed(1705160)
library(mvtnorm)
library(Matrix)
library(OpenMx)
options(mxCondenseMatrixSlots=TRUE)
mxOption(NULL,"Default optimizer","SLSQP")

#With more threads, the job will run more quickly, but will require more memory:
mxOption(NULL,"Number of Threads",2)

#Number of simulees (participants):
N <- 500

#Number of SNPs from which to construct GRM:
msnps <- 50000

#True parameter values:
truevals <- c(
	lac1=sqrt(2.0),
	lac2=sqrt(2.0),
	lac3=sqrt(2.0),
	lac4=sqrt(2.0),
	lec1=sqrt(2.0),
	lec2=sqrt(2.0),
	lec3=sqrt(2.0),
	lec4=sqrt(2.0),
	vau1=1.0,
	vau2=1.0,
	vau3=1.0,
	vau4=1.0,
	veu1=1.0,
	veu2=1.0,
	veu3=1.0,
	veu4=1.0
)

#Function to scale and rotate independent normal deviates:
z2mvnorm <- function(z,ev,const){
	R <- t(ev$vectors %*% (t(ev$vectors) * sqrt(ev$values*const)))
	retval <- matrix(z,nrow=1) %*% R
	return(retval)
}

#Function to obtain start values for heritabilities from Haseman-Elston regression:
HEstartvals <- function(yy, GRM, bound=TRUE){
	#Standardized phenotype:
	yy1 <- (yy-mean(yy))/sd(yy)
	#Sample size, inferred from order of GRM:
	N <- nrow(GRM)
	#Upper triangle of GRM:
	U <- matrix(NA_real_,N,N)
	U[!lower.tri(U,diag=T)] <- GRM[!lower.tri(GRM,diag=T)]
	#Vector of off-diagonal GRM coefficients:
	xx <- rep(NA_real_, N*(N-1)/2)
	#Vectors to hold cross-products of participants' standardized phenotypes:
	yy1.2 <- rep(NA_real_, N*(N-1)/2)
	#Loop to populate vectors to be used in HE regression:
	sofar <- 1
	for(i in 1:(N-1)){
		yy1.2[sofar:(sofar+N-1-i)] <- yy1[i]*yy1[(i+1):N]
		xx[sofar:(sofar+N-1-i)] <- U[i,][!is.na(U[i,])]
		sofar <- sofar+N-i
		#print(i)
	}
	#Should be no more NAs:
	if(is.na(sum(is.na(yy1.2)))){stop("NAs detected")}
	if(is.na(sum(is.na(xx)))){stop("NAs detected")}
	#HE regression; slope estimates observed-scale heritability:
	her <- lm(yy1.2~xx)$coefficients[2]
	if(bound){
		#Enforce bounds on observed-scale heritability estimates:
		her <- ifelse(her<0, 0, her)
		her <- ifelse(her>0.9999, 0.9999, her)
	}
	return(her)
}

#Generate msnps SNPs.  SNPs are in linkage equilibrium, with MAF in interval [0.05,0.5]:
snps <- matrix(NA_real_,nrow=N,ncol=msnps)
for(mi in 1:msnps){
	maf <- runif(1,min=0.05,max=0.5)
	snps[,mi] <- rbinom(n=N,size=2,prob=maf)
	snps[,mi] <- (snps[,mi]-mean(snps[,mi]))/sd(snps[,mi])
	#print(mi)
}
GRM <- snps%*%t(snps) / msnps #<--#Calculate GRM from SNPs.
ev <- eigen(GRM,symmetric=T) #<--Eigen-decompose the GRM.

#If you don't care whether or not the GRM is positive-definite, you can comment out this part.  It "bends" the GRM to the nearest 
#(in a least-squares sense) positive-definite matrix, and then eigen-decomposes it again:
if(!all(ev$values > .Machine$double.eps)){
	GRM <- as.matrix(nearPD(GRM)$mat)
	ev <- eigen(GRM,symmetric=T)
}

rm(snps); gc() 

# Generate phenotypic data: ##########
#First generate common and unique factor scores:
Ac <- as.vector(z2mvnorm(rnorm(N),ev,1.0))
A1 <- as.vector(z2mvnorm(rnorm(N),ev,truevals["vau1"]))
A2 <- as.vector(z2mvnorm(rnorm(N),ev,truevals["vau2"]))
A3 <- as.vector(z2mvnorm(rnorm(N),ev,truevals["vau3"]))
A4 <- as.vector(z2mvnorm(rnorm(N),ev,truevals["vau4"]))
rm(ev); gc()
Ec <- rnorm(N)
E1 <- rnorm(N,sd=sqrt(truevals["veu1"]))
E2 <- rnorm(N,sd=sqrt(truevals["veu2"]))
E3 <- rnorm(N,sd=sqrt(truevals["veu3"]))
E4 <- rnorm(N,sd=sqrt(truevals["veu4"]))

#Generate scores on covariates (in this script they are independent of the phenotypes):
X1 <- runif(N)
X2 <- runif(N)

#Now generate phenotype scores:
Y1 <- truevals["lac1"]*Ac + truevals["lec1"]*Ec + A1 + E1
Y2 <- truevals["lac2"]*Ac + truevals["lec2"]*Ec + A2 + E2
Y3 <- truevals["lac3"]*Ac + truevals["lec3"]*Ec + A3 + E3
Y4 <- truevals["lac4"]*Ac + truevals["lec4"]*Ec + A4 + E4

#Assemble data matrix:
widedata <- cbind(Y1,Y2,Y3,Y4,X1,X2)
colnames(widedata) <- c("y1","y2","y3","y4","x1","x2")
rm(Y1,Y2,Y3,Y4,X1,X2,Ac,Ec,A1,A2,A3,A4,E1,E2,E3,E4); gc()

#Get start values from residual variances & HE regression:
resvarcs <- rep(NA_real_,4)
h2starts <- rep(NA_real_,4)
for(i in 1:4){
	resvarcs[i] <- var(lm(widedata[,i]~widedata[,5]+widedata[,6])$residuals)
	h2starts[i] <- HEstartvals(widedata[,i],GRM,T)
}

# Begin assembling mxGREML independent-pathway model: #########

#This is the MxExpectationGREML object.  The value of argument 'Xvars' is telling OpenMx to regress all five phenotypes onto both covariates:
xpec <- mxExpectationGREML(V="V",yvars=c("y1","y2","y3","y4"),Xvars=list(c("x1","x2")),blockByPheno=T)

ipmod <- mxModel(
	"IndePathway",
	xpec,
	#sort=FALSE is CRITICALLY IMPORTANT!  It turns off OpenMx's automatic sorting of data rows.
	#We don't want to rearrange the rows, because they are already aligned with the rows and columns of the GRM:
	mxData(observed=widedata,type="raw",sort=FALSE),
	#Start values are such that the variance in each trait is one-quarter common A, common E, unique A, and unique E:
	mxMatrix(
		type="Full",nrow=4,ncol=1,free=c(T,T,T,T),
		values=sqrt(h2starts*resvarcs/2),
		labels=c("lac1","lac2","lac3","lac4"), name="LA"),
	mxMatrix(
		type="Full",nrow=4,ncol=1,free=c(T,T,T,T),
		values=sqrt((resvarcs - h2starts*resvarcs)/2),
		labels=c("lec1","lec2","lec3","lec4"), name="LE"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=h2starts[1]*resvarcs[1]/2,labels="vau1",lbound=0,name="Vau1"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=h2starts[2]*resvarcs[2]/2,labels="vau2",lbound=0,name="Vau2"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=h2starts[3]*resvarcs[3]/2,labels="vau3",lbound=0,name="Vau3"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=h2starts[4]*resvarcs[4]/2,labels="vau4",lbound=0,name="Vau4"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=(resvarcs[1] - h2starts[1]*resvarcs[1])/2,labels="veu1",lbound=0.0001,name="Veu1"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=(resvarcs[2] - h2starts[2]*resvarcs[2])/2,labels="veu2",lbound=0.0001,name="Veu2"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=(resvarcs[3] - h2starts[3]*resvarcs[3])/2,labels="veu3",lbound=0.0001,name="Veu3"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=(resvarcs[4] - h2starts[4]*resvarcs[4])/2,labels="veu4",lbound=0.0001,name="Veu4"),
	
	#Some matrices constant w/r/t free parameters:
	mxMatrix(type="Symm",nrow=N,free=F,values=GRM,name="A"),
	mxMatrix(type="Unit",nrow=N,ncol=1,name="Uno"),
	mxMatrix(type="Zero",nrow=N,ncol=N,name="Zip"),
	
	#Outer product of loadings onto common A & E:
	mxAlgebra(LA%*%t(LA),name="LALAT"),
	mxAlgebra(LE%*%t(LE),name="LELET"),
	#Heritable uniqueness matrix:
	mxAlgebra(
		rbind(
			cbind(A%x%Vau1, Zip, Zip, Zip),
			cbind(Zip, A%x%Vau2, Zip, Zip),
			cbind(Zip, Zip, A%x%Vau3, Zip),
			cbind(Zip, Zip, Zip, A%x%Vau4)), name="UniqueA"),
	#Nonshared environmental uniqueness matrix:
	mxAlgebra(
		rbind(
			cbind(vec2diag(Uno%x%Veu1), Zip, Zip, Zip),
			cbind(Zip, vec2diag(Uno%x%Veu2), Zip, Zip),
			cbind(Zip, Zip, vec2diag(Uno%x%Veu3), Zip),
			cbind(Zip, Zip, Zip, vec2diag(Uno%x%Veu4))),name="UniqueE"),
	mxAlgebra( (LALAT%x%A) + (LELET%x%vec2diag(Uno)) + UniqueA + UniqueE, name="V" ), #<--Covariance matrix
	
	mxFitFunctionGREML()
)
rm(GRM); gc()

#Clobbering the unfitted MxModel with the fitted MxModel object will NOT reduce peak memory demand, but it will allow R to free more memory
#after the call to mxRun() is complete:
ipmod <- mxRun(ipmod)
object.size(ipmod) #<--How much memory does the fitted MxModel take up?:
summary(ipmod,verbose=T)
coef(ipmod); truevals
ipmod$output$fit
ipmod$output$gradient
ipmod$output$hessian
