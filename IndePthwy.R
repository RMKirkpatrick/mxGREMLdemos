# Copyright 2019 by Robert M. Kirkpatrick
# Licensed under CC BY 4.0 <http://creativecommons.org/licenses/by/4.0/>

# Demonstration of fitting an independent-pathway model.
# This is still a work-in-progress; it does not converge in a reasonable amount of time.

set.seed(1705160)
library(mvtnorm)
library(Matrix)
library(OpenMx)
options(mxCondenseMatrixSlots=TRUE)
mxOption(NULL,"Default optimizer","NPSOL")
mxOption(NULL,"Analytic Gradients","Yes")

#With more threads, the job will run more quickly, but will require more memory:
mxOption(NULL,"Number of Threads",2)

#The next two options together tell NPSOL to write a log entry for each of its major iterations to a file in   
#the working directory, called 'fort.1' (change the value of "Print file" option to a different positive integer to change the extension 
#to 'fort'):
# mxOption(NULL,"Print level",20)
# mxOption(NULL,"Print file",1)

#The following option will tell NPSOL to check the analytic gradient elements against its numerical derivatives.  If the analytic derivatives
#clearly appear incorrect, NPSOL will terminate, with status code 7.  You need the "Print level" option above to see any details, though.
#Useful if you modify the derivatives of the covariance matrix, but otherwise not worth the added computational effort:
# mxOption(NULL,"Verify level",3)

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

# TODO: Get start values from residual variances & HE regression.

# Begin assembling mxGREML independent-pathway model: #########

#This is the MxExpectationGREML object.  The value of argument 'Xvars' is telling OpenMx to regress all five phenotypes onto both covariates:
xpec <- mxExpectationGREML(V="V",yvars=c("y1","y2","y3","y4"),Xvars=list(c("x1","x2")),blockByPheno=T)

#Custom compute plan to use NPSOL, with analytic gradients but no warm start:
plan <- mxComputeSequence(
	steps=list(
		#mxComputeNewtonRaphson(verbose=5L),
		#mxComputeOnce("fitfunction", "hessian"),
		#mxComputeNelderMead(centerIniSimplex=T,fTolProx=1e-7,xTolProx=1e-7),
		#mxComputeGradientDescent(engine="CSOLNP",useGradient=T,verbose=0L,tolerance=1e-7,maxMajorIter=2L),
		mxComputeGradientDescent(engine="SLSQP",useGradient=T,verbose=5L),
		mxComputeOnce("fitfunction", c("gradient","hessian")),
		mxComputeStandardError(),
		mxComputeHessianQuality(),
		mxComputeReportDeriv(),
		mxComputeReportExpectation()
	))

ipmod <- mxModel(
	"IndePathway",
	xpec,
	plan,
	#sort=FALSE is CRITICALLY IMPORTANT!  It turns off OpenMx's automatic sorting of data rows.
	#We don't want to rearrange the rows, because they are already aligned with the rows and columns of the GRM:
	mxData(observed=widedata,type="raw",sort=FALSE),
	#Start values are such that the variance in each trait is one-quarter common A, common E, unique A, and unique E:
	mxMatrix(
		type="Full",nrow=4,ncol=1,free=c(T,T,T,T),
		values=sqrt(apply(widedata,2,var,na.rm=T)[1:4]/4),
		labels=c("lac1","lac2","lac3","lac4"), name="LA"),
	mxMatrix(
		type="Full",nrow=4,ncol=1,free=c(T,T,T,T),
		values=sqrt(apply(widedata,2,var,na.rm=T)[1:4]/4),
		labels=c("lec1","lec2","lec3","lec4"), name="LE"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(widedata[,"y1"],na.rm=T)/4,labels="vau1",lbound=0,name="Vau1"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(widedata[,"y2"],na.rm=T)/4,labels="vau2",lbound=0,name="Vau2"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(widedata[,"y3"],na.rm=T)/4,labels="vau3",lbound=0,name="Vau3"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(widedata[,"y4"],na.rm=T)/4,labels="vau4",lbound=0,name="Vau4"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(widedata[,"y1"],na.rm=T)/4,labels="veu1",lbound=0.0001,name="Veu1"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(widedata[,"y2"],na.rm=T)/4,labels="veu2",lbound=0.0001,name="Veu2"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(widedata[,"y3"],na.rm=T)/4,labels="veu3",lbound=0.0001,name="Veu3"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(widedata[,"y4"],na.rm=T)/4,labels="veu4",lbound=0.0001,name="Veu4"),
	
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
	
	#Derivatives of V w/r/t free parameters:
	mxAlgebra(rbind(
		cbind(A%x%(2*LA[1,1]),A%x%LA[2,1],A%x%LA[3,1],A%x%LA[4,1]),
		cbind(A%x%LA[2,1],Zip,Zip,Zip),
		cbind(A%x%LA[3,1],Zip,Zip,Zip),
		cbind(A%x%LA[4,1],Zip,Zip,Zip)),name="dV_dlac1"),
	mxAlgebra(rbind(
		cbind(Zip,A%x%LA[1,1],Zip,Zip),
		cbind(A%x%LA[1,1],A%x%(2*LA[2,1]),A%x%LA[3,1],A%x%LA[4,1]),
		cbind(Zip,A%x%LA[3,1],Zip,Zip),
		cbind(Zip,A%x%LA[4,1],Zip,Zip)),name="dV_dlac2"),
	mxAlgebra(rbind(
		cbind(Zip,Zip,A%x%LA[1,1],Zip),
		cbind(Zip,Zip,A%x%LA[2,1],Zip),
		cbind(A%x%LA[1,1],A%x%LA[2,1],A%x%(2*LA[3,1]),A%x%LA[4,1]),
		cbind(Zip,Zip,A%x%LA[4,1],Zip)),name="dV_dlac3"),
	mxAlgebra(rbind(
		cbind(Zip,Zip,Zip,A%x%LA[1,1]),
		cbind(Zip,Zip,Zip,A%x%LA[2,1]),
		cbind(Zip,Zip,Zip,A%x%LA[3,1]),
		cbind(A%x%LA[1,1],A%x%LA[2,1],A%x%LA[3,1],A%x%(2*LA[4,1]))),name="dV_dlac4"),
	
	mxAlgebra(rbind(
		cbind(vec2diag(Uno%x%(2*LE[1,1])),vec2diag(Uno%x%LE[2,1]),vec2diag(Uno%x%LE[3,1]),vec2diag(Uno%x%LE[4,1])),
		cbind(vec2diag(Uno%x%LE[2,1]),Zip,Zip,Zip),
		cbind(vec2diag(Uno%x%LE[3,1]),Zip,Zip,Zip),
		cbind(vec2diag(Uno%x%LE[4,1]),Zip,Zip,Zip)),name="dV_dlec1"),
	mxAlgebra(rbind(
		cbind(Zip,vec2diag(Uno%x%LE[1,1]),Zip,Zip),
		cbind(vec2diag(Uno%x%LE[1,1]),vec2diag(Uno%x%(2*LE[2,1])),vec2diag(Uno%x%LE[3,1]),vec2diag(Uno%x%LE[4,1])),
		cbind(Zip,vec2diag(Uno%x%LE[3,1]),Zip,Zip),
		cbind(Zip,vec2diag(Uno%x%LE[4,1]),Zip,Zip)),name="dV_dlec2"),
	mxAlgebra(rbind(
		cbind(Zip,Zip,vec2diag(Uno%x%LE[1,1]),Zip),
		cbind(Zip,Zip,vec2diag(Uno%x%LE[2,1]),Zip),
		cbind(vec2diag(Uno%x%LE[1,1]),vec2diag(Uno%x%LE[2,1]),vec2diag(Uno%x%(2*LE[3,1])),vec2diag(Uno%x%LE[4,1])),
		cbind(Zip,Zip,vec2diag(Uno%x%LE[4,1]),Zip)),name="dV_dlec3"),
	mxAlgebra(rbind(
		cbind(Zip,Zip,Zip,vec2diag(Uno%x%LE[1,1])),
		cbind(Zip,Zip,Zip,vec2diag(Uno%x%LE[2,1])),
		cbind(Zip,Zip,Zip,vec2diag(Uno%x%LE[3,1])),
		cbind(vec2diag(Uno%x%LE[1,1]),vec2diag(Uno%x%LE[2,1]),vec2diag(Uno%x%LE[3,1]),vec2diag(Uno%x%(2*LE[4,1])))),name="dV_dlec4"),
	
	mxAlgebra(
		rbind(
			cbind(A, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip)), name="dV_dvau1"),
	mxAlgebra(
		rbind(
			cbind(Zip, Zip, Zip, Zip),
			cbind(Zip, A, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip)), name="dV_dvau2"),
	mxAlgebra(
		rbind(
			cbind(Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, A, Zip),
			cbind(Zip, Zip, Zip, Zip)), name="dV_dvau3"),
	mxAlgebra(
		rbind(
			cbind(Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, A)), name="dV_dvau4"),
	
	mxAlgebra(
		rbind(
			cbind(vec2diag(Uno), Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip)), name="dV_dveu1"),
	mxAlgebra(
		rbind(
			cbind(Zip, Zip, Zip, Zip),
			cbind(Zip, vec2diag(Uno), Zip, Zip),
			cbind(Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip)), name="dV_dveu2"),
	mxAlgebra(
		rbind(
			cbind(Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, vec2diag(Uno), Zip),
			cbind(Zip, Zip, Zip, Zip)), name="dV_dveu3"),
	mxAlgebra(
		rbind(
			cbind(Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, vec2diag(Uno))), name="dV_dveu4"),
	
	mxFitFunctionGREML(dV=c(
		lac1="dV_dlac1",lac2="dV_dlac2",lac3="dV_dlac3",lac4="dV_dlac4",
		lec1="dV_dlec1",lec2="dV_dlec2",lec3="dV_dlec3",lec4="dV_dlec4",
		vau1="dV_dvau1",vau2="dV_dvau2",vau3="dV_dvau3",vau4="dV_dvau4",
		veu1="dV_dveu1",veu2="dV_dveu2",veu3="dV_dveu3",veu4="dV_dveu4"))
)
rm(GRM); gc()

#Clobbering the unfitted MxModel with the fitted MxModel object will NOT reduce peak memory demand, but it will allow R to free more memory
#after the call to mxRun() is complete:
ipmod <- mxRun(ipmod)
summary(ipmod,verbose=T)
coef(ipmod); truevals
ipmod$output$fit
ipmod$output$gradient
ipmod$output$hessian
