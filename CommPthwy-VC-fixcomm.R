# Copyright 2019-2021 by Robert M. Kirkpatrick
# Licensed under CC BY 4.0 <http://creativecommons.org/licenses/by/4.0/>

# Model identification: common-factor's variance is fixed to 1.0 with an implicit constraint.
# Parameterization: direct-variance.

set.seed(170516)
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
	vac=0.5, #<--Common additive-genetic variance
	l1=sqrt(2.0), #<--Factor loading for phenotype #1
	l2=sqrt(2.0), #<--factor loading for phenotype #2
	l3=sqrt(2.0), #<--factor loading for phenotype #3
	l4=sqrt(2.0), #<--factor loading for phenotype #4
	l5=sqrt(2.0), #<--factor loading for phenotype #5
	va1=1.0, #<--Additive-genetic uniqueness of phenotype #1
	va2=1.0, #<--Additive-genetic uniqueness of phenotype #2
	va3=1.0, #<--Additive-genetic uniqueness of phenotype #3
	va4=1.0, #<--Additive-genetic uniqueness of phenotype #4
	va5=1.0, #<--Additive-genetic uniqueness of phenotype #5
	ve1=1.0, #<--Nonshared environmental uniqueness of phenotype #1
	ve2=1.0, #<--Nonshared environmental uniqueness of phenotype #2
	ve3=1.0, #<--Nonshared environmental uniqueness of phenotype #3
	ve4=1.0, #<--Nonshared environmental uniqueness of phenotype #4
	ve5=1.0 #<--Nonshared environmental uniqueness of phenotype #5
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
Ac <- as.vector(z2mvnorm(rnorm(N),ev,truevals["vac"]))
A1 <- as.vector(z2mvnorm(rnorm(N),ev,truevals["va1"]))
A2 <- as.vector(z2mvnorm(rnorm(N),ev,truevals["va2"]))
A3 <- as.vector(z2mvnorm(rnorm(N),ev,truevals["va3"]))
A4 <- as.vector(z2mvnorm(rnorm(N),ev,truevals["va4"]))
A5 <- as.vector(z2mvnorm(rnorm(N),ev,truevals["va5"]))
rm(ev); gc()
Ec <- rnorm(N,sd=sqrt(1-truevals["vac"]))
E1 <- rnorm(N,sd=sqrt(truevals["ve1"]))
E2 <- rnorm(N,sd=sqrt(truevals["ve2"]))
E3 <- rnorm(N,sd=sqrt(truevals["ve3"]))
E4 <- rnorm(N,sd=sqrt(truevals["ve4"]))
E5 <- rnorm(N,sd=sqrt(truevals["ve5"]))

#Generate scores on covariates (in this script they are independent of the phenotypes):
X1 <- runif(N)
X2 <- runif(N)

#Now generate phenotype scores:
Y1 <- truevals["l1"]*(Ac+Ec) + A1 + E1
Y2 <- truevals["l2"]*(Ac+Ec) + A2 + E2
Y3 <- truevals["l3"]*(Ac+Ec) + A3 + E3
Y4 <- truevals["l4"]*(Ac+Ec) + A4 + E4
Y5 <- truevals["l5"]*(Ac+Ec) + A5 + E5

#Assemble data matrix:
widedata <- cbind(Y1,Y2,Y3,Y4,Y5,X1,X2)
colnames(widedata) <- c("y1","y2","y3","y4","y5","x1","x2")
rm(Y1,Y2,Y3,Y4,Y5,X1,X2,Ac,Ec,A1,A2,A3,A4,A5,E1,E2,E3,E4,E5); gc()

# done generating data ###############

#Fit a non-biometric phenotypic factor model first.  It's a good way to get starting values for the loadings in the mxGREML model, and can alert 
#you to possible data-handling errors:  
phenofacmodel <- mxModel(
	name="PhenFactor",
	mxData(observed=widedata,type="raw"),
	#This is a conventional FIML analysis, so the regression coefficients are explicit free parameters:
	mxMatrix(type="Full",nrow=1,ncol=5,free=T,values=colMeans(widedata[,1:5],na.rm=T),
					 labels=c("b01","b02","b03","b04","b05"),name="B0"),
	mxMatrix(type="Full",nrow=1,ncol=5,free=T,values=0,labels=c("b11","b12","b13","b14","b15"),name="B1"),
	mxMatrix(type="Full",nrow=1,ncol=5,free=T,values=0,labels=c("b21","b22","b23","b24","b25"),name="B2"),
	mxMatrix(type="Full",nrow=1,ncol=2,free=F,labels=c("data.x1","data.x2"),name="Xdef"),
	mxAlgebra(B0 + (Xdef[1,1]%x%B1) + (Xdef[1,2]%x%B2), name="Mu",dimnames=list(NULL,c("y1","y2","y3","y4","y5"))),
	mxMatrix(type="Full",nrow=5,ncol=1,free=c(T,T,T,T,T),values=1,
					 labels=c("l1","l2","l3","l4","l5"),
					 name="Lambda",
					 dimnames=list(c("y1","y2","y3","y4","y5"),"F1")),
	#The common factor's variance, fixed to 1.0 to identify the scale:
	mxMatrix(type="Symm",nrow=1,free=F,values=1,labels="fvar",lbound=0,name="Fvar"),
	#The unique variances:
	mxMatrix(
		type="Diag",nrow=5,free=T,labels=c("u1","u2","u3","u4","u5"),lbound=0.0001,name="U",
		values=diag(cov(widedata[,1:5],use="pair"))),
	mxAlgebra(Lambda%*%Fvar%*%t(Lambda),name="Comm"),
	mxAlgebra(Comm+U,name="Sigmay",dimnames=list(c("y1","y2","y3","y4","y5"),c("y1","y2","y3","y4","y5"))),
	mxExpectationNormal(covariance="Sigmay",means="Mu"),
	mxFitFunctionML()
)
phenofacrun <- mxRun(phenofacmodel)
summary(phenofacrun)
iniest <- coef(phenofacrun) #<--#Store point estimates to use as start values in the mxGREML analysis

# Begin assembling mxGREML common-pathway model: #########

#This is the MxExpectationGREML object.  The value of argument 'Xvars' is telling OpenMx to regress all five phenotypes onto both covariates:
xpec <- mxExpectationGREML(V="V",yvars=c("y1","y2","y3","y4","y5"),Xvars=list(c("x1","x2")))
#To drop, say, x2 from the model:
# xpec <- mxExpectationGREML(V="V",yvars=c("y1","y2","y3","y4","y5"),Xvars=list("x1"))
#If you don't want to include ANY covariates, you'd instead do:
# xpec <- mxExpectationGREML(V="V",yvars=c("y1","y2","y3","y4","y5"))
#As an example, the following would regress y1 and y2 onto x1 and x2, regress y3 onto x1, regress y4 onto x2, and include no covariates for y5:
# xpec <- mxExpectationGREML(V="V",yvars=c("y1","y2","y3","y4","y5"),Xvars=list(c("x1","x2"),c("x1","x2"),"x1","x2",character(0)))

#Custom compute plan to use NPSOL, with analytic gradients but no warm start:
plan <- mxComputeSequence(
	steps=list(
		mxComputeGradientDescent(engine="NPSOL",verbose=5L),
		#^^^Note:  If you are running the R GUI under Windows, delete the 'verbose=5L' argument in the above.
		mxComputeOnce("fitfunction", c("gradient","hessian")),
		mxComputeStandardError(),
		mxComputeHessianQuality(),
		mxComputeReportDeriv(),
		mxComputeReportExpectation()
	))
#^^^On the IBG laptops at the 2020 Boulder Workshop, I ran this script with NPSOL and with
#the number of threads set to 4.  NPSOL appeared to enter a closed loop in which it repeatedly
#got a non-finite fitfunction value.  If that's the case for you, use SLSQP instead.

cpmod <- mxModel(
	"CommonPathway",
	xpec,
	plan,
	#sort=FALSE is CRITICALLY IMPORTANT!  It turns off OpenMx's automatic sorting of data rows.
	#We don't want to rearrange the rows, because they are already aligned with the rows and columns of the GRM:
	mxData(observed=widedata,type="raw",sort=FALSE),
	#vac is started at one-half the total common-factor variance:
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=1/2,labels="vac",
					 lbound=0, ubound=1,
					 name="Vac"),
	#The factor loadings are started at their estimates from the non-biometric factor model:
	mxMatrix(type="Full",nrow=5,ncol=1,free=c(T,T,T,T,T),values=c(iniest["l1"],iniest["l2"],iniest["l3"],iniest["l4"],iniest["l5"]),
					 labels=c("l1","l2","l3","l4","l5"), name="Lambda"),
	#The unique variance components are started each at one-half the corresponding values from the non-biometric model:
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=iniest["u1"]/2,labels="va1",lbound=0,name="Va1"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=iniest["u2"]/2,labels="va2",lbound=0,name="Va2"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=iniest["u3"]/2,labels="va3",lbound=0,name="Va3"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=iniest["u4"]/2,labels="va4",lbound=0,name="Va4"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=iniest["u5"]/2,labels="va5",lbound=0,name="Va5"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=iniest["u1"]/2,labels="ve1",lbound=0.0001,name="Ve1"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=iniest["u2"]/2,labels="ve2",lbound=0.0001,name="Ve2"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=iniest["u3"]/2,labels="ve3",lbound=0.0001,name="Ve3"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=iniest["u4"]/2,labels="ve4",lbound=0.0001,name="Ve4"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=iniest["u5"]/2,labels="ve5",lbound=0.0001,name="Ve5"),
	mxMatrix(type="Symm",nrow=N,free=F,values=GRM,name="A"),
	mxMatrix(type="Unit",nrow=N,ncol=1,name="Uno"),
	mxMatrix(type="Zero",nrow=N,ncol=N,name="Zip"),
	
	mxAlgebra( (A%x%Vac) + vec2diag(Uno%x%(1-Vac)), name="SigmaFac"), #<--Covariance matrix for all N individuals' common factors.
	mxAlgebra( Lambda%*%t(Lambda), name="LamLamT"), #<--Outer product of loadings vector with itself.
	#Heritable uniqueness matrix:
	mxAlgebra(
		rbind(
			cbind(A%x%Va1, Zip, Zip, Zip, Zip),
			cbind(Zip, A%x%Va2, Zip, Zip, Zip),
			cbind(Zip, Zip, A%x%Va3, Zip, Zip),
			cbind(Zip, Zip, Zip, A%x%Va4, Zip),
			cbind(Zip, Zip, Zip, Zip, A%x%Va5)), name="UniqueA"),
	#Nonshared environmental uniqueness matrix:
	mxAlgebra(
		rbind(
			cbind(vec2diag(Uno%x%Ve1), Zip, Zip, Zip, Zip),
			cbind(Zip, vec2diag(Uno%x%Ve2), Zip, Zip, Zip),
			cbind(Zip, Zip, vec2diag(Uno%x%Ve3), Zip, Zip),
			cbind(Zip, Zip, Zip, vec2diag(Uno%x%Ve4), Zip),
			cbind(Zip, Zip, Zip, Zip, vec2diag(Uno%x%Ve5))), name="UniqueE"),
	mxAlgebra( (LamLamT%x%SigmaFac) + UniqueA + UniqueE, name="V"), #<--Covariance matrix
	
	#Derivatives of V w/r/t free parameters:
	mxAlgebra(LamLamT%x%(A-vec2diag(Uno)), name="dV_dvac"),
	mxAlgebra(rbind(
		cbind(SigmaFac%x%(2*Lambda[1,1]),SigmaFac%x%Lambda[2,1],SigmaFac%x%Lambda[3,1],SigmaFac%x%Lambda[4,1],SigmaFac%x%Lambda[5,1]),
		cbind(SigmaFac%x%Lambda[2,1],Zip,Zip,Zip,Zip),
		cbind(SigmaFac%x%Lambda[3,1],Zip,Zip,Zip,Zip),
		cbind(SigmaFac%x%Lambda[4,1],Zip,Zip,Zip,Zip),
		cbind(SigmaFac%x%Lambda[5,1],Zip,Zip,Zip,Zip)),
		name="dV_dl1"),
	mxAlgebra(rbind(
		cbind(Zip,SigmaFac%x%Lambda[1,1],Zip,Zip,Zip),
		cbind(SigmaFac%x%Lambda[1,1],SigmaFac%x%(2*Lambda[2,1]),SigmaFac%x%Lambda[3,1],SigmaFac%x%Lambda[4,1],SigmaFac%x%Lambda[5,1]),
		cbind(Zip,SigmaFac%x%Lambda[3,1],Zip,Zip,Zip),
		cbind(Zip,SigmaFac%x%Lambda[4,1],Zip,Zip,Zip),
		cbind(Zip,SigmaFac%x%Lambda[5,1],Zip,Zip,Zip)),
		name="dV_dl2"),
	mxAlgebra(rbind(
		cbind(Zip,Zip,SigmaFac%x%Lambda[1,1],Zip,Zip),
		cbind(Zip,Zip,SigmaFac%x%Lambda[2,1],Zip,Zip),
		cbind(SigmaFac%x%Lambda[1,1],SigmaFac%x%Lambda[2,1],SigmaFac%x%(2*Lambda[3,1]),SigmaFac%x%Lambda[4,1],SigmaFac%x%Lambda[5,1]),
		cbind(Zip,Zip,SigmaFac%x%Lambda[4,1],Zip,Zip),
		cbind(Zip,Zip,SigmaFac%x%Lambda[5,1],Zip,Zip)),
		name="dV_dl3"),
	mxAlgebra(rbind(
		cbind(Zip,Zip,Zip,SigmaFac%x%Lambda[1,1],Zip),
		cbind(Zip,Zip,Zip,SigmaFac%x%Lambda[2,1],Zip),
		cbind(Zip,Zip,Zip,SigmaFac%x%Lambda[3,1],Zip),
		cbind(SigmaFac%x%Lambda[1,1],SigmaFac%x%Lambda[2,1],SigmaFac%x%Lambda[3,1],SigmaFac%x%(2*Lambda[4,1]),SigmaFac%x%Lambda[5,1]),
		cbind(Zip,Zip,Zip,SigmaFac%x%Lambda[5,1],Zip)),
		name="dV_dl4"),
	mxAlgebra(rbind(
		cbind(Zip,Zip,Zip,Zip,SigmaFac%x%Lambda[1,1]),
		cbind(Zip,Zip,Zip,Zip,SigmaFac%x%Lambda[2,1]),
		cbind(Zip,Zip,Zip,Zip,SigmaFac%x%Lambda[3,1]),
		cbind(Zip,Zip,Zip,Zip,SigmaFac%x%Lambda[4,1]),
		cbind(SigmaFac%x%Lambda[1,1],SigmaFac%x%Lambda[2,1],SigmaFac%x%Lambda[3,1],SigmaFac%x%Lambda[4,1],SigmaFac%x%(2*Lambda[5,1]))),
		name="dV_dl5"),
	mxAlgebra(
		rbind(
			cbind(A, Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip, Zip)), name="dV_dva1"),
	mxAlgebra(
		rbind(
			cbind(Zip, Zip, Zip, Zip, Zip),
			cbind(Zip, A, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip, Zip)), name="dV_dva2"),
	mxAlgebra(
		rbind(
			cbind(Zip, Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, A, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip, Zip)), name="dV_dva3"),
	mxAlgebra(
		rbind(
			cbind(Zip, Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, A, Zip),
			cbind(Zip, Zip, Zip, Zip, Zip)), name="dV_dva4"),
	mxAlgebra(
		rbind(
			cbind(Zip, Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip, A)), name="dV_dva5"),
	mxAlgebra(
		rbind(
			cbind(vec2diag(Uno), Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip, Zip)), name="dV_dve1"),
	mxAlgebra(
		rbind(
			cbind(Zip, Zip, Zip, Zip, Zip),
			cbind(Zip, vec2diag(Uno), Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip, Zip)), name="dV_dve2"),
	mxAlgebra(
		rbind(
			cbind(Zip, Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, vec2diag(Uno), Zip, Zip),
			cbind(Zip, Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip, Zip)), name="dV_dve3"),
	mxAlgebra(
		rbind(
			cbind(Zip, Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, vec2diag(Uno), Zip),
			cbind(Zip, Zip, Zip, Zip, Zip)), name="dV_dve4"),
	mxAlgebra(
		rbind(
			cbind(Zip, Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip, Zip),
			cbind(Zip, Zip, Zip, Zip, vec2diag(Uno))), name="dV_dve5"),
	#GREML fitfunction object:
	mxFitFunctionGREML(
		dV=c(vac="dV_dvac",l1="dV_dl1",l2="dV_dl2",l3="dV_dl3",l4="dV_dl4",l5="dV_dl5",
				 va1="dV_dva1",va2="dV_dva2",va3="dV_dva3",va4="dV_dva4",va5="dV_dva5",
				 ve1="dV_dve1",ve2="dV_dve2",ve3="dV_dve3",ve4="dV_dve4",ve5="dV_dve5"))
)
rm(GRM); gc()

#Clobbering the unfitted MxModel with the fitted MxModel object will NOT reduce peak memory demand, but it will allow R to free more memory
#after the call to mxRun() is complete:
cpmod <- mxRun(cpmod)
object.size(cpmod) #<--How much memory does the fitted MxModel take up?:
summary(cpmod,verbose=T)
coef(cpmod); truevals
cpmod$output$fit
cpmod$output$gradient
cpmod$output$hessian
