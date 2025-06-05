# Copyright 2025 by Robert M. Kirkpatrick
# Licensed under CC BY 4.0 <http://creativecommons.org/licenses/by/4.0/>

# This script demonstrates how to incorporate a regression ("one-headed path") between two
# endogenous phenotypes in a GREML model.  Specifically, it demonstrates three different ways
# of doing so.

library(mvtnorm)
library(Matrix)
library(OpenMx)
options(mxCondenseMatrixSlots=TRUE)
set.seed(476) 

#With more threads, the job will run more quickly, but will require more memory:
mxOption(NULL,"Number of Threads",2)
#mxOption(NULL,"Nudge zero starts","No")
mxOption(NULL,"Default optimizer","CSOLNP")

#Number of simulees (participants):
N <- 1500

#Number of SNPs from which to construct GRM:
msnps <- 50000

#True parameter values, for data generation:
truevals <- c(
	mu1=10,
	va1=5,
	ve1=5,
	b=1,
	mu2=1,
	va2=0.5,
	ve2=0.5
)

#Generate msnps SNPs.  SNPs are in linkage equilibrium, with MAF in interval [0.05,0.5]:
snps <- matrix(NA_real_,nrow=N,ncol=msnps)
for(mi in 1:msnps){
	maf <- runif(1,min=0.05,max=0.5)
	snps[,mi] <- rbinom(n=N,size=2,prob=maf)
	snps[,mi] <- (snps[,mi]-mean(snps[,mi]))/sd(snps[,mi])
	#print(mi)
}
GRM <- snps%*%t(snps) / msnps #<--#Calculate GRM from SNPs.

#If you don't care whether or not the GRM is positive-definite, you can change the `if(1)` below to `if(0)`.
#The part inside the curly braces "bends" the GRM to the nearest (in a least-squares sense) positive-definite matrix:
if(1){
	ev <- eigen(GRM,symmetric=T,only.values=T) #<--Eigen-decompose the GRM.
	if(!all(ev$values > .Machine$double.eps)){
		GRM <- as.matrix(nearPD(GRM)$mat)
	}
	rm(ev)
}

rm(snps); gc()

# Simulate two phenotypes:
y1 <- as.vector(rmvnorm(n=1,mean=rep(truevals["mu1"],N),sigma=(truevals["va1"]*GRM)+diag(truevals["ve1"],N)))
y2 <- as.vector(rmvnorm(n=1,mean=rep(truevals["mu2"],N),sigma=(truevals["va2"]*GRM)+diag(truevals["ve2"],N))) + truevals["b"]*y1
#^^^Notice the regression slope "b" multiplied by y1; thus, the true model includes a regression of y2 onto y1
widedata <- cbind(y1=y1,y2=y2) #<--Dataset, 1 row per participant.
lmod <- lm(y2~y1) #<--Useful for start values.


# Model #1 demonstrates the simplest (though perhaps the least interesting) way to model a one-headed path 
# between phenotypes--treat y1 as a fixed-effect covariate of y2, and model zero residual covariance
# between y1 & y2:
m1 <- mxModel(
	"EndogenousRegression1",
	mxData(observed=widedata,type="raw",sort=F), #<--N.B. `sort=FALSE`!
	# There are no covariates for y1; y1 is a covariate for y2; we are using ordinary maximum-likelihood:
	mxExpectationGREML(V="V",yvars=c("y1","y2"),Xvars=list(character(0),c("y1")),addOnes=T,REML=F),
	
	# The four variance components:
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(widedata[,"y1"])/2,labels="ve1",lbound=0.0001,name="Ve1"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(lmod$residuals)/2,labels="ve2",lbound=0.0001,name="Ve2"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(widedata[,"y1"])/2,labels="va1",name="Va1"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(lmod$residuals)/2,labels="va2",name="Va2"),
	
	# A column vector of 1s:
	mxMatrix(type="Unit",nrow=N,ncol=1,name="Uno"),
	# The GRM:
	mxMatrix(type="Symm",nrow=N,free=F,values=GRM,name="A"),
	
	# The covariance matrix for phenotype #1:
	mxAlgebra(A%x%Va1 + vec2diag(Uno%x%Ve1), name="V11"),
	# Since y1 is an exogenous covariate for y2, there is no residual covariance between y1 & y2:
	mxMatrix(type="Zero",nrow=N,ncol=N,name="V12"),
	# The covariance matrix for phenotype #2:
	mxAlgebra(A%x%Va2 + vec2diag(Uno%x%Ve2), name="V22"),
	
	# The model-expected covariance matrix, 'V':
	mxAlgebra(rbind(
		cbind(V11,V12),
		cbind(V12,V22)
	),name="V"),

	# The GREML fitfunction:
	mxFitFunctionGREML(infoMatType="expected")
	#^^^Although it's slower, we'll use `infoMatType="expected"`, to make Model #1 more comparable to the models we
	#will fit subsequently in this script.  The average-information-matrix approximation to the Hessian is only
	#available if the model uses an "implicit" model for the phenotypic means (as Model #1 does); the other models in 
	#this demo script use an "explicit" model for the means, and can therefore only use the expected-information
	#approximation.
)
m1 <- mxRun(m1)
if(m1$output$status$code > 1){
	m1$compute <- mxComputeSequence(
		steps=list(
			mxComputeGradientDescent(engine="CSOLNP", verbose=5L),
			mxComputeOnce('fitfunction', c('gradient','hessian')),
			mxComputeStandardError(),
			mxComputeHessianQuality(),
			mxComputeReportDeriv(),
			mxComputeReportExpectation()
		))
	m1 <- mxRun(m1)
}

summary(m1)


# Model #1 had an *implicit* model of the phenotypic mean, i.e. modeled in terms of fixed-effects covariates.
# Model #2, in contrast, has an *explicit* model of the phenotypic mean, that is, it will contain an MxAlgebra
# that represents the phenotypic mean vector, 'yhat'.  It still regresses y1 out of y2, leaving no 
# residual covariance between them.  Model #2 involves putting data into an MxMatrix, which is somewhat unusual in OpenMx:
m2 <- mxModel(
	"EndogenousRegression2",
	mxData(observed=widedata,type="raw",sort=F), #<--N.B. `sort=FALSE`!
	# No `Xvars`, but rather, a `yhat`:
	mxExpectationGREML(V="V",yvars=c("y1","y2"),REML=F,yhat="yhat"),
	
	# The four variance components:
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(widedata[,"y1"])/2,labels="ve1",lbound=0.0001,ubound=NA,name="Ve1"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(lmod$residuals)/2,labels="ve2",lbound=0.0001,ubound=NA,name="Ve2"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(widedata[,"y1"])/2,labels="va1",lbound=NA,ubound=NA,name="Va1"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(lmod$residuals)/2,labels="va2",lbound=NA,ubound=NA,name="Va2"),
	
	# The regression intercept:
	mxMatrix(
		type="Full",nrow=2*N,ncol=1,free=T,
		values=c(rep(mean(widedata[,"y1"]),N), rep(lmod$coef[1],N)),
		labels=c(rep("mu1",N), rep("mu2",N)),lbound=NA,ubound=NA,
		name="M0"),
	# The covariate--zero for y1, and y1 for y2:
	mxMatrix(
		type="Full",nrow=2*N,ncol=1,free=F,
		values=c(rep(0,N), widedata[,"y1"]),
		name="M1"),
	# The regression slope:
	mxMatrix(
		type="Full",nrow=1,ncol=1,free=T,
		values=lmod$coef[2],lbound=NA,ubound=NA,
		labels="b",name="B"),
	
	# A column vector of 1s:
	mxMatrix(type="Unit",nrow=N,ncol=1,name="Uno"),
	# A column vector of zeroes:
	mxMatrix(type="Zero",nrow=N,ncol=1,name="Zip"),
	# The GRM:
	mxMatrix(type="Symm",nrow=N,free=F,values=GRM,name="A"),
	
	# The algebra for 'yhat':
	mxAlgebra(M0 + B*M1, name="yhat"),
	
	# The covariance matrix for phenotype #1:
	mxAlgebra(A%x%Va1 + vec2diag(Uno%x%Ve1), name="V11"),
	# Since y1 is a covariate for y2, there is no residual covariance between y1 & y2:
	mxMatrix(type="Zero",nrow=N,ncol=N,name="V12"),
	# The covariance matrix for phenotype #2:
	mxAlgebra(A%x%Va2 + vec2diag(Uno%x%Ve2), name="V22"),
	
	# The model-expected covariance matrix, 'V':
	mxAlgebra(rbind(
		cbind(V11,V12),
		cbind(V12,V22)
	),name="V"),
	
	# The GREML fitfunction--notice the orthogonality of the mean and covariance parameters in this model:
	mxFitFunctionGREML()
	#^^^It doesn't matter what value we pass for argument `infoMatType`; only the expected-information approximation
	#is available if we are using an explicit means model.
)
m2 <- mxRun(m2)
summary(m2)


# Model #3 also uses an explicit means model, but it models the regression of y2 onto y1 
# in both the mean and covariance:
m3 <- mxModel(
	"EndogenousRegression3",
	mxData(observed=widedata,type="raw",sort=F), #<--N.B. `sort=FALSE`!
	# Again, no `Xvars`, but rather, a `yhat`:
	mxExpectationGREML(V="V",yvars=c("y1","y2"),REML=F,yhat="yhat"),
	
	# The four variance components:
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(widedata[,"y1"])/2,labels="ve1",lbound=0.0001,ubound=NA,name="Ve1"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(lmod$residuals)/2,labels="ve2",lbound=0.0001,ubound=NA,name="Ve2"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(widedata[,"y1"])/2,labels="va1",lbound=NA,ubound=NA,name="Va1"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(lmod$residuals)/2,labels="va2",lbound=NA,ubound=NA,name="Va2"),
	
	# The regression "intercept":
	mxMatrix(
		type="Full",nrow=2*N,ncol=1,free=T,
		values=c(rep(mean(widedata[,"y1"]),N), rep(lmod$coef[1],N)),
		labels=c(rep("mu1",N), rep("mu2",N)),lbound=NA,ubound=NA,
		name="M0"),
	# The "covariate"...
	# Recall that in Model #2, we regressed y2 onto y1.  But here, we instead regress y2 onto the MEAN of y1.
	# That leaves residual covariance between the two phenotypes...
	mxMatrix(
		type="Full",nrow=2*N,ncol=1,
		free=c(rep(F,N), rep(T,N)),
		values=c(rep(0,N), rep(mean(widedata[,"y1"]),N)),
		labels=c(rep(NA,N), rep("mu1",N)),
		name="M1"),
	# The regression slope:
	mxMatrix(
		type="Full",nrow=1,ncol=1,free=T,
		values=lmod$coef[2],lbound=NA,ubound=NA,
		labels="b",name="B"),
	
	# Column vector of 1s:
	mxMatrix(type="Unit",nrow=N,ncol=1,name="Uno"),
	# The GRM:
	mxMatrix(type="Symm",nrow=N,free=F,values=GRM,name="A"),
	
	# The algebra for 'yhat':
	mxAlgebra(M0 + B*M1, name="yhat"),
	
	# The covariance matrix for phenotype #1:
	mxAlgebra(A%x%Va1 + vec2diag(Uno%x%Ve1), name="V11"),
	# Now that we are regressing y2 onto the mean of y1, rather than onto y1 itself,
	# we must model the residual covariance between the phenotypes.  
	# Remember the path-tracing rules...:
	mxAlgebra(V11%x%B, name="V12"),
	# Now that we are regressing y2 onto the mean of y1, rather than onto y1 itself,
	# we must model the total variance of y2 in terms of (1) variance attributable to y1,
	# and (2) variance specific to y2.  Again, remember the path-tracing rules...:
	mxAlgebra(V11%x%(B^2) + A%x%Va2 + vec2diag(Uno%x%Ve2), name="V22"),
	
	# The model-expected covariance matrix, 'V':
	mxAlgebra(rbind(
		cbind(V11,V12),
		cbind(V12,V22)
	),name="V"),

	# The GREML fitfunction:
		mxFitFunctionGREML()
)
m3 <- mxRun(m3)
summary(m3)


# Same -2logL at the solution:
m1$output$fit; m2$output$fit; m3$output$fit
# Point estimates agree:
c(coef(m1),m1$expectation$b)-coef(m2)
coef(m3)-coef(m2)
coef(m3)-c(coef(m1),m1$expectation$b)
cor(c(coef(m1),m1$expectation$b),coef(m2))
cor(c(coef(m1),m1$expectation$b),coef(m3))
cor(coef(m2),coef(m3))
# Running times...:
m1$output$wallTime; m2$output$wallTime; m3$output$wallTime
#^^^Model #1 only has 4, not 7, explicit free parameters to optimize.

# Standard errors agree:
cbind(c(m1$output$standardErrors,sqrt(diag(m1$expectation$bcov))),m2$output$standardErrors,m3$output$standardErrors)
