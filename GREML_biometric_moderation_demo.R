# Copyright 2019-2021 by Robert M. Kirkpatrick
# Licensed under CC BY 4.0 <http://creativecommons.org/licenses/by/4.0/>

#This script demonstrates use of mxGREML to fit a continuous biometric-moderation model
#(a la Purcell, 2002).  It uses the so-called "univariate" form of the model, 
#under which the putative moderator is exogenous.  Construction of the "Gamma"
#and "Lambda" matrices is per Tahmasbi, Evan, Turkheimer, & Keller (2017 preprint),
# http://dx.doi.org/10.1101/191080 .

library(mvtnorm)
library(Matrix)
library(OpenMx)
options(mxCondenseMatrixSlots=TRUE)
set.seed(180114)

#With more threads, the job will run more quickly, but will require more memory:
mxOption(NULL,"Number of Threads",2)
mxOption(NULL,"Default optimizer","NPSOL")
mxOption(NULL,"Nudge zero starts","No")

#Number of simulees (participants):
N <- 2000

#Number of SNPs from which to construct GRM:
msnps <- 50000

#True parameter values:
truevals <- c(
	a0=8.68, #<--Main effect of A
	a1=2.79, #<--A moderation effect
	e0=6.12, #<--Main effect of E
	e1=0.83, #<--E moderation effect
	b=16.11 #<--Main effect of moderator
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
ev <- eigen(GRM,symmetric=T) #<--Eigen-decompose the GRM.

#If you don't care whether or not the GRM is positive-definite, you can comment out this part.  It "bends" the GRM to the nearest 
#(in a least-squares sense) positive-definite matrix:
if(!all(ev$values > .Machine$double.eps)){
	GRM <- as.matrix(nearPD(GRM)$mat)
}
rm(snps, ev); gc()

#Simulate data:
m <- runif(n=N)
Lambda <- outer(m,m)
LambdaGRM <- Lambda * GRM
Gamma <- matrix(NA_real_,N,N)
for(i in 1:N){
	Gamma[i,] <- m[i]
}
Gamma <- Gamma + t(Gamma)
GammaGRM <- Gamma*GRM
K1 <- diag(2*m)
K2 <- diag(m^2)
V <- (truevals["a0"]^2)*GRM + truevals["a0"]*truevals["a1"]*LambdaGRM + (truevals["a1"]^2)*GammaGRM + diag(rep(truevals["e0"]^2,N)) + 
	truevals["e0"]*truevals["e1"]*K1 + (truevals["e1"]^2)*K2
y <- as.vector(rmvnorm(n=1, mean=100+truevals["b"]*m, sigma=V))
rm(V); gc()

#Haseman-Elston for start values:
yc <- lm(y~m)$residuals
#Upper triangle of GRM:
U1 <- matrix(NA_real_,N,N)
U1[!lower.tri(U1,diag=T)] <- GRM[!lower.tri(GRM,diag=T)]
hex <- rep(NA_real_, N*(N-1)/2)
hey <- rep(NA_real_, N*(N-1)/2)
sofar <- 1
for(i in 1:(N-1)){
	hey[sofar:(sofar+N-1-i)] <- yc[i]*yc[(i+1):N]
	hex[sofar:(sofar+N-1-i)] <- U1[i,][!is.na(U1[i,])]
	sofar <- sofar+N-i
	#print(i)
}
rm(U1); gc()
her <- lm(hey~hex)$coefficients[2]


#Options to set if you are going to use NPSOL to verify analytic derivatives:
# mxOption(NULL,"Print level",20)
# mxOption(NULL,"Print file",2)
# mxOption(NULL,"Verify level",3)

#Compute plan; uses Newton-Raphson; fast, but not very robust:
plan <- mxComputeSequence(
	steps=list(
		mxComputeNewtonRaphson(verbose=5L),
		#mxComputeGradientDescent(engine="NPSOL",verbose=5L),
		mxComputeOnce("fitfunction", c("gradient","hessian")),
		mxComputeStandardError(),
		mxComputeHessianQuality(),
		mxComputeReportDeriv(),
		mxComputeReportExpectation()
	))
#^^^Note:  If you are running the R GUI under Windows, delete the 'verbose=5L' argument in the above.

#GREML model:
gremlmod <- mxModel(
	"GREML",
	plan,
	mxData(cbind(y=y,m=m),type="raw",sort=FALSE),
	mxExpectationGREML(V="V",yvars="y",Xvars="m"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=sqrt(her),labels="a0",name="A0"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0,labels="a1",name="A1"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=sqrt(max(var(yc,na.rm=T)-her,0)),labels="e0",name="E0"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0,labels="e1",name="E1"),
	
	mxMatrix(type="Symm",nrow=N,free=F,values=GRM,name="A"),
	mxMatrix(type="Symm",nrow=N,free=F,values=LambdaGRM,name="LambdaA"),
	mxMatrix(type="Symm",nrow=N,free=F,values=GammaGRM,name="GammaA"),
	mxMatrix(type="Full",nrow=N,ncol=1,free=F,values=diag(K1),name="K1d"),
	mxMatrix(type="Full",nrow=N,ncol=1,free=F,values=diag(K2),name="K2d"),
	mxMatrix(type="Unit",nrow=N,ncol=1,name="Uno"),
	
	mxAlgebra( (A0^2)%x%A + (A0*A1)%x%LambdaA + (A1^2)%x%GammaA + vec2diag(Uno%x%E0^2) + vec2diag(K1d%x%(E0*E1)) + vec2diag(K2d%x%(E1^2)),
						 name="V" ),
	
	mxAlgebra( (2*A0)%x%A + A1%x%LambdaA, name="dV_da0"),
	mxAlgebra( A0%x%LambdaA + (2*A1)%x%GammaA, name="dV_da1"),
	mxAlgebra( vec2diag(Uno%x%(2*E0)) + vec2diag(K1d%x%E1), name="dV_de0"),
	mxAlgebra( vec2diag(K1d%x%E0) + vec2diag(K2d%x%(2*E1)), name="dV_de1"),
	
	mxFitFunctionGREML(dV=c(a0="dV_da0", a1="dV_da1", e0="dV_de0", e1="dV_de1"))
)
rm(Gamma,GammaGRM,GRM,K1,K2,Lambda,LambdaGRM,plan); gc()
gremlmod <- mxRun(gremlmod)
gc()
object.size(gremlmod) #<--How much memory does the fitted MxModel take up?:

#If Newton-Raphson doesn't reach a good solution, try again with (possibly warm-started) NPSOL :
if( !(gremlmod$output$status$code %in% c(0,1)) ){
	ws <- try(chol(gremlmod$output$hessian))
	if("try-error" %in% class(ws)){ws <- NULL}
	gremlmod$compute <- mxComputeSequence(steps=list(
		mxComputeGradientDescent(engine="NPSOL",useGradient=T,verbose=5L,warmStart=ws),
		mxComputeOnce('fitfunction', c('gradient','hessian')),
		mxComputeStandardError(),
		mxComputeHessianQuality(),
		mxComputeReportDeriv(),
		mxComputeReportExpectation()
	))
	#^^^Note:  If you are running the R GUI under Windows, delete the 'verbose=5L' argument in the above.
	gremlmod <- mxRun(gremlmod)
	gc()
}

#If NPSOL doesn't reach a good solution, try again with SLSQP:
if( !(gremlmod$output$status$code %in% c(0,1)) ){
	gremlmod$compute <- mxComputeSequence(steps=list(
		mxComputeGradientDescent(engine="SLSQP",verbose=5L),
		mxComputeOnce('fitfunction', c('gradient','hessian')),
		mxComputeStandardError(),
		mxComputeHessianQuality(),
		mxComputeReportDeriv(),
		mxComputeReportExpectation()
	))
	#^^^Note:  If you are running the R GUI under Windows, delete the 'verbose=5L' argument in the above.
	gremlmod <- mxRun(gremlmod)
	gc()
}

summary(gremlmod,verbose=T)
