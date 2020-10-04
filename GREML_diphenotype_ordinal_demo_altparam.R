# Copyright 2019-2020 by Robert M. Kirkpatrick
# Licensed under CC BY 4.0 <http://creativecommons.org/licenses/by/4.0/>

# This script demonstrates a GREML analysis of 2 ordinal phenotypes (3 levels each).
# The model is parameterized in terms of the traits' latent-scale correlations and variance proportions.

library(mvtnorm)
library(Matrix)
library(OpenMx)
options(mxCondenseMatrixSlots=TRUE)
set.seed(171216)

#With more threads, the job will run more quickly, but will require more memory:
mxOption(NULL,"Number of Threads",2)

#Number of simulees (participants):
N <- 1000

#Number of SNPs from which to construct GRM:
msnps <- 50000

#True parameter values, for data generation:
truevals <- c(
	y1_tau1=qnorm(0.5), #<--Threshold #1 for trait 1; 0
	y1_tau2=qnorm(0.75), #<--Threshold #2 for trait 1; 0.6744898
	y2_tau1=qnorm(0.75), #<--Threshold #1 for trait 2; 0.6744898
	y2_tau2=qnorm(0.9), #<--Threshold #2 for trait 2; 1.281552
	y1_h2lat=0.3, #<--Latent-scale heritability of trait 1
	y2_h2lat=0.8, #<--Latent-scale heritability of trait 2
	latgencov=0.5*sqrt(0.3)*sqrt(0.8), #<--Traits' latent-scale genetic covariance; 0.244949
	latenvcov=0.25 #<--Traits' latent-scale nonshared environmental covariance; re = 0.6681531
)

#True parameter values, for observed scale:
truevals2 <- c(
	mu1=(0.5*0)+(0.25*1)+(0.25*2), #<--Observed mean of trait 1; 0.75
	mu2=(0.75*0)+(0.15*1)+(0.1*2), #<--Observed mean of trait 2; 0.35
	vp1=(0.5*0^2)+(0.25*1^2)+(0.25*2^2) - (0.75^2), #<--Observed variance of trait 1; 0.6875
	vp2=(0.75*0^2)+(0.15*1^2)+(0.1*2^2) - (0.35^2), #<--Observe variance of trait 2; 0.4275
	k1=sum(dnorm(c(qnorm(0.5),qnorm(0.75))))^2, #<--Trait 1's heritability conversion factor; 0.5136859
	k2=sum(dnorm(c(qnorm(0.75),qnorm(0.9))))^2, #<--Trait 2's heritability conversion factor; 0.2433201
	#Observed-scale heritability of trait #1, 0.2241539:
	y1_h2obs=0.3 * sum(dnorm(c(qnorm(0.5),qnorm(0.75))))^2 / ((0.5*0^2)+(0.25*1^2)+(0.25*2^2) - (0.75^2)),
	#Observed-scale heritability of trait #2, 0.4553359:
	y2_h2obs=0.8 * sum(dnorm(c(qnorm(0.75),qnorm(0.9))))^2 / ((0.75*0^2)+(0.15*1^2)+(0.1*2^2) - (0.35^2))
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
#Free some memory:
rm(snps, ev); gc()

#Covariance matrix for genetic liabilities:
varglat <- rbind(
	cbind(GRM*truevals["y1_h2lat"], GRM*truevals["latgencov"]),
	cbind(GRM*truevals["latgencov"], GRM*truevals["y2_h2lat"])
)
#Covariance matrix for nonshared environmental liabilities:
varelat <- matrix(data=c(1-truevals["y1_h2lat"], 0.5-truevals["latgencov"], 0.5-truevals["latgencov"], 1-truevals["y2_h2lat"]), nrow=2)
#Genetic liabilities:
glat <- rmvnorm(n=1, mean=rep(0,2*N), sigma=varglat)
#Nonshared environmental liabilities:
elat <- rmvnorm(n=N, mean=c(0,0), sigma=varelat)
#Trait 1, latent scale:
y1lat <- glat[1:N]+elat[,1]
#Trait 2, latent scale:
y2lat <- glat[(N+1):(2*N)]+elat[,2]
#Trait 1, observed scale:
y1obs <- (y1lat>truevals["y1_tau1"]) + (y1lat>truevals["y1_tau2"])
#Trait 2, observed scale:
y2obs <- (y2lat>truevals["y2_tau1"]) + (y2lat>truevals["y2_tau2"])
table(y1obs)
table(y2obs)
#Free some memory:
rm(elat, glat, varelat, varglat); gc()

# Estimate thresholds, using a conventional ordinal-FIML analyisis: ###
threshmod <- mxModel(
	"Threshold_Model",
	mxData(observed=data.frame(y1=mxFactor(y1obs,levels=c(0,1,2)),y2=mxFactor(y2obs,levels=c(0,1,2))),type="raw"),
	mxMatrix(
		type="Full",nrow=2,ncol=2,free=T,
		values=c(
			qnorm(mean(y1obs <= 0, na.rm=T)), qnorm(mean(y1obs <= 1, na.rm=T)), qnorm(mean(y2obs <= 0, na.rm=T)), qnorm(mean(y2obs <= 1, na.rm=T))),
		labels=c("y1_tau1","y1_tau2","y2_tau1","y2_tau2"),
		name="Tau"),
	mxMatrix(type="Stand",nrow=2,free=T,values=cor(y1obs,y2obs,use="pair"),labels="rho",lbound=-0.9999,ubound=0.9999,name="R"),
	mxMatrix(type="Zero",nrow=1,ncol=2,name="Mu"),
	mxFitFunctionML(),
	mxExpectationNormal(covariance="R",means="Mu",dimnames=c("y1","y2"),thresholds="Tau",threshnames=c("y1","y2"))
)
threshmod <- mxRun(threshmod)
summary(threshmod)
#Conversion factor for trait 1:
( k1 <- sum(dnorm(threshmod$Tau$values[,1]))^2 )
#Conversion factor for trait 2:
( k2 <- sum(dnorm(threshmod$Tau$values[,2]))^2 )
#Free some memory:
rm(threshmod); gc()


# Use Haseman-Elston to estimate heritabilities, for use as OpenMx start values ###
#Standardized phenotype 1:
yy1 <- (y1obs-mean(y1obs))/sd(y1obs)
#Standardized phenotype 2:
yy2 <- (y2obs-mean(y2obs))/sd(y2obs)
#Upper triangle of GRM:
U <- matrix(NA_real_,N,N)
U[!lower.tri(U,diag=T)] <- GRM[!lower.tri(GRM,diag=T)]
#Vector of off-diagonal GRM coefficients:
xx <- rep(NA_real_, N*(N-1)/2)
#Vectors to hold cross-products of participants' standardized phenotypes:
yy1.2 <- rep(NA_real_, N*(N-1)/2)
yy2.2 <- rep(NA_real_, N*(N-1)/2)
#Loop to populate vectors to be used in HE regression:
sofar <- 1
for(i in 1:(N-1)){
	yy1.2[sofar:(sofar+N-1-i)] <- yy1[i]*yy1[(i+1):N]
	yy2.2[sofar:(sofar+N-1-i)] <- yy2[i]*yy2[(i+1):N]
	xx[sofar:(sofar+N-1-i)] <- U[i,][!is.na(U[i,])]
	sofar <- sofar+N-i
	#print(i)
}
#Should be no more NAs:
sum(is.na(yy1.2))
sum(is.na(yy2.2))
sum(is.na(xx))
#HE regression for trait 1; slope estimates observed-scale heritability:
( her1 <- lm(yy1.2~xx)$coefficients[2] )
#HE regression for trait 2:
( her2 <- lm(yy2.2~xx)$coefficients[2] )
#Enforce bounds on observed-scale heritability estimates:
her1 <- ifelse(her1<0, 0, her1)
her1 <- ifelse(her1>0.9999, 0.9999, her1)
her2 <- ifelse(her2<0, 0, her2)
her2 <- ifelse(her2>0.9999, 0.9999, her2)
#Free some memory:
rm(U,xx,yy1,yy1.2,yy2,yy2.2); gc()


#With mxGREML, you have to treat phenotypes that are actually ordinal as though they're continuous:
widedata <- cbind(y1=as.numeric(y1obs), y2=as.numeric(y2obs))
#The only fixed effects are the two intercepts:
gremldat <- mxGREMLDataHandler(data=widedata,yvars=c("y1","y2"),blockByPheno=TRUE)
head(gremldat$yX)

#Options for verifying analytic derivatives with NPSOL:
# mxOption(NULL,"Print level",20)
# mxOption(NULL,"Print file",1)
# mxOption(NULL,"Verify level",3)

#Custom compute plan:
plan <- mxComputeSequence(
	steps=list(
		mxComputeGradientDescent(engine="SLSQP",useGradient=T,verbose=5L),
		mxComputeOnce("fitfunction", c("gradient","hessian")),
		mxComputeStandardError(),
		mxComputeHessianQuality(),
		mxComputeReportDeriv(),
		mxComputeReportExpectation()
))
#^^^Note:  If you are running the R GUI under Windows, delete the 'verbose=5L' argument in the above.

#mxGREML model:
gremlmod <- mxModel(
	"Diphenotype_Ordinal",
	mxExpectationGREML(V="V",dataset.is.yX=TRUE),
	mxData(observed=gremldat$yX, type="raw", sort=F),
	plan,
	#Trait 1's latent-scale heritability:
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=her1*var(y1obs,na.rm=T)/k1,labels="h2lat1",lbound=0,ubound=0.9999,name="H2lat1"),
	#Trait 2's latent-scale heritability:
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=her2*var(y2obs,na.rm=T)/k2,labels="h2lat2",lbound=0,ubound=0.9999,name="H2lat2"),
	#Conversion factors:
	mxMatrix(type="Full",nrow=1,ncol=1,free=F,values=k1,name="Konv1"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=F,values=k2,name="Konv2"),
	#Observed-scale phenotypic variances:
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(y1obs,na.rm=T),labels="vp1",lbound=0.0001,name="Vp1"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(y2obs,na.rm=T),labels="vp2",lbound=0.0001,name="Vp2"),
	#Genetic correlation:
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=cor(y1obs,y2obs,use="pair"),labels="ra",lbound=-0.9999,ubound=0.9999,name="Ra"),
	#Nonshared environmental correlation:
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=cor(y1obs,y2obs,use="pair"),labels="re",lbound=-0.9999,ubound=0.9999,name="Re"),
	
	#Fixed matrices for use in MxAlgebras:
	mxMatrix(type="Symm",nrow=N,free=F,values=GRM,name="A"),
	mxMatrix(type="Unit",nrow=N,ncol=1,name="Uno"),
	mxMatrix(type="Zero",nrow=N,ncol=N,name="Zip"),
	
	#Covariance matrix, V:
	mxAlgebra(
		rbind(
			cbind(A%x%(H2lat1*Konv1), Zip),
			cbind(A%x%(Ra*sqrt(H2lat1*Konv1*H2lat2*Konv2)), A%x%(H2lat2*Konv2))
		) +
			rbind(
				cbind(vec2diag(Uno%x%(Vp1-H2lat1*Konv1)), Zip),
				cbind(vec2diag(Uno%x%(Re*sqrt(Vp1-H2lat1*Konv1)*sqrt(Vp2-H2lat2*Konv2))), vec2diag(Uno%x%(Vp2-H2lat2*Konv2)))
			), name="V"),
	#^^^Note that everything above the diagonal in V and its matrix derivatives is ignored by the backend; here, those elements are
	#simply set to zero, to avoid unnecessary computational work.
	
	#The parameterization of this model makes the derivatives of V rather messy...:
	mxAlgebra(
		rbind(
			cbind(A%x%Konv1, Zip),
			cbind(A%x%(Ra*sqrt(Konv1*H2lat2*Konv2)*0.5/sqrt(H2lat1)), Zip)
		) +
			rbind(
				cbind(vec2diag(Uno%x%(-1*Konv1)), Zip),
				cbind(vec2diag(Uno%x%(Re*sqrt(Vp2-H2lat2*Konv2)*0.5/sqrt(Vp1-H2lat1*Konv1)*(-1*Konv1))), Zip)
			), name="dV_dh2lat1"),
	mxAlgebra(
		rbind(
			cbind(Zip, Zip),
			cbind(A%x%(Ra*sqrt(Konv1*H2lat1*Konv2)*0.5/sqrt(H2lat2)), A%x%Konv2)
		) +
			rbind(
				cbind(Zip, Zip),
				cbind(vec2diag(Uno%x%(Re*sqrt(Vp1-H2lat1*Konv1)*0.5/sqrt(Vp2-H2lat2*Konv2)*(-1*Konv2))), vec2diag(Uno%x%(-1*Konv2)))
			), name="dV_dh2lat2"),
	mxAlgebra(
		rbind(
			cbind(vec2diag(Uno), Zip),
			cbind(vec2diag(Uno%x%(Re*sqrt(Vp2-H2lat2*Konv2)*0.5/sqrt(Vp1-H2lat1*Konv1))), Zip)
		), name="dV_dvp1"),
	mxAlgebra(
		rbind(
			cbind(Zip, Zip),
			cbind(vec2diag(Uno%x%(Re*sqrt(Vp1-H2lat1*Konv1)*0.5/sqrt(Vp2-H2lat2*Konv2))), vec2diag(Uno))
		), name="dV_dvp2"),
	mxAlgebra(
		rbind(
			cbind(Zip, Zip),
			cbind(A%x%sqrt(H2lat1*Konv1*H2lat2*Konv2), Zip)
		), name="dV_dra"),
	mxAlgebra(
		rbind(
			cbind(Zip, Zip),
			cbind(vec2diag(Uno%x%(sqrt(Vp1-H2lat1*Konv1)*sqrt(Vp2-H2lat2*Konv2))), Zip)), 
		name="dV_dre"),
	
	mxFitFunctionGREML(dV=c(h2lat1="dV_dh2lat1",h2lat2="dV_dh2lat2",vp1="dV_dvp1",vp2="dV_dvp2",ra="dV_dra",re="dV_dre"))

)
#Free some memory:
rm(GRM,widedata,gremldat,plan); gc()
gremlmod <- mxRun(gremlmod)
gc()
#See results:
summary(gremlmod, verbose=T)
gremlmod$output$fit
gremlmod$output$gradient
gremlmod$output$hessian
