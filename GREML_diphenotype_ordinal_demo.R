# Copyright 2019-2021 by Robert M. Kirkpatrick
# Licensed under CC BY 4.0 <http://creativecommons.org/licenses/by/4.0/>

# This script demonstrates a GREML analysis of 2 ordinal phenotypes (3 levels each).
# The model is parameterized in terms of the traits' observed-scale variance-covariance components.

library(mvtnorm)
library(Matrix)
library(OpenMx)
options(mxCondenseMatrixSlots=TRUE)
library(polycor)
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
#Nonshared-environmental liabilities:
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

#Use Haseman-Elston to estimate heritabilities:
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
# mxOption(NULL,"Function precision",1e-7)

gremlmod <- mxModel(
	"Diphenotype_Ordinal",
	mxExpectationGREML(V="V",dataset.is.yX=TRUE),
	#Custom compute plan:
	mxComputeSequence(
		steps=list(
			#mxComputeGradientDescent(engine="NPSOL",verbose=5L),
			mxComputeNewtonRaphson(verbose=5L),
			mxComputeOnce("fitfunction", c("gradient","hessian")),
			mxComputeConfidenceInterval(
				plan=mxComputeGradientDescent(engine="SLSQP"),
				constraintType="ineq", verbose=5L
					),
			mxComputeStandardError(),
			mxComputeHessianQuality(),
			mxComputeReportDeriv(),
			mxComputeReportExpectation()
		)),
	#^^^Note:  If you are running the R GUI under Windows, delete the 'verbose=5L' arguments in the above.
	mxData(observed=gremldat$yX, type="raw", sort=F),
	
	#Trait 1's observed-scale additive-genetic variance:
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=her1*var(y1obs),labels="va1",name="Va1"),
	#Trait 2's observed-scale additive-genetic variance:
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=her2*var(y2obs),labels="va2",name="Va2"),
	#Trait 1's observed-scale nonshared-environmental variance:
	mxMatrix(type="Diag",nrow=N,free=T,values=(1-her1)*var(y1obs),labels="ve1",lbound=0,name="Ve1"),
	#Trait 2's observed-scale nonshared-environmental variance:
	mxMatrix(type="Diag",nrow=N,free=T,values=(1-her2)*var(y2obs),labels="ve2",lbound=0,name="Ve2"),
	#Observed-scale genetic covariance:
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.5*cov(y1obs,y2obs),labels="ca",name="Ca"),
	#Observed-scale nonshared-environmental covariance:
	mxMatrix(type="Diag",nrow=N,free=T,values=0.5*cov(y1obs,y2obs),labels="ce",name="Ce"),
	
	#Observed-to-latent conversion factors:
	mxMatrix(type="Full",nrow=1,ncol=1,free=F,values=k1,name="Konv1"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=F,values=k2,name="Konv2"),
	#Observed-scale phenotypic variances:
	mxAlgebra(Va1+Ve1[1,1], name="Vp1"),
	mxAlgebra(Va2+Ve2[1,1], name="Vp2"),
	#Observed-scale heritabilities:
	mxAlgebra(Va1 / Vp1, name="y1_h2obs"),
	mxAlgebra(Va2 / Vp2, name="y2_h2obs"),
	#Latent-scale heritabilities:
	mxAlgebra(Va1 / Konv1, name="y1_h2lat"),
	mxAlgebra(Va2 / Konv2, name="y2_h2lat"),
	#Genetic correlation:
	mxAlgebra(Ca / sqrt(Va1) / sqrt(Va2), name="ra"),
	
	mxMatrix(type="Symm",nrow=N,free=F,values=GRM,name="A"),
	mxMatrix(type="Iden",nrow=N,name="I"),
	mxMatrix(type="Zero",nrow=N,ncol=N,name="Zip"),
	
	#Model-expected covariance matrix, V:
	mxAlgebra(
		rbind(
			cbind(A%x%Va1, Zip),
			cbind(A%x%Ca, A%x%Va2)
		) +
			rbind(
				cbind(Ve1, Zip),
				cbind(Ce, Ve2)
			), name="V"),
	#^^^Note that everything above the diagonal in V is ignored by the backend; here, those elements are
	#simply set to zero, to avoid unnecessary computational work.
	
	#Derivatives of V w/r/t free parameters:
	mxMatrix(
		type="Symm",nrow=2*N,free=F,name="dV_dva1",
		values=rbind(
			cbind(GRM, matrix(0,N,N)),
			cbind(matrix(0,N,N), matrix(0,N,N))
		)),
	mxMatrix(
		type="Symm",nrow=2*N,free=F,name="dV_dva2",
		values=rbind(
			cbind(matrix(0,N,N), matrix(0,N,N)),
			cbind(matrix(0,N,N), GRM)
		)),
	mxMatrix(
		type="Full",nrow=2*N,ncol=2*N,free=F,name="dV_dca",
		values=rbind(
			cbind(matrix(0,N,N), matrix(0,N,N)),
			cbind(GRM, matrix(0,N,N))
		)),
	mxMatrix(
		type="Symm",nrow=2*N,free=F,name="dV_dve1",
		values=rbind(
			cbind(diag(N), matrix(0,N,N)),
			cbind(matrix(0,N,N), matrix(0,N,N))
		)),
	mxMatrix(
		type="Symm",nrow=2*N,free=F,name="dV_dve2",
		values=rbind(
			cbind(matrix(0,N,N), matrix(0,N,N)),
			cbind(matrix(0,N,N), diag(N))
		)),
	mxMatrix(
		type="Full",nrow=2*N,ncol=2*N,free=F,name="dV_dce",
		values=rbind(
			cbind(matrix(0,N,N), matrix(0,N,N)),
			cbind(diag(N), matrix(0,N,N))
		)),
	
	#GREML fitfunction:
	mxFitFunctionGREML(dV=c(va1="dV_dva1",va2="dV_dva2",ca="dV_dca",ve1="dV_dve1",ve2="dV_dve2",ce="dV_dce")),
	#Let's request CIs for the latent-scale parameters:
	mxCI(c("y1_h2lat","y2_h2lat","ra"))
)
rm(gremldat, GRM, widedata); gc()
gremlmod <- mxRun(gremlmod)
gc()
summary(gremlmod, verbose=T)
mxEval(y1_h2obs, gremlmod)
mxEval(y2_h2obs, gremlmod)
mxEval(y1_h2lat, gremlmod)
mxEval(y2_h2lat, gremlmod)
mxEval(ra, gremlmod)
