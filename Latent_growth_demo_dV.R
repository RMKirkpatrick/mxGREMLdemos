# Copyright 2019-2025 by Robert M. Kirkpatrick
# Licensed under CC BY 4.0 <http://creativecommons.org/licenses/by/4.0/>

#mxGREML script to fit a latent-growth model.
#Parameterization: direct variance components

library(mvtnorm)
library(Matrix)
library(OpenMx)
library(lme4)
options(mxCondenseMatrixSlots=TRUE)
set.seed(210930)

#With more threads, the job will run more quickly, but will require more memory:
mxOption(NULL,"Number of Threads",2)

#Number of simulees (participants):
N <- 500

#Number of SNPs from which to construct GRM:
msnps <- 50000

#True parameter values:
truevals <- c(
	val=0.5, #<--Additive-genetic variance in intercept
	vel=0.5, #<--Nonshared-environmental variance in intercept
	vas=0.5, #<--Additive-genetic variance in slope
	ves=0.5, #<--Nonshared-environmental variance in slope
	ca=-0.15, #<--Additive-genetic covariance between intercept and slope
	ce=-0.12, #<--Nonshared-environmental covariance between intercept and slope
	vau=1.0, #<--Additive-genetic uniqueness
	veu=1.0, #<--Nonshared-environmental uniqueness
	ALmean=3.25, #<--Mean heritable latent intercept
	ASmean=0.25, #<--Mean heritable latent slope
	ELmean=3.25, #<--Mean environmental latent intercept
	ESmean=0.25 #<--Mean environmental latent slope
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

rm(snps,ev); gc()

#Age (the time metric), centered at expected wave-3 age:
age1 <- round(9 + runif(N,-0.25,0.25) - 13, 2)
age2 <- round(11 + runif(N,-0.25,0.25) - 13, 2)
age3 <- round(13 + runif(N,-0.25,0.25) - 13, 2)
age4 <- round(15 + runif(N,-0.25,0.25) - 13, 2)
age5 <- round(17 + runif(N,-0.25,0.25) - 13, 2)

#Additive-genetic covariance matrix for common factors:
varFA <- rbind(
	cbind(GRM*truevals["val"],GRM*truevals["ca"]),
	cbind(GRM*truevals["ca"],GRM*truevals["vas"])
)

#Nonshared environmental covariance matrix for common factors:
varFE <- rbind(
	cbind(diag(truevals["vel"],N),diag(truevals["ce"],N)),
	cbind(diag(truevals["ce"],N),diag(truevals["ves"],N))
)
#Heritable latent intercepts & slopes:
ALAS <- rmvnorm(1, mean=c(rep(truevals["ALmean"],N),rep(truevals["ASmean"],N)), sigma=varFA)
#Environmental latent intercepts and slopes:
ELES <- rmvnorm(1, mean=c(rep(truevals["ELmean"],N),rep(truevals["ESmean"],N)), sigma=varFE)

#Phenotypes:
y1 <- ALAS[1:N] + ELES[1:N] + (ALAS[(N+1):(2*N)]+ELES[(N+1):(2*N)])*age1 + rmvnorm(1, mean=rep(0,N), sigma=GRM*truevals["vau"]) + 
	rnorm(N,mean=0,sd=sqrt(truevals["veu"]))
y2 <- ALAS[1:N] + ELES[1:N] + (ALAS[(N+1):(2*N)]+ELES[(N+1):(2*N)])*age2 + rmvnorm(1, mean=rep(0,N), sigma=GRM*truevals["vau"]) + 
	rnorm(N,mean=0,sd=sqrt(truevals["veu"]))
y3 <- ALAS[1:N] + ELES[1:N] + (ALAS[(N+1):(2*N)]+ELES[(N+1):(2*N)])*age3 + rmvnorm(1, mean=rep(0,N), sigma=GRM*truevals["vau"]) + 
	rnorm(N,mean=0,sd=sqrt(truevals["veu"]))
y4 <- ALAS[1:N] + ELES[1:N] + (ALAS[(N+1):(2*N)]+ELES[(N+1):(2*N)])*age4 + rmvnorm(1, mean=rep(0,N), sigma=GRM*truevals["vau"]) + 
	rnorm(N,mean=0,sd=sqrt(truevals["veu"]))
y5 <- ALAS[1:N] + ELES[1:N] + (ALAS[(N+1):(2*N)]+ELES[(N+1):(2*N)])*age5 + rmvnorm(1, mean=rep(0,N), sigma=GRM*truevals["vau"]) + 
	rnorm(N,mean=0,sd=sqrt(truevals["veu"]))
# Finished generating data ###


#Non-biometric REML analysis, to get start values:
y.all <- c(y1,y2,y3,y4,y5)
age.all <- c(age1,age2,age3,age4,age5)
id <- as.factor(rep(1:N,5))
remlmod <- lmer(y.all~1+age.all + (1+age.all|id))
( remlsumm <- summary(remlmod))
truevals

#Create wide data object:
widedata <- cbind(
	y1=as.vector(y1), age1=age1, y2=as.vector(y2), age2=age2, y3=as.vector(y3), age3=age3, y4=as.vector(y4), age4=age4, y5=as.vector(y5), age5=age5)
head(widedata)

#Remove unneeded objects:
rm(y.all,age.all,ALAS,ELES,varFA,varFE,y1,y2,y3,y4,y5)
gc()

#GREML expectation:
xpec <- mxExpectationGREML(
	V="V", yvars=c("y1","y2","y3","y4","y5"), 
	Xvars=list("age1","age2","age3","age4","age5"), 
	staggerZeroes=FALSE)

# #Options for verifying analytic derivatives with NPSOL:
# mxOption(NULL,"Print level",20)
# mxOption(NULL,"Print file",1)
# mxOption(NULL,"Verify level",3)
# mxOption(NULL,"Function precision",1e-7)

#Custom compute plan:
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

gremlmod <- mxModel(
	"LatentGrowth",
	xpec, #<--Expectation
	plan, #<--Compute plan
	#sort=FALSE is CRITICALLY IMPORTANT!  It turns off OpenMx's automatic sorting of data rows.
	#We don't want to rearrange the rows, because they are already aligned with the rows and columns of the GRM:
	mxData(observed=widedata,type="raw",sort=FALSE),
	
	#Matrices to hold free parameters relevant to the latent intercept and slope; variance components from non-biometric REML analysis are used as
	#start values:
	mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=0.5*VarCorr(remlmod)[[1]][1,1], labels="val", name="Val"),
	#Note the column vector:
	mxMatrix(type="Full", nrow=N, ncol=1, free=T, values=0.5*VarCorr(remlmod)[[1]][1,1], labels="vel", name="Vel"),
	mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=0.5*VarCorr(remlmod)[[1]][2,2], labels="vas", name="Vas"),
	mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=0.5*VarCorr(remlmod)[[1]][2,2], labels="ves", name="Ves"),
	mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=0.5*VarCorr(remlmod)[[1]][1,2], labels="ca", name="Ca"),
	mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=0.5*VarCorr(remlmod)[[1]][1,2], labels="ce", name="Ce"),
	
	#Matrices to hold unique-variance components; each is initialized at one-half the non-biometric REML residual variance:
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.5*attr(VarCorr(remlmod),"sc"),labels="vau",name="Vau"),
	#Note the diagonal matrix:
	mxMatrix(type="Diag",nrow=N,free=T,values=0.5*attr(VarCorr(remlmod),"sc"),labels="veu",lbound=0.0001,name="Veu"),
	
	#GRM:
	mxMatrix(type="Symm",nrow=N,ncol=N,free=F,values=GRM,name="A"),
	
	#Column vectors of ages:
	mxMatrix(type="Full",nrow=N,ncol=1,free=F,values=age1,name="Age1"),
	mxMatrix(type="Full",nrow=N,ncol=1,free=F,values=age2,name="Age2"),
	mxMatrix(type="Full",nrow=N,ncol=1,free=F,values=age3,name="Age3"),
	mxMatrix(type="Full",nrow=N,ncol=1,free=F,values=age4,name="Age4"),
	mxMatrix(type="Full",nrow=N,ncol=1,free=F,values=age5,name="Age5"),
	
	#Identity matrix:
	mxMatrix(type="Iden",nrow=N,name="I"),
	#Matrix of zeroes:
	mxMatrix(type="Zero",nrow=N,ncol=N,name="Zip"),
	
	#Algebras to be reused multiple times:
	mxAlgebra(A%x%Val, name="Aval"),
	mxAlgebra(Aval + vec2diag(Vel), name="AvalPlusVel"),
	
	#Blocks (submatrices) of the derivative of V w/r/t additive-genetic variance in slope:
	mxMatrix(type="Full", nrow=N, ncol=N, free=F, values=outer(age1,age1)*GRM, name="dV_dvas_11"),
	mxMatrix(type="Full", nrow=N, ncol=N, free=F, values=outer(age2,age1)*GRM, name="dV_dvas_21"),
	mxMatrix(type="Full", nrow=N, ncol=N, free=F, values=outer(age2,age2)*GRM, name="dV_dvas_22"),
	mxMatrix(type="Full", nrow=N, ncol=N, free=F, values=outer(age3,age1)*GRM, name="dV_dvas_31"),
	mxMatrix(type="Full", nrow=N, ncol=N, free=F, values=outer(age3,age2)*GRM, name="dV_dvas_32"),
	mxMatrix(type="Full", nrow=N, ncol=N, free=F, values=outer(age3,age3)*GRM, name="dV_dvas_33"),
	mxMatrix(type="Full", nrow=N, ncol=N, free=F, values=outer(age4,age1)*GRM, name="dV_dvas_41"),
	mxMatrix(type="Full", nrow=N, ncol=N, free=F, values=outer(age4,age2)*GRM, name="dV_dvas_42"),
	mxMatrix(type="Full", nrow=N, ncol=N, free=F, values=outer(age4,age3)*GRM, name="dV_dvas_43"),
	mxMatrix(type="Full", nrow=N, ncol=N, free=F, values=outer(age4,age4)*GRM, name="dV_dvas_44"),
	mxMatrix(type="Full", nrow=N, ncol=N, free=F, values=outer(age5,age1)*GRM, name="dV_dvas_51"),
	mxMatrix(type="Full", nrow=N, ncol=N, free=F, values=outer(age5,age2)*GRM, name="dV_dvas_52"),
	mxMatrix(type="Full", nrow=N, ncol=N, free=F, values=outer(age5,age3)*GRM, name="dV_dvas_53"),
	mxMatrix(type="Full", nrow=N, ncol=N, free=F, values=outer(age5,age4)*GRM, name="dV_dvas_54"),
	mxMatrix(type="Full", nrow=N, ncol=N, free=F, values=outer(age5,age5)*GRM, name="dV_dvas_55"),
	
	#Blocks (submatrices) of the derivative of V w/r/t additive-genetic covariance between intercept and slope:
	mxMatrix(type="Full", nrow=N, ncol=N, free=F, values=diag(age1)%*%GRM + GRM%*%diag(age1), name="dV_dca_11"),
	mxMatrix(type="Full", nrow=N, ncol=N, free=F, values=diag(age2)%*%GRM + GRM%*%diag(age1), name="dV_dca_21"),
	mxMatrix(type="Full", nrow=N, ncol=N, free=F, values=diag(age2)%*%GRM + GRM%*%diag(age2), name="dV_dca_22"),
	mxMatrix(type="Full", nrow=N, ncol=N, free=F, values=diag(age3)%*%GRM + GRM%*%diag(age1), name="dV_dca_31"),
	mxMatrix(type="Full", nrow=N, ncol=N, free=F, values=diag(age3)%*%GRM + GRM%*%diag(age2), name="dV_dca_32"),
	mxMatrix(type="Full", nrow=N, ncol=N, free=F, values=diag(age3)%*%GRM + GRM%*%diag(age3), name="dV_dca_33"),
	mxMatrix(type="Full", nrow=N, ncol=N, free=F, values=diag(age4)%*%GRM + GRM%*%diag(age1), name="dV_dca_41"),
	mxMatrix(type="Full", nrow=N, ncol=N, free=F, values=diag(age4)%*%GRM + GRM%*%diag(age2), name="dV_dca_42"),
	mxMatrix(type="Full", nrow=N, ncol=N, free=F, values=diag(age4)%*%GRM + GRM%*%diag(age3), name="dV_dca_43"),
	mxMatrix(type="Full", nrow=N, ncol=N, free=F, values=diag(age4)%*%GRM + GRM%*%diag(age4), name="dV_dca_44"),
	mxMatrix(type="Full", nrow=N, ncol=N, free=F, values=diag(age5)%*%GRM + GRM%*%diag(age1), name="dV_dca_51"),
	mxMatrix(type="Full", nrow=N, ncol=N, free=F, values=diag(age5)%*%GRM + GRM%*%diag(age2), name="dV_dca_52"),
	mxMatrix(type="Full", nrow=N, ncol=N, free=F, values=diag(age5)%*%GRM + GRM%*%diag(age3), name="dV_dca_53"),
	mxMatrix(type="Full", nrow=N, ncol=N, free=F, values=diag(age5)%*%GRM + GRM%*%diag(age4), name="dV_dca_54"),
	mxMatrix(type="Full", nrow=N, ncol=N, free=F, values=diag(age5)%*%GRM + GRM%*%diag(age5), name="dV_dca_55"),
	
	#Model-expected covariance matrix, V, specified blockwise.  Notice that blocks entirely above the diagonal can be left as zero;
	#the mxGREML backend ignores the values of elements above the diagonal of V and its derivatives:
	mxAlgebra(rbind(
		cbind(
			#Block 1,1:
			AvalPlusVel + dV_dvas_11%x%vas + vec2diag((Age1*Age1)%x%Ves) + vec2diag((Age1+Age1)%x%(Ce)) + dV_dca_11%x%Ca + 
				A%x%Vau + Veu,
			#Block 1,2 thru 1,5:
			Zip,Zip,Zip,Zip),
		cbind(
			#Block 2,1:
			AvalPlusVel + dV_dvas_21%x%vas + vec2diag((Age2*Age1)%x%Ves) + vec2diag((Age2+Age1)%x%(Ce)) + dV_dca_21%x%Ca,
			#Block 2,2:
			AvalPlusVel + dV_dvas_22%x%vas + vec2diag((Age2*Age2)%x%Ves) + vec2diag((Age2+Age2)%x%(Ce)) + dV_dca_22%x%Ca + 
				A%x%Vau + Veu,
			#Block 2,3 thru 2,5:
			Zip, Zip, Zip),
		cbind(
			#Block 3,1:
			AvalPlusVel + dV_dvas_31%x%vas + vec2diag((Age3*Age1)%x%Ves) + vec2diag((Age3+Age1)%x%(Ce)) + dV_dca_31%x%Ca,
			#Block 3,2:
			AvalPlusVel + dV_dvas_32%x%vas + vec2diag((Age3*Age2)%x%Ves) + vec2diag((Age3+Age2)%x%(Ce)) + dV_dca_32%x%Ca,
			#Block 3,3:
			AvalPlusVel + dV_dvas_33%x%vas + vec2diag((Age3*Age3)%x%Ves) + vec2diag((Age3+Age3)%x%(Ce)) + dV_dca_33%x%Ca + 
				A%x%Vau + Veu,
			#Block 3,4 & 3,5:
			Zip, Zip),
		cbind(
			#Block 4,1:
			AvalPlusVel + dV_dvas_41%x%vas + vec2diag((Age4*Age1)%x%Ves) + vec2diag((Age4+Age1)%x%(Ce)) + dV_dca_41%x%Ca,
			#Block 4,2:
			AvalPlusVel + dV_dvas_42%x%vas + vec2diag((Age4*Age2)%x%Ves) + vec2diag((Age4+Age2)%x%(Ce)) + dV_dca_42%x%Ca,
			#Block 4,3:
			AvalPlusVel + dV_dvas_43%x%vas + vec2diag((Age4*Age3)%x%Ves) + vec2diag((Age4+Age3)%x%(Ce)) + dV_dca_43%x%Ca,
			#Block 4,4:
			AvalPlusVel + dV_dvas_44%x%vas + vec2diag((Age4*Age4)%x%Ves) + vec2diag((Age4+Age4)%x%(Ce)) + dV_dca_44%x%Ca + 
				A%x%Vau + Veu,
			#Block 4,5:
			Zip),
		cbind(
			#Block 5,1:
			AvalPlusVel + dV_dvas_51%x%vas + vec2diag((Age5*Age1)%x%Ves) + vec2diag((Age5+Age1)%x%(Ce)) + dV_dca_51%x%Ca,
			#Block 5,2:
			AvalPlusVel + dV_dvas_52%x%vas + vec2diag((Age5*Age2)%x%Ves) + vec2diag((Age5+Age2)%x%(Ce)) + dV_dca_52%x%Ca,
			#Block 5,3:
			AvalPlusVel + dV_dvas_53%x%vas + vec2diag((Age5*Age3)%x%Ves) + vec2diag((Age5+Age3)%x%(Ce)) + dV_dca_53%x%Ca,
			#Block 5,4:
			AvalPlusVel + dV_dvas_54%x%vas + vec2diag((Age5*Age4)%x%Ves) + vec2diag((Age5+Age4)%x%(Ce)) + dV_dca_54%x%Ca,
			#Block 5,5:
			AvalPlusVel + dV_dvas_55%x%vas + vec2diag((Age5*Age5)%x%Ves) + vec2diag((Age5+Age5)%x%(Ce)) + dV_dca_55%x%Ca + 
				A%x%Vau + Veu)
	), name="V"),
	
	#Derivative of V w/r/t additive-genetic variance in intercept:
	mxAlgebra(
		rbind(
			cbind(A,Zip,Zip,Zip,Zip),
			cbind(A,A,Zip,Zip,Zip),
			cbind(A,A,A,Zip,Zip),
			cbind(A,A,A,A,Zip),
			cbind(A,A,A,A,A)
		), name="dV_dval"),
	
	#Derivative of V w/r/t nonshared-environmental variance in intercept:
	mxAlgebra(
		rbind(
			cbind(I,Zip,Zip,Zip,Zip),
			cbind(I,I,Zip,Zip,Zip),
			cbind(I,I,I,Zip,Zip),
			cbind(I,I,I,I,Zip),
			cbind(I,I,I,I,I)
		), name="dV_dvel"),
	
	#Derivative of V w/r/t additive-genetic variance in slope:
	mxAlgebra(
		rbind(
			cbind(dV_dvas_11, Zip, Zip, Zip, Zip),
			cbind(dV_dvas_21, dV_dvas_22, Zip, Zip, Zip),
			cbind(dV_dvas_31, dV_dvas_32, dV_dvas_33, Zip, Zip),
			cbind(dV_dvas_41, dV_dvas_42, dV_dvas_43, dV_dvas_44, Zip),
			cbind(dV_dvas_51, dV_dvas_52, dV_dvas_53, dV_dvas_54, dV_dvas_55)
		), name="dV_dvas"),
	
	#Derivative of V w/r/t additive-genetic covariance between slope and intercept:
	mxAlgebra(
		rbind(
			cbind(dV_dca_11, Zip, Zip, Zip, Zip),
			cbind(dV_dca_21, dV_dca_22, Zip, Zip, Zip),
			cbind(dV_dca_31, dV_dca_32, dV_dca_33, Zip, Zip),
			cbind(dV_dca_41, dV_dca_42, dV_dca_43, dV_dca_44, Zip),
			cbind(dV_dca_51, dV_dca_52, dV_dca_53, dV_dca_54, dV_dca_55)
		), name="dV_dca"),
	
	#Derivative of V w/r/t nonshared-environmental variance in slope:
	mxAlgebra(rbind(
		cbind(vec2diag(Age1*Age1), Zip, Zip, Zip, Zip),
		cbind(vec2diag(Age2*Age1), vec2diag(Age2*Age2), Zip, Zip, Zip),
		cbind(vec2diag(Age3*Age1), vec2diag(Age3*Age2), vec2diag(Age3*Age3), Zip, Zip),
		cbind(vec2diag(Age4*Age1), vec2diag(Age4*Age2), vec2diag(Age4*Age3), vec2diag(Age4*Age4), Zip),
		cbind(vec2diag(Age5*Age1), vec2diag(Age5*Age2), vec2diag(Age5*Age3), vec2diag(Age5*Age4), vec2diag(Age5*Age5))
	), name="dV_dves"),
	
	#Derivative of V w/r/t nonshared-environmental covariance between slope and intercept:
	mxAlgebra(rbind(
		cbind(vec2diag(Age1+Age1), Zip, Zip, Zip, Zip),
		cbind(vec2diag(Age2+Age1), vec2diag(Age2+Age2), Zip, Zip, Zip),
		cbind(vec2diag(Age3+Age1), vec2diag(Age3+Age2), vec2diag(Age3+Age3), Zip, Zip),
		cbind(vec2diag(Age4+Age1), vec2diag(Age4+Age2), vec2diag(Age4+Age3), vec2diag(Age4+Age4), Zip),
		cbind(vec2diag(Age5+Age1), vec2diag(Age5+Age2), vec2diag(Age5+Age3), vec2diag(Age5+Age4), vec2diag(Age5+Age5))
	), name="dV_dce"),
	
	#Derivative of V w/r/t additive-genetic unique variance:
	mxAlgebra(
		rbind(
			cbind(A,Zip,Zip,Zip,Zip),
			cbind(Zip,A,Zip,Zip,Zip),
			cbind(Zip,Zip,A,Zip,Zip),
			cbind(Zip,Zip,Zip,A,Zip),
			cbind(Zip,Zip,Zip,Zip,A)
		), name="dV_dvau"),
	
	#Derivative of V w/r/t nonshared-environmental unique variance:
	mxAlgebra(
		rbind(
			cbind(I,Zip,Zip,Zip,Zip),
			cbind(Zip,I,Zip,Zip,Zip),
			cbind(Zip,Zip,I,Zip,Zip),
			cbind(Zip,Zip,Zip,I,Zip),
			cbind(Zip,Zip,Zip,Zip,I)
		), name="dV_dveu"),
	
	#GREML fitfunction object:
	mxFitFunctionGREML(dV=c(
		val="dV_dval",vel="dV_dvel",vas="dV_dvas",ves="dV_dves",ca="dV_dca",ce="dV_dce",vau="dV_dvau",veu="dV_dveu"))
)

#Remove unneeded objects:
rm(widedata,GRM); gc()

#Clobbering the unfitted MxModel with the fitted MxModel object will not reduce peak memory demand, but it will allow R to free more memory
#after the call to mxRun() is complete:
gremlmod <- mxRun(gremlmod)
gc()
object.size(gremlmod) #<--How much memory does the fitted MxModel take up?:

#If Newton-Raphson doesn't reach a good solution, try again with SLSQP:
if( !(gremlmod$output$status$code %in% c(0,1)) ){
	gremlmod$compute <- mxComputeSequence(steps=list(
		mxComputeGradientDescent(engine="SLSQP",verbose=5L),
		#^^^Note:  If you are running the R GUI under Windows, delete the 'verbose=5L' argument in the above.
		mxComputeOnce('fitfunction', c('gradient','hessian')),
		mxComputeStandardError(),
		mxComputeHessianQuality(),
		mxComputeReportDeriv(),
		mxComputeReportExpectation()
	))
	gremlmod <- mxRun(gremlmod)
	gc()
}

summary(gremlmod, verbose=T)
truevals
gremlmod$output$fit
gremlmod$output$gradient
gremlmod$output$hessian

