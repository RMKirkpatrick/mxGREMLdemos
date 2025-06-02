# Copyright 2019-2025 by Robert M. Kirkpatrick
# Licensed under CC BY 4.0 <http://creativecommons.org/licenses/by/4.0/>

#This script is adapted from inst/models/nightly/GREML_monophenotype_1GRM.R in the OpenMx source repository.

#This script demonstrates the use of the GREML feature in a simple but semi-realistic example.
#It first simulates a genomic-relatedness matrix (GRM), a phenotype, and a null covariate.  Then, it
#fits a simple GREML model to estimate additive-genetic variance, residual variance, and heritability.


require(OpenMx)
options(mxCondenseMatrixSlots=TRUE)  #<--Saves memory
mxOption(NULL,"Default optimizer","SLSQP")
#With more threads, the job will run more quickly, but will require more memory:
mxOption(NULL,"Number of Threads",2)
require(mvtnorm)


#Generate data:
set.seed(476)
A <- matrix(0,1000,1000)  #<--Empty GRM
A[lower.tri(A)] <- runif(499500, -0.025, 0.025)
A <- A + t(A)
diag(A) <- runif(1000,0.95,1.05) #<--GRM now complete
y <- t(rmvnorm(1,sigma=A*0.5))  #<--Phenotype 'y' has a "population" variance of 1 and h2 of 0.5 
y <- y + rnorm(1000,sd=sqrt(0.5))
x <- rnorm(1000) #<--Covariate 'x' is actually independent of the phenotype.
#Merge variables into data matrix:
dat <- cbind(y,x)
colnames(dat) <- c("y","x") #<--Column names

#The GREML expectation tells OpenMx that the model-expected covariance matrix is named 'V', that the one 
#phenotype is has column label 'y' in the dataset, that the one covariate has column label 'x' in the dataset,
#and that a lead column of ones needs to be appended to the 'X' matrix (for the intercept):
ge <- mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T, REML=F)

#The GREML fitfunction tells OpenMx that the derivative of 'V' with respect to free parameter 
#'va'(the additive-genetic variance) is a matrix named 'A', and that the derivative of 'V' w/r/t free parameter
#'#'ve' is a matrix named 'I'.  At runtime, the GREML fitfunction will use these derivatives to help with 
#'#optimization and compute standard errors:
gff <- mxFitFunctionGREML(dV=c(va="A",ve="I"))

#This is a custom compute plan.  It is necessary because we want to use the Newton-Raphson optimizer, which
#can use analytic first and second derivatives of the GREML fitfunction to speed up convergence.  It looks
#especially messy here because we want a profile-likelihood confidence interval for the heritability:
plan <- mxComputeSequence(steps=list(
	mxComputeNewtonRaphson(),
	mxComputeOnce('fitfunction', c('gradient','hessian')),
	mxComputeConfidenceInterval(
		plan=mxComputeGradientDescent(nudgeZeroStarts=FALSE, maxMajorIter=150)),
	mxComputeStandardError(),
	mxComputeReportDeriv(),
	mxComputeReportExpectation()
))

#The MxData object.  N.B. use of 'sort=FALSE' is CRITICALLY IMPORTANT, because the rows and columns of dataset
#'dat' and the rows and columns of GRM 'A' are already properly ordered:
mxdat <- mxData(observed = dat, type="raw", sort=FALSE)

#We will create some of the necessary objects inside the mxModel() statement.  We mainly want to avoid creating 
#more copies of the GRM than we need to:
testmod <- mxModel(
	"GREML_1GRM_1trait", #<--Model name
	#1x1 matrix containing residual variance component:
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = var(y)/2, labels = "ve", lbound = 0.0001, 
					 name = "Ve"),
	#1x1 matrix containing additive-genetic variance component:
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = var(y)/2, labels = "va", name = "Va"),
	#1000x1000 identity matrix--the "relatedness matrix" for the residuals:
	mxMatrix("Iden",nrow=1000,name="I"),
	#The GRM:
	mxMatrix("Symm",nrow=1000,free=F,values=A,name="A"),
	#The model-expected covariance matrix:
	mxAlgebra((A%x%Va) + (I%x%Ve), name="V"),
	#An MxAlgebra for the heritability:
	mxAlgebra(Va/(Va+Ve), name="h2"),
	mxCI("h2"), #<--Request confidence interval for heritability
	mxdat, #<--MxData object
	ge, #<--GREML expectation
	gff, #<--GREML fitfunction
	plan #<--Custom compute plan
)

testrun <- mxRun(testmod,intervals = T) #<--Run model
summary(testrun) #<--Model summary

#Obtain SE of h2 from delta-method approximation (e.g., Lynch & Walsh, 1998, Appendix 1):
scm <- testrun$output$vcov #<--Sampling covariance matrix for ve and va
pointest <- testrun$output$estimate #<--Point estimates of ve and va
h2se <- sqrt(
	(pointest[2]/(pointest[1]+pointest[2]))^2 * (
		(scm[2,2]/pointest[2]^2) - (2*scm[1,2]/pointest[1]/(pointest[1]+pointest[2])) + 
			(sum(scm)*(pointest[1]+pointest[2])^-2)
))
#Compare:
mxEval(h2,testrun,T)[1,1] + 2*c(-h2se,h2se)
testrun$output$confidenceIntervals


#Diagonalize the problem: ###
eigenA <- eigen(A) #<--Eigen decomposition of the GRM
#We "rotate out" the dependence among participants by premultiplying the 'y' and 'X' matrices by the 
#eigenvectors of the GRM:
yrot <- t(eigenA$vectors) %*% y
xrot <- t(eigenA$vectors) %*% cbind(1,x)
datrot <- cbind(yrot,xrot,eigenA$values)
colnames(datrot) <- c("y","x0","x1","EigVal")
#Make a new MxModel:
testmod2 <- mxModel(
	"GREMLtest_1GRM_1trait_diagonalized",
	mxData(observed = datrot, type="raw", sort=FALSE),
	mxMatrix(
		type = "Full", nrow = 1, ncol=1, free=T, values = var(y)/2, labels = "ve", lbound = 0.0001, 
		name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = var(y)/2, labels = "va", name = "Va"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=mean(y),labels="b0",name="B0"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.1,labels="b1",name="B1"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=F,labels="data.EigVal",name="Eig"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=F,labels="data.x0",name="X0"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=F,labels="data.x1",name="X1"),
	mxAlgebra(B0*X0 + B1*X1,name="Mu"),
	mxAlgebra(Ve + Va*Eig,name="Sigma"),
	mxExpectationNormal(covariance="Sigma",means="Mu",dimnames="y"),
	mxFitFunctionML(),
	mxAlgebra(Va/(Va+Ve), name="h2"),
	mxCI("h2")
)

testrun2 <- mxRun(testmod2,intervals=T)
#Results are substantially equivalent to those from the previous MxModel:
summary(testrun2)

# Reparametrize the problem in terms of total variance and heritability: ###

gff3 <- mxFitFunctionGREML(dV=c(h2="dVdH2",vp="dVdVp")) #<--Need new fitfunction object
#Need new compute plan:
plan3 <- mxComputeSequence(steps=list(
	mxComputeNewtonRaphson(),
	mxComputeOnce('fitfunction', c('gradient','hessian')),
	mxComputeConfidenceInterval(
		plan=mxComputeGradientDescent(nudgeZeroStarts=FALSE, maxMajorIter=150)),
	mxComputeStandardError(),
	mxComputeReportDeriv(),
	mxComputeReportExpectation()
))

testmod3 <- mxModel(
	"GREML_1GRM_1trait_altparam", #<--Model name
	#1x1 matrix containing heritability:
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.5, labels = "h2", lbound = 0.0001, ubound=0.9999,
					 name = "H2"),
	#1x1 matrix containing total phenotypic variance:
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = var(y), labels = "vp", name = "Vp"),
	#1000x1000 identity matrix--the "relatedness matrix" for the residuals:
	mxMatrix("Iden",nrow=1000,name="I"),
	#The GRM:
	mxMatrix("Symm",nrow=1000,free=F,values=A,name="A"),
	#MxAlgebra for additive-genetic variance:
	mxAlgebra(H2*Vp, name="Va"),
	#MxAlgebra for residual variance:
	mxAlgebra((1-H2)*Vp, name="Ve"),
	#The model-expected covariance matrix:
	mxAlgebra( (A%x%Va) + (I%x%Ve), name="V"),
	#MxAlgebras for derivatives of V w/r/t free parameters:
	mxAlgebra((A-I)%x%Vp, name="dVdH2"),
	mxAlgebra(I + (A-I)%x%H2, name="dVdVp"),
	mxCI("h2"), #<--Request confidence interval for heritability
	mxdat, #<--MxData object
	ge, #<--GREML expectation
	gff3, #<--GREML fitfunction
	plan3 #<--Custom compute plan
)

testrun3 <- mxRun(testmod3, intervals = T)
summary(testrun3)

#Compare:
mxEval(h2,testrun3,T)[1,1] + 2*c(-0.07866897,0.07866897) #<--0.07866897 is the SE of h2
testrun3$output$confidenceIntervals
