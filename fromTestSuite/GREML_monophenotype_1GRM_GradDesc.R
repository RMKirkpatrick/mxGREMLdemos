# Copyright 2019-2025 by Robert M. Kirkpatrick
# Licensed under CC BY 4.0 <http://creativecommons.org/licenses/by/4.0/>

#This script is adapted from inst/models/nightly/GREML_monophenotype_1GRM_GradDesc.R in the OpenMx source repository.

#This script demonstrates the use of the GREML feature in a simple but semi-realistic example;
#the syntax is further simplified by using the default gradient-descent optimizer instead of Newton-Raphson.
#It first simulates a genomic-relatedness matrix (GRM), a phenotype, and a null covariate.  Then, it
#fits a simple GREML model to estimate additive-genetic variance, residual variance, and heritability.


require(OpenMx)
mxOption(NULL,"Default optimizer","SLSQP")
#With more threads, the job will run more quickly, but will require more memory:
mxOption(NULL,"Number of Threads",2)
options(mxCondenseMatrixSlots=TRUE)  #<--Saves memory
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
ge <- mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T,REML=FALSE)
#^^^Ordinary rather than restricted maximum-likelihood, in order to enable comparison with 
#the other script in this subdirectory.

#The GREML fitfunction object:
gff <- mxFitFunctionGREML()

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
	gff #<--GREML fitfunction
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

