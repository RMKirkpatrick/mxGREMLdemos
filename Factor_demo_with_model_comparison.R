# Copyright 2019-2021 by Robert M. Kirkpatrick
# Licensed under CC BY 4.0 <http://creativecommons.org/licenses/by/4.0/>

# This is an adaptation of a script previously used at the 2017 AGES Workshop at Virginia Commonwealth University.

#This script demonstrates the use of the GREML feature in a factor model.  There are three 
#phenotypes and 1 common factor; for the sake of simplicity, the only heritable variance is common variance (i.e., in the common factor).
#The script first loads a genomic-relatedness matrix (GRM) and a file containing the phenotype and a covariate.  Then, it 
#fits a simple GREML model to estimate additive-genetic variance, residual variance, and heritability.
#The GRM is 500x500, and was computed from 50,000 simulated SNPs in linkage equilibrium.

require(OpenMx)
options(mxCondenseMatrixSlots=TRUE)  #<--Saves memory
mxOption(NULL,"Analytic Gradients","Yes")
#You need to set R's working directory to the directory containing the data files for this demo.
#(i.e., YOU MUST CHANGE THE NEXT LINE TO REFLECT WHERE, ON YOUR COMPUTER, YOU'VE PLACED THE DATA FILES):
setwd("/home/rmk/OpenMx_dev/GREML_demos/repo/mxGREMLdemos/AGES2017/data")

N <- 500 #<--Total number of participants.
#^^^We're using a small sample size for the sake of making the MxModel run quickly.  500 individuals is too small
#to have any real power with a model like this one.


# Load and check data: ##################################################################

#Load the GRM.  Argument 'prefix' to omxReadGRMBin() should be everything in the filename and path of one of the GRM's
#that precedes the ".grm" :
grmstuff <- omxReadGRMBin(prefix="AGES2017FactorModGRM",returnList=T)
closeAllConnections()
str(grmstuff)
#^^^grmstuff will be a list of 4 elements:
#     $diag is a numeric vector containing the GRM's diagonal elements.
#     $off is a numeric vector containing the GRM's off-diagonal elements.
#     $id is a dataframe that contains the family and individual IDs corresponding to the rows and columns of the GRM.
#     $N is the number of markers used to compute the GRM.

#Use the elements of grmstuff to make a proper GRM:
#Initialize:
GRM <- matrix(
	0.0,nrow=N,ncol=N,
	dimnames=list(grmstuff$id$V1+grmstuff$id$V2,grmstuff$id$V1+grmstuff$id$V2))
#^^^We're giving the GRM dimnames based on IDs.
#Populate the upper triangle of the GRM, excluding diagonal, with grmstuff$off:
GRM[!lower.tri(GRM,diag=T)] <- grmstuff$off
#Populate the lower triangle, excluding diagonal, of the GRM:
GRM <- GRM + t(GRM)
#Populate the diagonal of the GRM:
diag(GRM) <- grmstuff$diag
#Look at upper left corner of GRM:
GRM[1:5,1:5]

#Load phenotype file:
phenfile <- read.csv("GREMLFactorModData.csv",header=T)
#The columns are family ID, individual ID, phenotype #1, phenotype #2, phenotype #3, a covariate relevant to all
#3 phenotypes ('x0'), a covariate relevant to phenotype #1, a covariate relevant to phenotype #2, and a covariate
#relevant to phenotype #3:
head(phenfile)
#There is one missing value on each phenotype:
colSums(is.na(phenfile))
#Give phenfile rownames based on IDs:
rownames(phenfile) <- phenfile$famid + phenfile$id
#Check to make sure that there are no individuals represented in the GRM but not in the phenotype file (there won't
#be, since the data are simulated, but this is a good check to make when working with real data):
all(rownames(GRM) %in% rownames(phenfile)) #<--Should return TRUE.
#And check vice versa:
all(rownames(phenfile) %in% rownames(GRM)) #<--Should also be TRUE.
#In fact, the rows of phenfile are matched with those of the GRM:
all(rownames(phenfile)==rownames(GRM)) #<--TRUE.


# Use mxGREMLDataHandler(): ################################################

#We give the data-handler our dataset, tell it which columns are the phenotype, and tell it which covariates go
#with each phenotype:
gremldat <- mxGREMLDataHandler(
	data=phenfile,yvars=c("y1","y2","y3"),
	Xvars=list(c("x0","x1"),c("x0","x2"),c("x0","x3")),
	addOnes=T,blockByPheno=T,staggerZeroes=T #<--FYI, these are the defaults for these arguments.
)
#We can see that three datapoints have been dropped from the data, due to missing data:
str(gremldat)
#This is what the processed dataset looks like.  In the column names of matrix X, the part before the underscore
#identifies the phenotype being regressed against the covariate, and the part after the underscore identifies
#the covariate ('1' just refers to the constant, for the regression intercept):
head(gremldat$yX)
#Here we can see where phenotype #1 stops and phenotype #2 starts in y:
gremldat$yX[497:502,]
#Here we can see where phenotype #2 stops and phenotype #3 starts in y:
gremldat$yX[996:1001,]
casesToDrop <- gremldat$casesToDrop

#Obtain OLS residual variance from a regression of y onto X (for use with start values):
orv <- var(lm(gremldat$yX[,1] ~ gremldat$yX[,-1])$residual)


# Create MxData object, expectation, fitfunction, and compute plan: ##########################################

#The MxData object.  We will use gremldat$yX as the raw dataset to be used in our MxModel:
mxdat <- mxData(observed=gremldat$yX,type="raw",sort=FALSE)
#^^^N.B. the 'sort=FALSE' argument! 

#We  tell the GREML espectation that the name of the model-expected covariance matrix is "V", and that the dataset
#being input is y horizontally adhered to X (and therefore, it should not call the data-handler at runtime):
ge <- mxExpectationGREML(V="V",dataset.is.yX=TRUE,casesToDropFromV=casesToDrop)

#We tell the GREML fitfunction the names of the first partial derivatives of 'V' w/r/t each free parameter:
gff <- mxFitFunctionGREML(
	dV=c(l1="dV_dl1",l2="dV_dl2",l3="dV_dl3",va="dV_dva",vu1="dV_dvu1",vu2="dV_dvu2",vu3="dV_dvu3"))

#Custom compute plan, to use Newton-Raphson and see what the optimizer is doing in real time (via the 
# 'verbose' argument):
plan <- mxComputeSequence(
	steps=list(
		mxComputeNewtonRaphson(verbose=5L),
		mxComputeOnce('fitfunction', c('gradient','hessian')),
		mxComputeStandardError(),
		mxComputeHessianQuality(),
		mxComputeReportDeriv(),
		mxComputeReportExpectation()
	))
#^^^Note:  If you are running the R GUI under Windows, delete the 'verbose=5L' argument in the above.



# Create the MxModel for the factor model: #########################################################
#(We will create a lot of objects inside the mxModel() statement.  This helps to save memory.)

factorMod <- mxModel(
	"GREMLfactordemo",
	mxdat,
	ge,
	gff,
	plan,
	#This will be our matrix of factor loadings:
	mxMatrix(type="Full",nrow=3,ncol=1,free=T,values=sqrt(orv/2),
					 labels=c("l1","l2","l3"),
					 name="Lambda"),
	#1x1 matrix to hold heritability of the common factor.  Note that, to identify the model, we restrict this
	#parameter to the interval [0,1]:
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.5,labels="va",lbound=0,ubound=1,name="Va"),
	#Column vector containing unique variances.  To set start values, we guess that (conditional on X) about half the 
	#variance in each phenotype is unique:
	mxMatrix(type="Full",nrow=3*N,ncol=1,free=T,values=orv/2,
					 labels=c(rep("vu1",N),rep("vu2",N),rep("vu3",N)),lbound=0.0001,name="Vu"),
	#The GRM:
	mxMatrix(type="Symm",nrow=N,free=F,values=GRM,name="A"),
	#NxN matrix of zeroes:
	mxMatrix(type="Zero",nrow=N,ncol=N,name="Zip"),
	#NxN identity matrix:
	mxMatrix(type="Iden",nrow=N,name="I"),
	#Covariance matrix of the common factor:
	mxAlgebra( (A%x%Va) + (I%x%(1-Va)), name="SigmaFac"),
	#Outer product of Lambda with itself:
	mxAlgebra( Lambda%*%t(Lambda), name="LamLamT"),
	#The model-expected covariance matrix, per the Fundamental Theorem of Factor Analysis:
	mxAlgebra(LamLamT%x%SigmaFac + vec2diag(Vu), name="V"),
	#First partial derivative of V w/r/t va:
	mxAlgebra(LamLamT%x%(A-I), name="dV_dva"),
	#First partial derivative of V w/r/t vu1:
	mxAlgebra(rbind(
		cbind(I,Zip,Zip),
		cbind(Zip,Zip,Zip),
		cbind(Zip,Zip,Zip)
	), name="dV_dvu1"),
	#First partial derivative of V w/r/t vu2:
	mxAlgebra(rbind(
		cbind(Zip,Zip,Zip),
		cbind(Zip,I,Zip),
		cbind(Zip,Zip,Zip)
	), name="dV_dvu2"),
	#First partial derivative of V w/r/t vu3:
	mxAlgebra(rbind(
		cbind(Zip,Zip,Zip),
		cbind(Zip,Zip,Zip),
		cbind(Zip,Zip,I)
	), name="dV_dvu3"),
	#First partial derivative of LamLamT w/r/t l1:
	mxAlgebra(rbind(
		cbind(2*Lambda[1,1],Lambda[2,1],Lambda[3,1]),
		cbind(Lambda[2,1],0,0),
		cbind(Lambda[3,1],0,0)),
		name="dLamLamT_dl1"),
	#First partial derivative of V w/r/t l1:
	mxAlgebra(dLamLamT_dl1%x%SigmaFac, name="dV_dl1"),
	#First partial derivative of LamLamT w/r/t l2:
	mxAlgebra(rbind(
		cbind(0,Lambda[1,1],0),
		cbind(Lambda[1,1],2*Lambda[2,1],Lambda[3,1]),
		cbind(0,Lambda[3,1],0)),
		name="dLamLamT_dl2"),
	#First partial derivative of V w/r/t l2:
	mxAlgebra(dLamLamT_dl2%x%SigmaFac, name="dV_dl2"),
	#First partial derivative of LamLamT w/r/t l3:
	mxAlgebra(rbind(
		cbind(0,0,Lambda[1,1]),
		cbind(0,0,Lambda[2,1]),
		cbind(Lambda[1,1],Lambda[2,1],2*Lambda[3,1])),
		name="dLamLamT_dl3"),
	#First partial derivative of V w/r/t l3:
	mxAlgebra(dLamLamT_dl3%x%SigmaFac, name="dV_dl3")
)

#Run and view summary:
factorRun <- mxRun(factorMod)
object.size(factorRun) #<--How much memory does the fitted MxModel take up?:
summary(factorRun)

## Now, we fit an unstructured model, using direct-symmetric parameterization.  ##############################
# This is the most "saturated" model we could fit in this context.

#The next two options together tell NPSOL to write a log entry for each of its major iterations to a file in   
#the working directory, called 'fort.2' (change the value of "Print file" option to a different positive integer to change the extension 
#to 'fort'):
# mxOption(NULL,"Print level",20)
# mxOption(NULL,"Print file",2)

#The following option will tell NPSOL to check the analytic gradient elements against its numerical derivatives.  If the analytic derivatives
#clearly appear incorrect, NPSOL will terminate, with status code 7.  You need the "Print level" option above to see any details, though.
#Useful if you modify the derivatives of the covariance matrix, but otherwise not worth the added computational effort:
# mxOption(NULL,"Verify level",3)
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


directVCMod <- mxModel(
	"DirectSymm",
	mxdat,
	ge,
	plan,
	#The GRM:
	mxMatrix(type="Symm",nrow=N,free=F,values=GRM,name="A"),
	#NxN identity matrix:
	mxMatrix(type="Iden",nrow=N,name="I"),
	# 3x3 within-person additive-genetic covariance matrix:
	mxMatrix(
		type="Symm",nrow=3,ncol=3,free=T,
		values=c(
			orv/2,0,0,
			0,orv/2,0,
			0,0,orv/2),
		labels=c(
			"a11","a12","a13",
			"a12","a22","a23",
			"a13","a23","a33"),name="covA"),
	# 3x3 within-person nonshared-environmental covariance matrix:
	mxMatrix(
		type="Symm",nrow=3,ncol=3,free=T,
		values=c(
			orv/2,0,0,
			0,orv/2,0,
			0,0,orv/2),
		labels=c(
			"e11","e12","e13",
			"e12","e22","e23",
			"e13","e23","e33"),name="covE"),
	#Model-expected covariance matrix, V:
	mxAlgebra( (covA%x%A) + (covE%x%I), name="V" ),
	#^^^This is not the most computationally efficient way to compute V, but it's done here for didactic clarity.
	
	#Derivatives of 'covA' w/r/t the parameters that appear in it:
	mxMatrix(
		type="Symm",nrow=3,ncol=3,free=F,values=c(
			1,0,0,
			0,0,0,
			0,0,0
		),name="dcovA_da11"),
	mxMatrix(
		type="Symm",nrow=3,ncol=3,free=F,values=c(
			0,1,0,
			1,0,0,
			0,0,0
		),name="dcovA_da12"),
	mxMatrix(
		type="Symm",nrow=3,ncol=3,free=F,values=c(
			0,0,1,
			0,0,0,
			1,0,0
		),name="dcovA_da13"),
	mxMatrix(
		type="Symm",nrow=3,ncol=3,free=F,values=c(
			0,0,0,
			0,1,0,
			0,0,0
		),name="dcovA_da22"),
	mxMatrix(
		type="Symm",nrow=3,ncol=3,free=F,values=c(
			0,0,0,
			0,0,1,
			0,1,0
		),name="dcovA_da23"),
	mxMatrix(
		type="Symm",nrow=3,ncol=3,free=F,values=c(
			0,0,0,
			0,0,0,
			0,0,1
		),name="dcovA_da33"),
	#Derivatives of 'covE' w/r/t the parameters that appear in it:
	mxMatrix(
		type="Symm",nrow=3,ncol=3,free=F,values=c(
			1,0,0,
			0,0,0,
			0,0,0
		),name="dcovE_de11"),
	mxMatrix(
		type="Symm",nrow=3,ncol=3,free=F,values=c(
			0,1,0,
			1,0,0,
			0,0,0
		),name="dcovE_de12"),
	mxMatrix(
		type="Symm",nrow=3,ncol=3,free=F,values=c(
			0,0,1,
			0,0,0,
			1,0,0
		),name="dcovE_de13"),
	mxMatrix(
		type="Symm",nrow=3,ncol=3,free=F,values=c(
			0,0,0,
			0,1,0,
			0,0,0
		),name="dcovE_de22"),
	mxMatrix(
		type="Symm",nrow=3,ncol=3,free=F,values=c(
			0,0,0,
			0,0,1,
			0,1,0
		),name="dcovE_de23"),
	mxMatrix(
		type="Symm",nrow=3,ncol=3,free=F,values=c(
			0,0,0,
			0,0,0,
			0,0,1
		),name="dcovE_de33"),
	
	#Derivatives of V w/r/t free parameters:
	mxAlgebra(dcovA_da11%x%A, name="dV_da11"),
	mxAlgebra(dcovA_da12%x%A, name="dV_da12"),
	mxAlgebra(dcovA_da13%x%A, name="dV_da13"),
	mxAlgebra(dcovA_da22%x%A, name="dV_da22"),
	mxAlgebra(dcovA_da23%x%A, name="dV_da23"),
	mxAlgebra(dcovA_da33%x%A, name="dV_da33"),
	mxAlgebra(dcovE_de11%x%I, name="dV_de11"),
	mxAlgebra(dcovE_de12%x%I, name="dV_de12"),
	mxAlgebra(dcovE_de13%x%I, name="dV_de13"),
	mxAlgebra(dcovE_de22%x%I, name="dV_de22"),
	mxAlgebra(dcovE_de23%x%I, name="dV_de23"),
	mxAlgebra(dcovE_de33%x%I, name="dV_de33"),
	#^^^NOTICE THAT ALL OF THEM ARE CONSTANT W/R/T FREE PARAMETERS!
	#That is, they can be calculated once, and need not be recalculated as the optimizer changes free parameter values.
	
	#GREML Fitfunction
	mxFitFunctionGREML(
		dV=c(
			a11="dV_da11",a12="dV_da12",a13="dV_da13",a22="dV_da22",a23="dV_da23",a33="dV_da33",
			e11="dV_de11",e12="dV_de12",e13="dV_de13",e22="dV_de22",e23="dV_de23",e33="dV_de33"))
)
directVCRun <- mxRun(directVCMod)
summary(directVCRun)

## Now, we fit the unstructured model, using Cholesky parameterization.  ##############################

cholMod <- mxModel(
	"Cholesky",
	mxdat,
	ge,
	plan,
	#The GRM:
	mxMatrix(type="Symm",nrow=N,free=F,values=GRM,name="A"),
	#NxN identity matrix:
	mxMatrix(type="Iden",nrow=N,name="I"),
	#Lower-triangular matrix containing additive-genetic paths:
	mxMatrix(
		type="Lower",nrow=3,ncol=3,free=T,
		values=c(sqrt(orv/2),0,0,sqrt(orv/2),0,sqrt(orv/2)),
		labels=c("a11","a12","a13","a22","a23","a33"),name="a"),
	#Lower-triangular matrix containing nonshared-environmental paths:
	mxMatrix(
		type="Lower",nrow=3,ncol=3,free=T,
		values=c(sqrt(orv/2),0,0,sqrt(orv/2),0,sqrt(orv/2)),
		labels=c("e11","e12","e13","e22","e23","e33"),name="e"),
	# 3x3 within-person additive-genetic covariance matrix:
	mxAlgebra(a%*%t(a),name="covA"),
	# 3x3 within-person nonshared-environmental covariance matrix:
	mxAlgebra(e%*%t(e),name="covE"),
	#Model-expected covariance matrix, V:
	mxAlgebra( (covA%x%A) + (covE%x%I), name="V" ),
	#^^^This is not the most computationally efficient way to compute V, but it's done here for didactic clarity.
	
	#Derivatives of 'covA' w/r/t the parameters that appear in it:
	mxAlgebra(
		rbind(
			cbind(2*a11,a12,a13),
			cbind(a12,0,0),
			cbind(a13,0,0)
		), name="dcovA_da11"),
	mxAlgebra(
		rbind(
			cbind(0,a11,0),
			cbind(a11,2*a12,a13),
			cbind(0,a13,0)
		), name="dcovA_da12"),
	mxAlgebra(
		rbind(
			cbind(0,0,a11),
			cbind(0,0,a12),
			cbind(a11,a12,2*a13)
		), name="dcovA_da13"),
	mxAlgebra(
		rbind(
			cbind(0,0,0),
			cbind(0,2*a22,a23),
			cbind(0,a23,0)
		), name="dcovA_da22"),
	mxAlgebra(
		rbind(
			cbind(0,0,0),
			cbind(0,0,a22),
			cbind(0,a22,2*a23)
		), name="dcovA_da23"),
	mxAlgebra(
		rbind(
			cbind(0,0,0),
			cbind(0,0,0),
			cbind(0,0,2*a33)
		), name="dcovA_da33"),
	#Derivatives of 'covE' w/r/t the parameters that appear in it:
	mxAlgebra(
		rbind(
			cbind(2*e11,e12,e13),
			cbind(e12,0,0),
			cbind(e13,0,0)
		), name="dcovE_de11"),
	mxAlgebra(
		rbind(
			cbind(0,e11,0),
			cbind(e11,2*e12,e13),
			cbind(0,e13,0)
		), name="dcovE_de12"),
	mxAlgebra(
		rbind(
			cbind(0,0,e11),
			cbind(0,0,e12),
			cbind(e11,e12,2*e13)
		), name="dcovE_de13"),
	mxAlgebra(
		rbind(
			cbind(0,0,0),
			cbind(0,2*e22,e23),
			cbind(0,e23,0)
		), name="dcovE_de22"),
	mxAlgebra(
		rbind(
			cbind(0,0,0),
			cbind(0,0,e22),
			cbind(0,e22,2*e23)
		), name="dcovE_de23"),
	mxAlgebra(
		rbind(
			cbind(0,0,0),
			cbind(0,0,0),
			cbind(0,0,2*e33)
		), name="dcovE_de33"),
	
	#N.B. that all of the derivatives of 'covA' and 'covE' CONTAIN FREE PARAMETERS!
	#That is in contrast to the direct-symmetric specification above.
	
	#Derivatives of V w/r/t free parameters:
	mxAlgebra(dcovA_da11%x%A, name="dV_da11"),
	mxAlgebra(dcovA_da12%x%A, name="dV_da12"),
	mxAlgebra(dcovA_da13%x%A, name="dV_da13"),
	mxAlgebra(dcovA_da22%x%A, name="dV_da22"),
	mxAlgebra(dcovA_da23%x%A, name="dV_da23"),
	mxAlgebra(dcovA_da33%x%A, name="dV_da33"),
	mxAlgebra(dcovE_de11%x%I, name="dV_de11"),
	mxAlgebra(dcovE_de12%x%I, name="dV_de12"),
	mxAlgebra(dcovE_de13%x%I, name="dV_de13"),
	mxAlgebra(dcovE_de22%x%I, name="dV_de22"),
	mxAlgebra(dcovE_de23%x%I, name="dV_de23"),
	mxAlgebra(dcovE_de33%x%I, name="dV_de33"),
	#^^^Because the 'dcov*' algebras contain free parameters, these derivatives of V must be recalculated whenever the 
	#parameters change during optimization!
	
	#GREML fitfunction:
	mxFitFunctionGREML(
		dV=c(
			a11="dV_da11",a12="dV_da12",a13="dV_da13",a22="dV_da22",a23="dV_da23",a33="dV_da33",
			e11="dV_de11",e12="dV_de12",e13="dV_de13",e22="dV_de22",e23="dV_de23",e33="dV_de33"))
)
cholRun <- mxRun(cholMod)

if(cholRun$output$status$code > 1){
	cholRun$compute <- mxComputeSequence(
		steps=list(
			#mxComputeGradientDescent(engine="NPSOL", verbose=5L),
			mxComputeGradientDescent(engine="CSOLNP", verbose=5L),
			mxComputeOnce('fitfunction', c('gradient','hessian')),
			mxComputeStandardError(),
			mxComputeHessianQuality(),
			mxComputeReportDeriv(),
			mxComputeReportExpectation()
		))
	cholRun <- mxRun(cholRun)
}

if(cholRun$output$status$code > 1){
	mxOption(NULL,"Analytic Gradients","No")
	cholRun$compute <- mxComputeSequence(
		steps=list(
			mxComputeGradientDescent(engine="SLSQP", verbose=5L),
			mxComputeNumericDeriv(),
			mxComputeStandardError(),
			mxComputeHessianQuality(),
			mxComputeReportDeriv(),
			mxComputeReportExpectation()
		))
	cholRun <- mxRun(cholRun)
}

summary(cholRun)

#It took at least 2 tries to get the Cholesky model to converge, and its solution still has a worse -2logL
#than the direct-symmetric model:
cholRun$output$fit
directVCRun$output$fit

#The Cholesky model's running time is 5 times slower than the direct-symmetric model's,
#and that's ignoring the failed attempt(s) to fit the Cholesky!:
cholRun$output$wallTime
directVCRun$output$wallTime

#Finally, compare the factor model to the "saturated" model:
mxCompare(directVCRun, factorRun) #<--Non-significant deterioration of fit; factor model has smaller AIC.

#The reason why the direct-symmetric model can achieve a smaller -2logL than the Cholesky model
#is that the Cholesky is subject to an implicit constraint, to which the direct-symmetric model
#is not subject--the Cholesky parameterization forces `covA` and `covE` to be positive-definite:
eigen(mxEval(covA,cholRun,T))$values #<--All eigenvalues are positive.
eigen(mxEval(covA,directVCRun,T))$values #<--Not all eigenvalues are positive.
