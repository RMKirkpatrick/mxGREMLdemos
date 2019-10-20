# Copyright 2019 by Robert M. Kirkpatrick
# Licensed under CC BY 4.0 <http://creativecommons.org/licenses/by/4.0/>

# This is an adaptation of a script previously used at the 2017 AGES Workshop at Virginia Commonwealth University.

#This script demonstrates the use of the GREML feature in a simple but realistic example.
#It first loads a genomic-relatedness matrix (GRM) and a file containing the phenotype and a covariate.  Then, it 
#fits a simple GREML model to estimate additive-genetic variance, residual variance, and heritability.
#The GRM is 1000x1000, and was computed from 50,000 simulated SNPs in linkage equilibrium.
#Note that this analysis is simple enough that it could be done in GCTA (and GCTA would do it more quickly).

require(OpenMx)
options(mxCondenseMatrixSlots=TRUE)  #<--Saves memory
#You need to set R's working directory to the directory containing the data files for this demo.
#(i.e., YOU MUST CHANGE THE NEXT LINE TO REFLECT WHERE, ON YOUR COMPUTER, YOU'VE PLACED THE DATA FILES):
setwd("/home/rmk/OpenMx_dev/GREML_demos/repo/mxGREMLdemos/AGES2017/data")

#Total number of participants:
N <- 1000

# Load and check data: ##################################################################

#Load the GRM's data.  Argument 'prefix' to ReadGRMBin() should be everything in the filename and path of one of the GRM's
#that precedes the ".grm" :
grmstuff <- omxReadGRMBin(prefix="AGES2017GRM",returnList=T)
closeAllConnections()
str(grmstuff)
#^^^grmstuff will be a list of 4 elements:
#     $diag is a numeric vector containing the GRM's diagonal elements.
#     $off is a numeric vector containing the GRM's off-diagonal elements.
#     $id is a dataframe that contains the family and individual IDs orresponding to the rows and columns of the GRM.
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
phenfile <- read.csv("GREMLphenodata.csv",header=T)
#The columns are family ID, individual ID, phenotype, and covariate:
head(phenfile)
#There is one missing value for the phenotype and for the covariate:
sum(is.na(phenfile$y))
sum(is.na(phenfile$x))
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

#Give the data-handler the dataset, and tell it that the one phenotype has column label 'y', 
#that the one covariate has column label 'x', and that a lead column of ones needs to be appended to the 'X' 
#matrix (for the intercept):
gremldat <- mxGREMLDataHandler(data=phenfile, yvars="y", Xvars="x", addOnes=T)
#We can see that participants #7 and #77 in phenfile were dropped due to missing data:
str(gremldat)
#We can see that the first column of gremldata$yX is phenotype vector y, and that the remaining columns constitute
#the X matrix--a regression constant and the covariate:
head(gremldat$yX)
#Number of datapoints, post-data-handling:
Np <- nrow(gremldat$yX)
#Drop rows and columns #7 and #77 from the GRM:
GRM <- GRM[-c(7,77),-c(7,77)]


# Create MxData object, expectation, fitfunction, and compute plan: ##########################################

#The MxData object.  We will use gremldat$yX as the raw dataset to be used in our MxModel:
mxdat <- mxData(observed=gremldat$yX, type="raw", sort=FALSE)
#^^^N.B. the 'sort=FALSE' argument!  I CANNOT EMPHASIZE ENOUGH HOW IMPORTANT THAT IS!!!  By default, OpenMx 
#automatically sorts the rows of a raw dataset at runtime, to try to achieve a perfomance boost.  However, we 
#currently have the rows of y and X sorted exactly the same as in the GRM, so we do not want OpenMx to do any 
#sorting!

#We  tell the GREML espectation that the name of the model-expected covariance matrix is "V", and that the dataset
#being input is y horizontally adhered to X (and therefore, it should not call the data-handler at runtime):
ge <- mxExpectationGREML(V="V",dataset.is.yX=TRUE)
#^^^Note: mxExpectationGREML() has argument 'casesToDropFromV', to which we could have passed 
#gremldat$casesToDrop.  In that case, we could define our covariance matrix and its derivatives as 
#full-size NxN matrices, and allow the backend to use the value of 'casesToDropFromV' to automatically 
#trim rows and columns in those matrices that correspond to missing observations (i.e., participants #7 and #77).
#However, that backend  matrix-resizing process carries a performance cost.

#The GREML fitfunction tells OpenMx that the derivative of 'V' with respect to free parameter 
# 'va'(the additive-genetic variance) is a matrix named 'A', and that the derivative of 'V' w/r/t free parameter
# 've' is a matrix named 'I'.  At runtime, the GREML fitfunction will use these derivatives to help with 
#optimization and to compute standard errors:
gff <- mxFitFunctionGREML(dV=c(va="A",ve="I"))

#This is a custom compute plan.  It is necessary because we want to use the Newton-Raphson optimizer, which
#can use analytic first and second derivatives of the GREML fitfunction to speed up convergence:
plan <- mxComputeSequence(
	#List of steps.  The first is the most notable.  It tells OpenMx to use the Newton-Raphson optimizer, and to 
	#display what it's doing in real time (via the verbose argument):
	steps=list(
		mxComputeNewtonRaphson(verbose=5L),
		mxComputeOnce("fitfunction", c("gradient","hessian")),
		mxComputeStandardError(),
		mxComputeHessianQuality(),
		mxComputeReportDeriv(),
		mxComputeReportExpectation()
))
#^^^Note:  If you are running the R GUI under Windows, delete the 'verbose=5L' argument in the above.


# Create the MxModel, and run it: ###############################################################################
#(We will create a lot of objects inside the mxModel() statement.  This helps to save memory.)

testmod <- mxModel(
	"GREML_1GRM_1trait", #<--Model name
	#1x1 matrix containing nonshared-environmental (residual) variance component:
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = var(gremldat$yX[,"y"])/2, labels = "ve", 
					 lbound = 0.0001, name = "Ve"),
	#1x1 matrix containing additive-genetic variance component:
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = var(gremldat$yX[,"y"])/2, labels = "va", 
					 lbound=0, name = "Va"), #<--Note the lower bound on zero.
	#998x998 identity matrix--the "relatedness matrix" for the residuals:
	mxMatrix("Iden",nrow=Np,name="I"),
	#The GRM:
	mxMatrix("Symm",nrow=Np,free=F,values=GRM,name="A"),
	#The model-expected covariance matrix:
	mxAlgebra((A%x%Va) + (I%x%Ve), name="V"),
	#An MxAlgebra for the heritability:
	mxAlgebra(Va/(Va+Ve), name="h2"),
	mxdat, #<--MxData object
	ge, #<--GREML expectation
	gff, #<--GREML fitfunction
	plan #<--Custom compute plan
)

#Run and view summary:
testrun <- mxRun(testmod)
summary(testrun)
mxEval(h2,testrun,T)

#If we had added the appropriate additional step to the compute plan, we could have requested 
#confidence intervals for va, ve, and/or h2 (but not for the regression coefficients).
#We have asymptotic-ML-theory standard errors for va and ve.  It is possible to obtain an SE for h2
#using a delta-method approximation, though its accuracy is quesionable given the fairly small N:
scm <- chol2inv(chol(testrun$output$hessian/2)) #<--Sampling covariance matrix for ve and va
pointest <- testrun$output$estimate #<--Point estimates of ve and va
h2se <- sqrt(
	(pointest[2]/(pointest[1]+pointest[2]))^2 * (
		(scm[2,2]/pointest[2]^2) - (2*scm[1,2]/pointest[1]/(pointest[1]+pointest[2])) + 
			(sum(scm)*(pointest[1]+pointest[2])^-2)
	))
print(h2se)

#What follows is a somewhat questionable likelihood-ratio test of the null hypothesis that the genetic variance is zero.
#It's questionable because the "full" model was fitted by restricted ML, but the "reduced" model was fitted by OLS, and 
#OLS is equivalent to ordinary ML (provided that the #errors are IID normal with expectation zero).  We use the "full" model's
#ordinary ML fitfunction value, but that's not the exact fitfunction that was optimized to fit it.
#The rationale is that, if there is no genetic variance, then our model reduces to an OLS regression of the phenotype onto the covariate:
lmod <- lm(y~x,data=as.data.frame(gremldat$yX))
#Obtain a likelihood-ratio test statistic:
LRT <- as.numeric(-2*logLik(lmod) - testrun$fitfunction$MLfit)
#Recall that we set a lower bound of zero on the genetic variance. Zero is also the value of the genetic variance
#under the null hypothesis.  We are effectively doing a one-sided hypothesis test.  Under the null hypothesis,
#the LRT follows a mixture distribution, the 0.95 quantile of which is about 2.707:
LRT > 2.707 #<--Significant, at significance level 0.05.
1-(0.5+0.5*pchisq(LRT,1)) #<--One-sided p-value.


