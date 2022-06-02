# Copyright 2019-2022 by Robert M. Kirkpatrick
# Licensed under CC BY 4.0 <http://creativecommons.org/licenses/by/4.0/>

# This is an adaptation of a script previously used at the 2017 AGES Workshop at Virginia Commonwealth University.

#This script demonstrates the use of the GREML feature in a factor model.  There are three 
#phenotypes and 1 common factor; for the sake of simplicity, the only heritable variance is common variance (i.e., in the common factor).
#The script first loads a genomic-relatedness matrix (GRM) and a file containing the phenotype and a covariate.  Then, it 
#fits a simple GREML model to estimate additive-genetic variance, residual variance, and heritability.
#The GRM is 500x500, and was computed from 50,000 simulated SNPs in linkage equilibrium.

require(OpenMx)
options(mxCondenseMatrixSlots=TRUE)  #<--Saves memory
mxOption(NULL,"Default optimizer","SLSQP")
#More threads means faster running time, but at the cost of higher memory demand:
mxOption(NULL,"Number of threads",2)
#You need to set R's working directory to the directory containing the data files for this demo.
#(i.e., YOU MUST CHANGE THE NEXT LINE TO REFLECT WHERE, ON YOUR COMPUTER, YOU'VE PLACED THE DATA FILES):
setwd("./AGES2017/data")

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


# Create the MxModel: ################################################################################
#(We will create a lot of objects inside the mxModel() statement.  This helps to save memory.)

factorMod <- mxModel(
	"GREMLfactordemo",
	mxdat,
	ge,
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
	mxFitFunctionGREML()
)

#Once the MxModel object has been created, we can delete objects in R's workspace which have been copied into the 
#MxModel, or which we simply don't need anymore.  This saves memory:
rm(GRM,ge,gremldat,mxdat,grmstuff)
#This tells R to do "garbage collection":
gc()

#Run and view summary:
factorRun <- mxRun(factorMod)
object.size(factorRun) #<--How much memory does the fitted MxModel take up?:
summary(factorRun)






