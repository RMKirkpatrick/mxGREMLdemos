# Copyright 2019-2025 by Robert M. Kirkpatrick
# Licensed under CC BY 4.0 <http://creativecommons.org/licenses/by/4.0/>

# IMPORTANT!  This script is meant to be run under a Unix-like operating system (e.g., Linux/GNU, macOS, Solaris).
# It will not work properly under Microsoft Windows!

#This script demonstrates the use of OpenMx's GREML feature.
#The script first loads a genomic-relatedness matrix (GRM) and a file containing the phenotypes.  Then, it 
#fits an unstructured GREML model to estimate the additive-genetic and nonshared-environment covariance matrices,
#as well as heritability coefficients and genetic correlations.
#The GRM is 1000x1000, and was computed from 50,000 simulated SNPs in linkage equilibrium.
#After fitting the MxModel, the script invokes GCTA from inside R, to fit the same model via GCTA.
#OpenMx and GCTA obtain substantially equivalent results.

require(OpenMx)
options(mxCondenseMatrixSlots=TRUE)  #<--Saves memory
#More threads means faster running time, but at the cost of higher memory demand:
mxOption(NULL,"Number of threads",2)
#You need to set R's working directory to the directory containing the data files for this demo.
#(i.e., YOU MUST CHANGE THE NEXT LINE TO REFLECT WHERE, ON YOUR COMPUTER, YOU'VE PLACED THE DATA FILES):
setwd("./data")

N <- 1000 #<--Total number of participants.
#^^^We're using a small sample size for the sake of making the MxModel run quickly.

# Load and check data: ##################################################################

#Load the GRM.  Argument 'prefix' to omxReadGRMBin() should be everything in the filename and path of one of the GRM's
#that precedes the ".grm" :
grmstuff <- omxReadGRMBin(prefix="DiphenoGRM",returnList=T)
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
phenfile <- read.csv("DiphenoData.csv",header=T)
#The columns are family ID, individual ID, phenotype #1, and phenotype #2.
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

#We give the data-handler our dataset:
gremldat <- mxGREMLDataHandler(
	data=phenfile,yvars=c("y1","y2"),
	Xvars=list(),
	addOnes=T,blockByPheno=T,staggerZeroes=T #<--FYI, these are the defaults for these arguments.
)
#We can see that 2 datapoints have been dropped from the data, due to missing data:
str(gremldat)
#This is what the processed dataset looks like.  Note that the data-handler created 2 'x' columns,
#which contain the constants for the regression intercepts:
head(gremldat$yX)
#Here we can see where phenotype #1 stops and phenotype #2 starts in y:
gremldat$yX[997:1002,]
casesToDrop <- gremldat$casesToDrop

#Obtain phenotypic variance & covariance for start values:
vpi <- var(gremldat$yX[,1],na.rm=T)
cpi <- cov(phenfile[,3:4],use="pair")[1,2]

# Create MxData object, expectation, and compute plan: ##########################################

#The MxData object.  We will use gremldat$yX as the raw dataset to be used in our MxModel:
mxdat <- mxData(observed=gremldat$yX,type="raw",sort=FALSE)
#^^^N.B. the 'sort=FALSE' argument! 

#We tell the GREML expectation that the name of the model-expected covariance matrix is "V", and that the dataset
#being input is y horizontally adhered to X (and therefore, it should not call the data-handler at runtime):
ge <- mxExpectationGREML(V="V",dataset.is.yX=TRUE,casesToDropFromV=casesToDrop)

#Custom compute plan, to use Newton-Raphson and see what the optimizer is doing in real time (via the 
# 'verbose' argument):
plan <- mxComputeSequence(
	steps=list(
		mxComputeNewtonRaphson(verbose=5L),
		#mxComputeGradientDescent(engine="NPSOL",verbose=5L),
		mxComputeOnce('fitfunction', c('gradient','hessian')),
		mxComputeStandardError(),
		mxComputeHessianQuality(),
		mxComputeReportDeriv(),
		mxComputeReportExpectation()
	))
#^^^Note:  If you are running the R GUI under Windows, delete the 'verbose=5L' argument in the above.



# Create the MxModel: ################################################################################
#(We will create a lot of objects inside the mxModel() statement.  This helps to save memory.)

gremlmod <- mxModel(
	"Diphenotype_Continuous",
	mxdat, #<--MxData object.
	ge, #<--GREML expectation.
	plan, #<--Custom compute plan.
	
	# 3 MxMatrices to hold additive-genetic parameters:
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=vpi/2,labels="va1",name="Va1"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=vpi/2,labels="va2",name="Va2"),
	mxMatrix(
		type="Full",nrow=1,ncol=1,free=T,values=cpi/2,labels="ca",name="Ca"),
	# 3 MxMatrices to hold non-shared environmental parametes;
	# note that each is an order-N column vector:
	mxMatrix(type="Full",nrow=N,ncol=1,free=T,values=vpi/2,labels="ve1",name="Ve1"),
	mxMatrix(type="Full",nrow=N,ncol=1,free=T,values=vpi/2,labels="ve2",name="Ve2"),
	mxMatrix(
		type="Full",nrow=N,ncol=1,free=T,values=cpi/2,labels="ce",name="Ce"),
	
	#MxAlgebras for the two traits' heritabilities, and their genetic correlation:
	mxAlgebra(Va1/(Va1+Ve1[1,1]), name="h2_1"),
	mxAlgebra(Va2/(Va2+Ve2[1,1]), name="h2_2"),
	mxAlgebra(Ca/sqrt(Va1)/sqrt(Va2), name="rg"),
	
	#The GRM:
	mxMatrix(type="Symm",nrow=N,free=F,values=GRM,name="A"),
	#An order-N column vector of ones, to be used in computing the derivatives of V:
	mxMatrix(type="Unit",nrow=N,ncol=1,name="Uno"),
	#An NxN matrix of zeroes, to be used in computing the derivatives of V:
	mxMatrix(type="Zero",nrow=N,ncol=N,name="Zip"),
	
	#The model-expected covariance matrix, 'V':
	mxAlgebra(
		rbind(
			cbind(A%x%Va1, A%x%Ca),
			cbind(A%x%Ca, A%x%Va2)
		)
		+ rbind(
			cbind(vec2diag(Ve1), vec2diag(Ce)),
			cbind(vec2diag(Ce), vec2diag(Ve2))
		), name="V"),
	
	#MxAlgebras representing the first partial derivatives of 'V' w/r/t parameters:
	mxAlgebra(
		rbind(
			cbind(A,Zip),
			cbind(Zip,Zip)
		), name="dV_dva1"),
	mxAlgebra(
		rbind(
			cbind(Zip,A),
			cbind(A,Zip)
		), name="dV_dca"),
	mxAlgebra(
		rbind(
			cbind(Zip,Zip),
			cbind(Zip,A)
		), name="dV_dva2"),
	mxAlgebra(
		rbind(
			cbind(vec2diag(Uno),Zip),
			cbind(Zip,Zip)
		), name="dV_dve1"),
	mxAlgebra(
		rbind(
			cbind(Zip,vec2diag(Uno)),
			cbind(vec2diag(Uno),Zip)
		), name="dV_dce"),
	mxAlgebra(
		rbind(
			cbind(Zip,Zip),
			cbind(Zip,vec2diag(Uno))
		), name="dV_dve2"),
	
	#GREML fitfunction:
	mxFitFunctionGREML(
		dV=c(va1="dV_dva1",ca="dV_dca",va2="dV_dva2",ve1="dV_dve1",ce="dV_dce",ve2="dV_dve2"))
)

gremlmod <- mxRun(gremlmod)
( gremlsumm <- summary(gremlmod) )

#Save phenotype file to disk in format GCTA can read:
write.table(phenfile,"phenfile.dat",row.names=F,col.names=F)

#Invoke GCTA from inside R:
cmd <- system(
	intern=F,
	command="gcta --thread-num 2 --reml-bivar 1 2 --reml-no-constrain --reml-bivar-no-constrain --pheno phenfile.dat --grm DiphenoGRM --out diph")

if(cmd==0){ #<--If the above call to system() did not end in an error.
	#Read in GCTA output file:
	hsq <- read.table("diph.hsq",header=T,nrows=11)
	#OpenMx & GCTA point estimates of variance components are close:
	coef(gremlmod); hsq[1:6,2]
	#OpenMx & GCTA standard errors of variance components are close:
	gremlmod$output$standardErrors; hsq[1:6,3]
	#OpenMx's & GCTA's heritability estimates and their SEs are close:
	c(mxEval(h2_1,gremlmod,T), mxSE(h2_1,gremlmod)); hsq[9,]
	c(mxEval(h2_2,gremlmod,T), mxSE(h2_2,gremlmod)); hsq[10,]
	#OpenMx's & GCTA's genetic correlation estimates SEs are close:
	c(mxEval(rg,gremlmod,T), mxSE(rg,gremlmod)); hsq[11,]
	
	#The mxGREML fitfunction reports -2 times the natural logarithm of the normalized likelihood:
	gremlmod$output$fit
	#To compare OpenMx's & GCTA's -2logLs at the solution, we need to add the normalizing constant to the value GCTA
	#reports, -5449.22, and then multiply by -2:
	-2*(-5449.22 - 0.5*gremlsumm$observedStatistics*log(2*pi))
	#^^^They match.
}
