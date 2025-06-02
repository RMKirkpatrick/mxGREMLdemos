Copyright 2019-2025 by Robert M. Kirkpatrick
Licensed under CC BY 4.0 <http://creativecommons.org/licenses/by/4.0/>

This repository is a collection of R scripts that demonstrate use of the 
mxGREML feature in the R package 'OpenMx'.  The scripts in subdirectory 
`fromTestSuite` are adapted from scripts in the OpenMx project's nightly test
suite.  The scripts in subdirectory `AGES2017` are adapted from example 
scripts used at the 2017 AGES Workshop at Virginia Commonwealth University.

Some scripts use MxAlgebras that represent the first partial derivative of
the model-expected covariance matrix, which in turn the OpenMx backend
can use to calculate the first and second partial derivatives of the
GREML fitfunction.  Other scripts use "semi-analytic" derivatives
(as described in Kirkpatrick et al., 2021), meaning that the OpenMx backend
automatically calculates numerical derivatives of the model-expected
covariance matrix, which are then used to analytically calculate derivatives
of the GREML fitfunction.  If the repository contains two versions of a
script--one using fully analytic derivatives and the other using semi-analytic
derivatives--then the version using fully analytic derivatives has 'dV' in its
filename.

Newcomers to the mxGREML feature are encouraged to first examine the scripts
in the `AGES2017` subdirectory, since those scripts are the most didactically
oriented in the repository.

Reference:
Kirkpatrick RM, Pritikin JN, Hunter MD, & Neale MC.  (2021).  Combining
structural-equation modeling with genomic-relatedness-matrix restricted maximum
likelihood in OpenMx.  Behavior Genetics, 51(3): 331-342.
doi: 10.1007/s10519-020-10037-5
