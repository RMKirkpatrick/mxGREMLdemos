Copyright 2019-2021 by Robert M. Kirkpatrick
Licensed under CC BY 4.0 <http://creativecommons.org/licenses/by/4.0/>

This repository is a collection of R scripts that demonstrate use of the 
mxGREML feature in the R package 'OpenMx'.  The scripts in subdirectory 
`fromTestSuite` are adapted from scripts in the OpenMx project's nightly test
suite.  The scripts in subdirectory `AGES2017` are adapted from example 
scripts used at the 2017 AGES Workshop at Virginia Commonwealth University.
Most scripts use MxAlgebras that represent the first partial derivative of
the model-expected covariance matrix, which in turn the OpenMx backend
can use to calculate the first and second partial derivatives of the
GREML fitfunction.  This use of analytic derivatives is best practice in terms
computational speed.  In contrast, scripts with 'NodV' in their filename use
SLSQP to calculate numeric derivatives in a multithreaded fashion, and do not
use analytic derivatives.

Newcomers to the mxGREML feature are encouraged to first examine the scripts
in the `AGES2017` subdirectory, since those scripts are the most didactically
oriented in the repository.
