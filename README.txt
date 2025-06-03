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

In early June of 2025, the scripts in the repository were updated to reflect 
the current recommended best practices in doing mxGREML analyses.  The current
advice, concerning several topics in alphabetical order, is as follows:
* CSOLNP:  CSOLNP should generally be the optimizer of first choice for cases 
where Newton-Raphson (q.v.) is not recommended, though there are a few purposes
where SLSQP (q.v.) might be preferred over CSOLNP.  Compared to SLSQP, CSOLNP
tends to be "safer", whereas SLSQP tends to be faster.  That is, CSOLNP can
generally be more trusted than SLSQP to reach a good solution in 1 fit attempt,
but when both optimizers reach substantially equivalent solutions, SLSQP
typically takes less time to get there.
* `infoMatType`:  `infoMatType` is an argument to `mxFitFunctionGREML()`.  It 
dictates whether the fitfunction will use the average-information or expected-
information matrix as an approximation to the Hessian.  Its default, 
`infoMatType="average"`, is recommended only when the model-expected covariance 
matrix 'V' is linear in the free parameters, or equivalently, when the 
derivatives of 'V' are constant with respect to the free parameters.  If that
condition is not met, then `infoMatType="expected"` is recommended because it 
will be more accurate (though if time is of the essence, 
`infoMatType="average"` will be faster, and may be "accurate enough" for some 
users).
* Newton-Raphson (N-R):  In cases where N-R works well, it is very fast and very
accurate.  N-R should generally be the optimizer of first choice when the
model-expected covariance matrix 'V' is linear in the free parameters,
or equivalently, when the derivatives of 'V' are constant with respect to the
free parameters.  Even if that condition is not met, N-R can sometimes be
useful nonetheless, but it often is not.  On the other hand, even if that
condition IS met, N-R can fail if, say, one or more parameters has an active 
`lbound` or `ubound` near the solution.  N.B. that, as of OpenMx v2.22.7, N-R
cannot be used with semi-analytic derivatives (i.e., it can only be used with
fully analytic derivatives) due to a bug.  That bug has been repaired in the
OpenMx source repository, and will be repaired in the next OpenMx release.
N-R is incompatible with MxConstraints.
* NPSOL:  NPSOL has repeatedly exhibited pathological behavior when optimizing
mxGREML models.  Specifically, it sometimes appears to get stuck in an infinite
loop, within which it does not make any discernible progress toward a solution.
Consequently, NPSOL is no longer recommended for routine use with mxGREML
models, and as of this writing, this repository contains no demo scripts that
use it.  NPSOL remains potentially useful in mxGREML analyses for two niche
purposes: warm starts (q.v.), and verifying the correctness of the user's
analytic derivatives.
* SLSQP:  For unclear reasons, SLSQP sometimes stops its optimization of an
mxGREML model "too early" (i.e., when it has not reached a local minimum), but
neither it nor OpenMx's interface to it detects this optimization failure.
Unless and until this potential bug is resolved, CSOLNP (q.v.) is recommended 
over SLSQP for most mxGREML use.  However, there are at least 3 particular
circumstances where SLSQP is justifiably be preferred over CSOLNP: 
(1) optimization of profile-likelihood confidence intervals (i.e., calls to
`mxCI()`), (2) optmization involving MxConstraints, and (3) optimization using 
only numeric derivatives (which is slow, and not generally advised).
* Warm starts:  If the Hessian is positive-definite at the start values, then
NPSOL can accept the Hessian's Cholesky factorization as a "warm start", and
thereby begin optimization with information about the fitfunction's local 
curvature.  With mxGREML models, providing NPSOL with a warm start is only
recommended when the model-expected covariance matrix 'V' is linear in the free
parameters (or equivalently, when the derivatives of 'V' are constant with 
respect to the free parameters) and default argument `infoMatType="average"`
is in use.

Reference:
Kirkpatrick RM, Pritikin JN, Hunter MD, & Neale MC.  (2021).  Combining
structural-equation modeling with genomic-relatedness-matrix restricted maximum
likelihood in OpenMx.  Behavior Genetics, 51(3): 331-342.
doi: 10.1007/s10519-020-10037-5
