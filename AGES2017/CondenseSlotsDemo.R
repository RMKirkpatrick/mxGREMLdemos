# Copyright 2019-2022 by Robert M. Kirkpatrick
# Licensed under CC BY 4.0 <http://creativecommons.org/licenses/by/4.0/>

# This script demonstrates how the "condensed matrix slots" feature works.

require(OpenMx)
options(mxCondenseMatrixSlots=FALSE) #<--Set this option to its on-load default...

#Consider the following MxMatrix.  It has no free parameters...:
Fu <- mxMatrix(type = "Full", nrow = 2, ncol=3, free = FALSE, values = 1:6, name="Fu",
							 dimnames = list(c("r1","r2"),c("c1","c2","c3")))
str(Fu) #<--We can get a summary of its structure with the str() command.
Fu@free #<--Access the literal value of its 'free' slot.  There are no free parameters, so it's a matrix of FALSE.
Fu@labels #<--Access the literal value of its 'labels' slot.  All NA.
Fu@lbound #<--Access the literal value its 'lbound' slot.  All NA.
Fu@ubound #<--Access the literal value its 'ubound' slot.  All NA.

#Now, imagine we want to work with a 10,000 x 10,000 MxMatrix that contains no free parameters.  That MxMatrix object
#would contain a 10,000 x 10,000 matrix of FALSE, and three 10,000 x 10,000 matrices of NA.
#That's a waste of memory!

#Function mxMatrix() has an argument 'condenseSlots' (which defaults to the value of global 
#option 'mxCondenseMatrixSlots').  Let's change the value of that option and reconstruct Fu:
options(mxCondenseMatrixSlots=TRUE)
Fu <- mxMatrix(type = "Full", nrow = 2, ncol=3, free = FALSE, values = 1:6, name="Fu",
							 dimnames = list(c("r1","r2"),c("c1","c2","c3")))
#Notice that now, slot 'free' is a 1x1 matrix of FALSE, and 'labels', 'lbound', and 'ubound' are 1x1 matrices of NA:
str(Fu)
#We can access the literal values of those slots and see:
Fu@free
Fu@labels
Fu@lbound
Fu@ubound
#BUT, notice what happens if we access those slots using $ instead of @:
Fu$free
Fu$labels
Fu$lbound
Fu$ubound
#^^^The matrices in those four slots are automatically decondensed to full size and given their dimension names.

#We can override the global option 'mxCondenseMatrixSlots' by supplying a value to argument 'condenseSlots' in 
#mxMatrix():
Fu <- mxMatrix(type = "Full", nrow = 2, ncol=3, free = FALSE, values = 1:6, name="Fu", condenseSlots=FALSE,
							 dimnames = list(c("r1","r2"),c("c1","c2","c3")))
str(Fu) #<--And those four slots contain full-size matrices once again.

#If we modify the appropriate slot WITH THE $ ACCESSOR, we can toggle condensed slots on and off (THIS WILL NOT 
#WORK WITH THE @ ACCESSOR):
Fu$.condenseSlots <- TRUE
str(Fu)
Fu$.condenseSlots <- FALSE
str(Fu)
Fu$.condenseSlots <- TRUE
str(Fu)

#Using the $ accessor to give the MxMatrix either a labeled element or free-parameter element prevents it from
#condensing the 'labels' and 'free' matrices:
Fu$free[2,2] <- TRUE
str(Fu)
Fu$free[2,2] <- FALSE
str(Fu)
Fu$free[2,2] <- TRUE
str(Fu)
Fu$free[2,2] <- FALSE
Fu$labels[2,2] <- "foo"
str(Fu)
Fu$labels[2,2] <- NA
str(Fu)
Fu$labels[2,2] <- "foo"
str(Fu)
Fu$labels[2,2] <- NA

