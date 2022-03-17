# FortranFilterSLP
FORTRAN implementation of FilterSLP

`original`
---------

This version dates from (at the latest) 4 July 2013. 
* It won't compile with modern `f90` (specifically `gcc version 9.4.0 - Ubuntu 9.4.0-1ubuntu1~20.04`) due to some FORTRAN usage that's now illegal.

Specifically, two logical statements `/=` in lines 224 and 229 of `rd_prob_da.f` can no longer be compiled

`master`
-------

This includes 
* Edits to allow it to be compiled with modern `f90`. Specifically the uses of `/=` have been replaced by `.neqv.`

As such, this represents the "best" version of the FORTRAN implementation of FilterSLP that @jajhall can come up with

Note
* It calls EMSOL (see https://github.com/ERGO-Code/EMSOL)
* To compile may require the`-fno-range-check` flag to cope with the mechanism used in `EMSV.INC` to set the 32nd bit in a 32-bit integer
* No version of the FORTRAN implementation of FilterSLP is supported in any way
* Anyone using a FORTRAN implementation of FilterSLP is strongly advised to continue using the version that they have!

