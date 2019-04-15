### Stokes eigenvalues

This directory contains the code used to generate
the figures in the paper "A boundary integral equation 
approach to computing eigenvalues of the Stokes operator".
This repository is meant to be a snapshot of the
code as used in the paper.
The software development for this paper led to some
useful code which has been pulled out and will
be maintained more professionally under a different
name (more to come on this!)

The primary files of interest are located in matlab/example
but compiling the mexfiles is necessary.

### Compiling the mexfiles


From this top-level folder, run

```
make mexfiles
```
to generate the mexfiles used by the MATLAB
routines. You may have to edit the flags for
compiling. The current version requires gfortran
and mex.

### Running the examples

The examples were run with MATLAB R2018a.

See instructions in the matlab/example folder


### Licenses

The software in this repository contained
in the src, mwrap, and matlab directories
(except for the matlab/external subdirectory)
is available under the terms of the FreeBSD
3 clause license unless otherwise noted
(see LICENSE.md) for details. Many of the files
src folder are from FMMLIB2D, which is subject to
the terms of src/COPYING. 

The software in the "external" folders is available
under different licenses and is housed here for archiving
purposes.
