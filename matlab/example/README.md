# examples for Stokes eigenvalue paper

In this folder are the files used to generate the
figures in the paper "A boundary integral equation 
approach to computing eigenvalues of the Stokes 
operator".

## Running the examples

First, you must build the mexfiles used by these
routines by doing

```
make mexfiles
```
in the top-level directory (two levels up).

These routines also rely on the FLAM library
available at https://github.com/klho/FLAM
and the chebfun library available at
http://www.chebfun.org/. 
(We include a local copy of the versions used 
to generate the figures in the top-level directory
under "extensions".)