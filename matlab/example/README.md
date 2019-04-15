# examples for Stokes eigenvalue paper

In this folder are the files used to generate the
figures in the paper "A boundary integral equation 
approach to computing eigenvalues of the Stokes 
operator".

## Requirements

The examples were run with MATLAB R2018a. We rely on
two external packages (included in this repo) as described
below. We also make use of the parallel computing
toolbox in MATLAB.

For the mexfiles, you will need a fortran and mex
compiler. Using a reasonably new gfortran version
(we used gcc v 5.4.0)
and the mex shipped with R2018a should work 

## Before running the examples

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

## Running the examples


### Barbell example

Run example_barbell_001.m then
example_barbell_001_big_plots.m,
which will take a few hours.

The post processing script is best run interactively
to see what information is produced for plotting.
See example_barbell_001_big_plots_post.m

### Many inclusions example

Run example_many_holes_004.m then
example_many_holes_004_big_plots.m.
Run example_many_holes_005.m
then example_many_holes_005_big_plots.m.
This will take a few hours.

The post processing script is best run interactively
to see what information is produced for plotting.
See example_many_holes_004_and_005_big_plots_post.m


### Annulus eigenvalue convergence studies
Run example_annulus_eig_conv.m.
This will take around 30 mins.

On output, figure 1 plots the error 
in convergence to an eigenvalue for
the combined field representation,
figure 2 plots the error in convergence 
to an eigenvalue for the double layer representation,
and figure 3 plots the error in convergence
to a spurious eigenvalue for the double layer
representation. 
In all the three plots, the dashed line indicates
the expected convergence rate of $N^{-20}$. 
Figure 4,5 plot the normalized determinant on the interval
[14,15] for combined field and double layer 
representation reprectively.

The results are also stored in mat-files/annulus_eig_conv.mat
which can be loaded for post processing.


### Annulus quadrature convergence studies

Run example_annulus_solv_conv.m. 
This will take a few minutes

On output, figure 1 plots the error in the vorticity
at a random point in the interior as a function of 
N, for the combined field representation 
and figure 2 plots the error for the double 
layer representation. 
The dashed line in both the plots indicates the 
expected convergence rate of $N^{-20}$.
The results are also stored as res-files/solv-conv-comb.pdf
and res-files/solv-conv-dl.pdf

The results are also stored in mat-files/annulus_solv_conv.mat
which can be loaded for post processing

### Annulus speed tests
Run example_annulus_speedtest.m, followed by
python plot_speed.py
This test will approximately take an hour.

The result is the time taken to compute the determinant
as a function of $N$ for three oscillatory stokes
parameters. The result is also stored as res-files/speed_res.pdf

The results are also stored in example_annulus_speedtest.mat
which can be loaded for post-processing


