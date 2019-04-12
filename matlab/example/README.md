# examples for Stokes eigenvalue paper

In this folder are the files used to generate the
figures in the paper "A boundary integral equation 
approach to computing eigenvalues of the Stokes 
operator".

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

