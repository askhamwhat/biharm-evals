### Biharmonic eigenvalue problem

This directory is early days. 

The test values for the difference kernel functions 
can be generated from the files mathematica/HelmBHTest.nb
and mathematica/HelmBHDerTests.nb. They are already included
in tests/test-data

### Running tests

The test directory contains makefiles for running 
many tests. Most accept a command line argument
for the system (either linux or mac), e.g.

make testhbhstokesgreenid SYSTEM=mac

will compile with the correct flags for mac.

