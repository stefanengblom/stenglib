# stenglib
Stefan Engblom's Matlab libraries - packages for daily use.

## License statement for stenglib

You may download all of **stenglib** and
use, modify and redistribute it in any way you like. A
redistributor must fully attribute the authorship and make a good
effort to cite the original location of the software. A researcher
making *critical* use of the software in research is requested to 
acknowledge this in publications related to the
research. A company may use the code in software products provided
that the original location and the author is clearly cited.

All code provided here comes with absolutely **no warranty**
and **no support** whatsoever is given.

There are a lot of freeware available on the net. Do **not**
download unless you agree to the above license.

## About stenglib

**stenglib** is loosely divided into 5 sub-packages, with few dependencies in between them:
* Tensor
* Fast
* Scicomp
* Utils
* Misc

For contact details, see [stenglib.org](http://www.stenglib.org).
I welcome bug reports and comments. Please *do not* send support requests.

## Tensor

Originally, I made the **Tensor** package
because I had the need to easily, efficiently and consistently
manage multi-dimensional arrays in Matlab. *Examples:* given a
matrix and a vector, how do you scale each row in the matrix by
the vector? How can you multiply a 3-D array with a matrix? The
package is useful to anyone who writes code for (pseudo-) spectral
methods, FEM, or who uses multi-dimensional arrays or tensor
notation a lot.

* **tndims** Number of dimensions.
  (Depend: 0, status: stable)
  [tndims.m](Tensor/tndims.m)
  [tndims.c](Tensor/source/tndims.c)

* **tsize** Size of array.
  (Depend: 0, status: stable)
  [tsize.m](Tensor/tsize.m)
  [tsize.c](Tensor/source/tsize.c)

* **tsum** Tensor summation.
  (Depend: 0, status: stable)
  [tsum.m](Tensor/tsum.m)
  [tsum.c](Tensor/source/tsum.c)

* **tprod** Tensor product.
  Based on a concept by D. Bertilsson, [COMSOL](http://www.comsol.com).
  (Depend: 0, status: stable)
  [tprod.m](Tensor/tprod.m)
  [tprod.c](Tensor/source/tprod.c)
  
There is also a [make.m](Tensor/source/make.m) available.
It will work on several, but not all, platforms.

## Fast

The routines in the **Fast** package exist
because some things just take too much time in Matlab. *Examples:*
replicate a data set in different dimensions (a.k.a. repmat),
assemble a sparse matrix, or evaluate set operations. These
routines should be of general interest to programmer in the
scientific computing community.

* **frepmat** Fast replication of array.
  (Depend: Tensor/{tndims,tsize} (weakly), status: stable)
  [frepmat.m](Fast/frepmat.m)
  [mexfrepmat.c](Fast/sourcee/mexfrepmat.c)

* **fsparse** Fast assembly of sparse matrix.
  (Depend: 0, status: stable but not completely settled) 
  [fsparse.m](Fast/fsparse.m)
  [fsparse.c](Fast/source/fsparse.c)
  
There is now a *parallel* **fsparse** version available. A
paper describing the algorithm is *S. Engblom, D. Lukarski:
Fast Matlab compatible sparse assembly on multicore computers*,
in *Parallel Comput.* 56:1--17 (2016) [(doi)](http://dx.doi.org/10.1016/j.parco.2016.04.001).
*Fact:* the **fsparse**-code has been selected as the base for the sparse assembly routines
in [PARALUTION](http://www.paralution.com).

* **clenshaw** Evaluation of 3-term recurrences.
  (Depend: 0, status: stable)
  [clenshaw.m](Fast/clenshaw.m)
  [clenshaw.c](Fast/source/clenshaw.c)

* **fsetop** Fast set operations based on hashing. Based on
  a concept by [P.-O. Persson](http://www.mit.edu/~persson) and a
  hash-function by [P. Hsieh](http://www.azillionmonkeys.com/qed/hash.html)
  (Depend: 0, status: stable)
  [fsetop.m](Fast/fsetop.m)
  [fsetop.c](Fast/source/fsetop.c)
  
* **sppmul** Sparse pattern multiply.
  (Depend: 0, status: stable) Download:
  [sppmul.m](Fast/sppmul.m)
  [sppmul.c](Fast/source/sppmul.c)

* **powerseries** Sum power series.
  (Depend: 0, status: stable)
  [powerseries.m](Fast/powerseries.m)
  [powerseries.c](Fast/source/powerseries.c)

As before there is a [make.m](Fast/source/make.m) available which you will probably have to modify.

## Scicomp

In **Scicomp** I've assembled some solvers from different areas within scientific computing: two solvers for
nonlinear problems and an implementation of the Nelder-Mead simplex algorithm. I have also put three routines 
for Gaussian quadratures with respect to discrete measures in this package.

* **rtsafe** Scalar nonlinear solver.
  (Depend: 0, status: stable)
  [rtsafe.m](Scicomp/rtsafe.m)
  
* **ainsolve** Solver for large sets of nonlinear equations.
  (Depend: 0, status: experimental)
  [ainsolve.m](Scicomp/ainsolve.m)

* **nmsimplex** Nelder-Meads simplex-algorithm.
  (Depend: Fast/frepmat, Utils/parseopts, status: stable)
  [nmsimplex.m](Scicomp/nmsimplex.m)

* **ode1s** ODE-solver: split-step Euler method with digital control.
  (Depend: Fast/{frepmat,fsetop}, Utils/{parseopts,report}, status: stable, but not finalized)
  [ode1s.m](Scicomp/ode1s.m)

* **gaussqd** Gaussian quadrature for discrete measures.
  (Depend: Tensor/tprod, Fast/fsparse, status: stable)
  [gaussqd.m](Scicomp/gaussqd.m)

* **gausspd** Gaussian polynomial for discrete measures.
  (Depend: Tensor/{tsum,tprod}, Fast/clenshaw, status: stable)
  [gausspd.m](Scicomp/gausspd.m)
  
* **poch** Pochhammer's function.
  (Depend: 0, status: stable)
  [poch.m](Scicomp/poch.m)

## Utils

In the package **Utils** I've collected various
routines for performing everyday tasks. Examples include
generating LaTeX-arrays from matrices, .gif-animations, small
perfect hash functions and removing files ending with a tilde (!).

* **arr2latex** LaTeX-table from matrix.
(Depend: Fast/frepmat, Utils/parseopts, status: stable)
[arr2latex.m](Utils/arr2latex.m)

* **assignopts** Assign options.
(Depend: 0, status: stable)
[assignopts.m](Utils/assignopts.m)

* **parseopts** Parse options.
(Depend: Fast/fsetop, status: stable)
[parseopts.m](Utils/parseopts.m)

* **matmerge** Merge .mat-files.
(Depend: Fast/fsetop, Tensor/{tsize,tndims}, status: stable but kind of a hack)
[matmerge.m](Utils/matmerge.m)

* **consistency** Local truncation error of linear multistep method.
(Depend: 0, status: stable)
[consistency.m](Utils/consistency.m)

* **stability** Plot of stability region for linear multistep method.
(Depend: 0, status: stable)
[stability.m](Utils/stability.m)

* **perfecthash** Perfect hash-function from strings.
(Depend: Tensor/tsum, Fast/fsparse, status: this is a hack, really)
[perfecthash.m](Utils/perfecthash.m)

* **movie2gif** Create GIF animation from MATLAB movie.
(Depend: imwrite-function in Matlab, status: stable)
[movie2gif.m](Utils/movie2gif.m)

* **report** Report progress of solver.
(Depend: 0, status: stable)
[report.m](Utils/report.m)

* **runtest** General test facility.
(Depend: 0, status: stable)
[runtest.m](Utils/runtest.m)

* **rmtilde** Remove files ending with a tilde ('~').
(Depend: 0, status: stable)
[rmtilde.m](Utils/rmtilde.m)

* **connect** Connectivity information for points.
(Depend: 0, status: stable)
[connect.m](Utils/connect.m)

* **ndop** General N-dimensional operator.
(Depend: Tensor/{tndims,tsum,tprod}, Fast/frepmat, status: stable)
[ndop.m](Utils/ndop.m)

* **spblock** Sparse matrix from blocks.
(Depend: Fast/fsparse, status: stable)
[spblock.m](Utils/spblock.m)

* **spym** Visualize magnitude of elements in matrix.
(Depend: 0, status: this is a hack)
[spym.m](Utils/spym.m)

## Misc

I'll put amusing routines in **Misc**. For now,
this humble package contains a sudoku solver and a short function
which I personally believe is the most beautiful Matlab-code ever
written.

* **sudsolve** Solve Sudoku puzzle.
(Depend: 0, status: stable)
[sudsolve.m](Misc/sudsolve.m)

* **repeat** Repeat elements of vector.
(Depend: 0, status: stable)
[repeat.m](Misc/repeat.m)
