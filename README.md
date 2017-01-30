# stenglib
Stefan Engblom's Matlab libraries

##Tensor##

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
  [tndims.m](Tensor/tensor.m)
  [tndims.c](Tensor/source/tensor.c)

* **tsize** Size of array.
  (Depend: 0, status: stable)
  [tsize.m](Tensor/tsize.m)
  [tsize.c](Tensor/source/tsize.c)

* **tsum** Tensor summation.
  (Depend: 0, status: stable)
  [tsum.m](Tensor/tsum.m)
  [tsum.c](Tensor/source/tsum.c)

* **tprod** Tensor product. Based on a concept by
      D. Bertilsson, [COMSOL](http://www.comsol.com).
  (Depend: 0, status: stable)
  [tprod.m](Tensor/tprod.m)
  [tprod.c](Tensor/source/tprod.c)
  
There is also a [make.m](Tensor/source/make.m) available.
It will work on several, but not all, platforms.

##Fast##

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
in *Parallel Comput.* 56:1--17 (2016) [(doi)](http://dx.doi.org/10.1016/j.parco.2016.04.001). *Fact:*
the **fsparse**-code has been selected as the base for the sparse assembly routines
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
