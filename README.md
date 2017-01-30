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
