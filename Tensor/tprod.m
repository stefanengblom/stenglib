function c = tprod(a,b,ia,ib)
%TPROD Tensor product.
%   C = TPROD(A,B,IA,IB) computes a general tensor product of the arrays A
%   and B. The supported tensor product can be described as a general
%   multiplication of the elements in A and B where some indices may be
%   equal and/or summed over. The mapping from input indices to output
%   indices, as well as how to sum, is described by the vectors IA and IB.
%
%   This function is best explained by an example. Let A be a
%   4-dimensional array and let B be a 3-dimensional array. Then
%     C = TPROD(A,B,[2 -1 1 -2],[-2 2 -1])
%   creates a 2-dimensional array (matrix) C in the following way.
%   First, the product
%     D(j2,j1,i1,i2) = A(i2,j1,i1,j2)*B(j2,i2,j1)
%   is formed. This is a 4-dimensional array D, where j1, j2, i1, i2
%   denote index variables. D is the tensor (outer) product of A and B,
%   followed by a permutation of the indices and setting some indices
%   equal. Of course, it is assumed that SIZE(A,1) = SIZE(B,2), SIZE(A,2)
%   = SIZE(B,3), and SIZE(A,4) = SIZE(B,1). Second, we sum over the index
%   variables corresponding to negative numbers (j1 and j2):
%     C(i1,i2) = sum of D(j2,j1,i1,i2)
%   where the indices j1 and j2 (independently) run through all their
%   possible values.
%
%   The arguments are assumed to have the following format:
%     - A, B, and C are (real or complex) arrays of doubles.
%
%     - IA and IB are vectors of doubles, containing nonzero integers. The
%     length of IA(IB) has to be equal to the number of dimensions of
%     A(B). However, since MATLAB automatically removes trailing singleton
%     dimensions, A is padded with singleton dimensions if the number of
%     dimensions of A is less than the length of IA (and similarly for B).
%     Unlike MATLAB, TPROD supports 0-D and 1-D tensors (scalars and
%     column vectors), so the index vectors may contain less than 2
%     elements.
%
%     - The numbers in IA(IB) have to be distinct.
%
%     - If a number occurs both in IA and IB, it is required that the
%     corresponding dimensions in A and B have the same size.
%
%     - If a negative number occurs in A(B), it must also occur in B(A).
%
%     - It is assumed that the union of the numbers in IA and IB together
%     with 0 form a contiguous sequence of integers.
%
%   For optimal performance of TPROD you should try to have the vectors
%   IA, IB and [IA IB] in order.
%
%   Examples:
%     % the tensor (outer) product of the matrices A and B
%     A = rand(3,4); B = rand(2,5);
%     C1 = tprod(A,B,[1 2],[3 4]);
%
%     % the Kronecker product of the matrices A and B
%     C2 = reshape(tprod(A,B,[2 4],[1 3]),size(A).*size(B));
%
%     % the ordinary matrix product of the matrices A and B
%     A = rand(3,4); B = rand(4,2);
%     C3 = tprod(A,B,[1 -1],[-1 2]);
%
%     % the element-wise product of the matrices A and B
%     A = rand(3,5); B = rand(3,5);
%     C4 = tprod(A,B,[1 2],[1 2]);
%
%     % the transpose of the matrix A
%     C5 = tprod(A,1,[2 1],[]);
%
%     % the sum of all entries in the matrix A
%     C6 = tprod(A,ones(size(A)),[-1 -2],[-1 -2]);
%
%     % the sum of all diagonal entries in the matrix A (the trace)
%     C7 = tprod(A,eye(size(A)),[-1 -2],[-1 -2]);
%
%   See also KRON, MTIMES, TIMES.

% S. Engblom 2005-04-04
% Based on a concept by D. Bertilsson, COMSOL AB.

error('.MEX-file not found on path.');
