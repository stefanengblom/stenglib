function c = tsum(a,b,ia,ib)
%TSUM Tensor summation.
%   C = TSUM(A,B,IA,IB) computes a general sum of the arrays A and B in
%   the order specified by the vectors IA and IB.
%
%   This syntax is best explained by an example. Let A be a 3-dimensional
%   array and let B be a 2-dimensional array. Then
%     C = TSUM(A,B,[1 3 2],[3 1])
%   creates a 3-dimensional array C by forming the sum
%     C(i1,i2,i3) = A(i1,i3,i2)+B(i3,i1).
%   Of course, it is assumed that SIZE(A,1) = SIZE(B,2) and SIZE(A,2) =
%   SIZE(B,1). The output C has the dimensions
%   SIZE(A,1)-by-size(A,3)-by-SIZE(A,2).
%
%   The arguments are assumed to have the following format:
%     - A, B, and C are (real or complex) arrays of doubles.
%
%     - IA and IB are vectors of doubles, containing positive
%     integers. The length of IA(IB) has to be equal to the number of
%     dimensions of A(B). However, since MATLAB automatically removes
%     trailing singleton dimensions, A is padded with singleton dimensions
%     if the number of dimensions of A is less than the length of IA (and
%     similarly for B). Unlike MATLAB, TSUM supports 0-D and 1-D tensors
%     (scalars and column vectors), so the index vectors may contain less
%     than 2 elements.
%
%     - The numbers in IA(IB) have to be distinct.
%
%     - If a number occurs both in IA and IB, it is required that the
%     corresponding dimensions in A and B have the same size.
%
%     - It is assumed that the union of the numbers in IA and IB form a
%     contiguous sequence of integers.
%
%   For optimal performance of TSUM you should try to have the vectors IA,
%   IB and [IA IB] in order.
%
%   C = TSUM(A,IA) sums the array A along the dimensions specified in the
%   ordered vector IA. This is an extension of the MATLAB-function SUM,
%   which only allows for summation along one dimension. After summing,
%   the dimensions specified in IA become singleton dimensions of C; you
%   may efficiently remove them by using the MATLAB-function SQUEEZE.
%
%   As a further extension of MATLABs SUM-function, negative numbers in IA
%   may be used to force backward summation of the corresponding
%   dimensions.
%
%   Examples:
%     % the "outer sum" of the vectors v and w
%     v = rand(3,1); w = rand(4,1);
%     C1 = tsum(v,w,[1],[2]);
%
%     % the symmetric part of the matrix A
%     A = rand(4);
%     C2 = 0.5*tsum(A,A,[1 2],[2 1]);
%
%     % the transpose of the matrix A
%     C3 = tsum(A,0,[2 1],[]);
%
%     % add i to row i of A
%     C4 = tsum(A,1:4,[1 2],[3 1]);
%
%     % sum dimensions 1 and 3 of B
%     B = rand(3,4,5);
%     C5 = tsum(B,[1 3]);
%
%     % sum backwards
%     C6 = tsum([1:100000].^-4,-2);
%
%   See also TPROD, SUM, SQUEEZE.

% S. Engblom 2005-08-26

error('.MEX-file not found on path.');
