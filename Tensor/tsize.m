%TSIZE Size of array.
%   SIZ = TSIZE(A) returns the size of the array A. Unlike the MATLAB
%   function SIZE, the size of SIZ itself is 1-by-TNDIMS(A) (see this
%   function), so that scalars and column vectors are treated in a
%   consistent way.
%
%   SIZ = TSIZE(A,DIMS) returns the size along the dimensions DIMS. Unlike
%   the corresponding syntax using SIZE, DIMS may be a vector.
%
%   Examples:
%     [tsize(rand(2,2)) size(rand(2,2))]
%     [tsize(rand(2,1)) size(rand(2,1))]
%     [tsize(rand(1,1)) size(rand(1,1))]
%
%     tsize(rand(4,3,2,1),1:4)
%
%   See also TNDIMS, SIZE.

% S. Engblom 2005-04-10

error('.MEX-file not found on path.');
