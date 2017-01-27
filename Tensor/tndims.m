%TNDIMS Number of dimensions.
%   N = TNDIMS(A) returns the number of dimensions of the array A. TNDIMS
%   works just like the MATLAB function NDIMS except that scalars and
%   column vectors are detected correctly, see the examples below. This is
%   useful since the equivalent construction using NDIMS and SIZE is quite
%   ugly.
%
%   Examples:
%     [tndims(rand(2,2)) ndims(rand(2,2))]
%     [tndims(rand(2,1)) ndims(rand(2,1))]
%     [tndims(rand(1,1)) ndims(rand(1,1))]
%
%   See also NDIMS, SIZE, TSIZE.

% S. Engblom 2005-04-10

error('.MEX-file not found on path.');
