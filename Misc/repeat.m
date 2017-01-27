function c = repeat(a,b)
%REPEAT Repeat elements of vector.
%   C = REPEAT(A,B) returns A repeated by B, see the example
%   below. Both arguments are assumed to be row vectors of the same
%   length and B must contain nonnegative integers only.
%
%   No error-checking is performed.
%
%   Example:
%     a = [0 1 pi inf nan -nan];
%     b = [1 4 0   2   3    0 ];
%     c = repeat(a,b) % returns c = [0 1 1 1 1 inf inf nan nan nan]

% S. Engblom 2005-10-14

% can't beat this!
ix = cumsum(sparse(1,cumsum([1 b]),1));
c = a(ix(1:end-1));
