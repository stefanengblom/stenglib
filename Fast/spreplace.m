function T = spreplace(S,x)
%SPREPLACE Sparse pattern replace.
%   T = SPREPLACE(S,X), where S is a sparse matrix and X a vector,
%   replaces the nonzero elements of S with X. SPREPLACE is thus
%   formally equivalent to the code
%      [i,j,s] = find(S);
%      T = sparse(i,j,x,size(S,1),size(S,2));
%   It is required that NNZ(S) = PROD(SIZE(X)).
%
%   Example:
%      S = sprand(30,30,0.2);
%      F = spreplace(S,exp(nonzeros(S)));
%
%   See also SPPMUL, FSPARSE.

% S. Engblom 2016-11-23

error('.MEX-file not found on path.');
