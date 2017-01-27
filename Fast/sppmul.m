function y = sppmul(S,v,x)
%SPPMUL Sparse pattern multiply.
%   Y = SPPMUL(S,V,X), where S is a sparse matrix and X is a full
%   matrix, computes the product Y = S*X but with the nonzero elements
%   of S replaced by those in the vector V. SPPMUL is thus formally
%   equivalent to the code
%      [i,j,s] = find(S);
%      Y = sparse(i,j,v,size(S,1),size(S,2))*X;
%   Of course, the intermediate sparse matrix is never actually 
%   constructed. It is required that NNZ(S) = PROD(SIZE(V)) and
%   that SIZE(S,2) = SIZE(X,1).

% S. Engblom 2005-10-27

error('.MEX-file not found on path.');
