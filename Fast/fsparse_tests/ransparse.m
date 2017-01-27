function [ii,jj,ss,siz] = ransparse(siz,nnz_row,nrep)
%RANSPARSE Test input for sparse matrix assembly.
%   [I,J,S,siz] = RANSPARSE(siz,NNZ_ROW,NREP) creates randomized input
%   data (I,J,S,siz) for assembly via SPARSE(I,J,S,siz,siz). The
%   resulting sparse matrix has approximately NNZ_ROW nonzeros per row
%   and each nonzero entry results from about NREP collisions.
%
%   Example:
%     [i,j,s,siz] = ransparse(100,15,8);
%     S = sparse(i,j,s,siz,siz);
%
%     nnz_row = mean(sum(S > 0,2))
%     [foo,bar,val] = find(S);
%     nrep = mean(val)

rng('default')
rng(123)           % init/seed the generator

nnz = nnz_row*siz; % number of nonzeros
ii = repmat((1:siz)',[1 nnz_row]);
jj = ceil(rand(siz,nnz_row)*siz);
ii = repmat(ii(:),[1 nrep]);
jj = repmat(jj(:),[1 nrep]); % (some jj's might be the same)

p = randperm(numel(ii));
ii = ii(p);
jj = jj(p);
ss = ones(size(ii));
