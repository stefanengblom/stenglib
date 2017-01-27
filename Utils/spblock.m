function S = spblock(i,j,v,mn)
%SPBLOCK Sparse matrix from blocks.
%   S = SPBLOCK(I,J,V,[M N]) assembles an M-by-N sparse matrix S from
%   the index vectors I and J and the cell vector V containing sparse
%   matrices. The assembly is performed in blocks so that for each
%   index k, the block V{k} is placed with the upper left corner in
%   (I(k),J(k)) in S. Multiple indices are as usual summed together,
%   see FSPARSE. Also, single index or blocks are automatically
%   expanded so as to match the other inputs.
%
%   S = SPBLOCK(I,J,V) does the same thing except that the size of the
%   result is determined automatically as the extreme bottom right corner of
%   any of the translated matrices in V.
%
%   Example:
%     A = sprand(3,3,0.5);
%     S1 = spblock([1 4 7],[1 4 7],{A,-2*A,A});
%     S2 = kron([1 0 0; 0 -2 0; 0 0 1],A);      % does the same thing
%
%     % a different example
%     B = sprand(3,2,0.8);
%     S3 = spblock([1 1 4],[1 4 4],{A,B,A});
%     S4 = [A B sparse(3,1); sparse(3,3) A];    % ditto
%
%   See also FSPARSE, FREPMAT, KRON, NDOP, CAT.

% S. Engblom 2012-01-24 (Revision)
% S. Engblom 2007-01-29 (Minor revision)
% S. Engblom 2004-05-24

if nargin ~= 3 & nargin ~= 4
  error('spblock:e2','Unknown syntax.');
end

% normalize input
i = i(:)-1; isz = size(i,1);
j = j(:)-1; jsz = size(j,1);
v = v(:);   vsz = size(v,1);
len = max([isz jsz vsz]);

if isz ~= len
  if isz ~= 1
    error('spblock:e1', ...
          'Index vectors must match the number of sparse blocks.');
  end
  i = frepmat(i,len);
end
if jsz ~= len
  if jsz ~= 1
    error('spblock:e1', ...
          'Index vectors must match the number of sparse blocks.');
  end
  j = frepmat(j,len);
end
if vsz ~= len
  if vsz ~= 1
    error('spblock:e3', ...
          'The sparse blocks must match the index vectors.');
  end
  v = frepmat(v,len);
end

% determine (an upper bound to) the number of nonzeros
nz = [0; cumsum(cellfun(@nnz,v))];
sz = [cellfun('size',v,1) cellfun('size',v,2)];
tnz = nz(end);

% allocate indices and values
ii = zeros(tnz,1);
jj = zeros(tnz,1);
ss = zeros(tnz,1);

% find on all blocks and build indices and values
for k = 1:len
  inds = nz(k)+1:nz(k+1);
  [ii(inds),jj(inds),ss(inds)] = find(v{k});
  ii(inds) = ii(inds)+i(k);
  jj(inds) = jj(inds)+j(k);
end

% the size of the resulting sparse matrix
if nargin < 4
  sz = sz+[i j];
  mn = max(sz);
end

% assemble the result
S = fsparse(ii,jj,ss,mn);
