function varargout = ndop(M,ix,sz,varargin)
%NDOP General N-dimensional operator.
%   S = NDOP(M,IX,SZ) constructs a sparse matrix S representing the
%   N-dimensional computational molecule M in a way such that the result of
%   M acting on a vector X is given by the sparse matrix-vector product S*X.
%
%   M is an N-dimensional full matrix with M(IX(1),...,IX(N)) being regarded
%   as the target point of the molecule.
%
%   SZ is a vector of length N specifying the size of the N-dimensional
%   rectangular computational domain.
%
%   S = NDOP(M,IX,SZ,L,R) is a syntax for variable coefficients and
%   additionally performs left and right scaling by vectors L and R. The
%   result is the same as if S had been multiplied from the left and right
%   by diagonal matrices with L and R as diagonals. Empty arguments should
%   be used when one-sided scaling is requested.
%
%   S = NDOP(M,IX,SZ,IP) is a syntax for boundary conditions and uses the
%   cell-vector IP of length N to specify the points where the molecule
%   should be applied. If IP{j} is empty, all points in that dimension are
%   mapped -- similar to MATLABs colon-notation.
%
%   S = NDOP(M,IX,SZ,[]) is a useful shorthand for mapping all points using
%   periodicity at all boundaries, see the remark on boundaries below.
%
%   S = NDOP(M,IX,SZ,IP,L,R) additionally performs scaling as before. Again,
%   empty arguments should be used for one-sided scaling.
%
%   [II,JJ,SS] = NDOP(M,IX,SZ,IP,L,R,I,J,S) is a special fast syntax for
%   boundary conditions. It does not assemble the resulting sparse matrix
%   but rather keeps it in the triplet-format for later assembly via a call
%   to FSPARSE. The result is concatenated with the obligatory inputs I, J
%   and S.
%
%   Boundaries: in the first case, boundary points in X are mapped to zero;
%   boundary points in X are all those points such that applying M would
%   cause an 'index out of bounds'-error. In the 'boundary syntax' however,
%   periodicity is used in order to be able to handle periodic boundary
%   conditions.
%
%   On return, S is a sparse square matrix of size PROD(SZ). The enumeration
%   of the solution vector X follows that of M, that is, the weights for
%   X(j-1) and X(j+1) when computing X(j) is given by
%   M(IX(1)-1,IX(2),...,IX(N)) and M(IX(1)+1,IX(2),...,IX(N))
%   respectively. Hence, "matrix-coordinates" rather than Cartesian
%   coordinates are used.
%
%   Example:
%     % classical 2nd order accurate Laplace-operator
%     % on the (N-by-N) unit square
%     N = 50; h = 1/(N-1);
%     M = [0 1 0; 1 -4 1; 0 1 0]/h^2;
%     S = ndop(M,[2 2],[N N]);
%
%     % Dirichlet BCs at all boundaries
%     BC = ndop(1,[1 1],[N N],{[1 N] []})+...   % lower/upper
%          ndop(1,[1 1],[N N],{[2:N-1] [1 N]}); % left/right
%
%     % solve -Laplace(u) = 1 with u = 0 on the boundary
%     f = ones(N); f([1 N],:) = 0; f(:,[1 N]) = 0;
%     u = (-S+BC)\f(:);
%     figure, surf(reshape(u,N,N));
%
%   See also FSPARSE, SPBLOCK.

% S. Engblom 2007-01-29 (Major revision)
% S. Engblom 2004-10-29 (Revision, changed to frepmat)
% S. Engblom 2003-03-26

% determine dimensionality
ix = ix(:); sz = sz(:);
ndim = size(ix,1);
if size(ix,1) ~= size(sz,1)
  error('ndop:e1','Ambiguous dimension of input.');
elseif tndims(M) > ndim
  error('ndop:e2','Molecule has incorrect dimension.');
end

% scaling arguments
if nargin == 5
  L = varargin{1}(:);
  R = varargin{2}(:);
elseif nargin == 6 || nargin == 9
  L = varargin{2}(:);
  R = varargin{3}(:);
else
  L = [];
  R = [];
end

% find on M
[imm,foo,ss] = find(M(:));
imm = imm(:); % fix for bug in find when M is zero scalar
ss = ss(:).';

% compute indices
imr = zeros(size(imm,1),ndim);
str = [1 cumprod(tsize(M,1:ndim))];
imm = imm-1;
for j = ndim:-1:1,
  imr(:,j) = floor(imm/str(j))+1-ix(j); % use ix as origin
  imm = rem(imm,str(j));
end

% stride for the (imagined) N-dimensional solution block
str = [1; cumprod(sz)];

% check of scaling
if ~isempty(L) && numel(L) ~= str(end)
  error('ndop:e3','Wrong size of scaling argument.');
end
if ~isempty(R) && numel(R) ~= str(end)
  error('ndop:e3','Wrong size of scaling argument.');
end

% normal syntax: apply molecule only where possible
if nargout <= 1 && (nargin == 3 || nargin == 5)
  % compute relative indices to DOFs coupled to each DOF (there will typically
  % be a zero somewhere in relind representing the target point of the
  % molecule)
  relind = imr*str(1:end-1);

  % build rowindices by finding all interior points: the location of the
  % target point together with the index of the extreme points of the
  % molecule determine which DOFs are inner DOFs
  first = max([-min(imr,[],1); zeros(1,ndim)],[],1)';
  last = sz-1+min([-max(imr,[],1); zeros(1,ndim)],[],1)';
  first = first.*str(1:end-1);
  last = last.*str(1:end-1);

  ii = 1;
  for j = 1:ndim
    % the rank of ii increases by 1
    ii = tsum(ii(:),[first(j):str(j):last(j)]',1,2);
  end

  % final assembly
  ii = ii(:);
  jj = tsum(ii,relind,1,2);
% 'boundary' or 'periodic' syntax
elseif nargout <= 1 && (nargin == 4 || nargin == 6) || ...
      (nargout == 3 && nargin == 9)
  ip = varargin{1}(:);
  if isempty(ip)
    % support of all-periodic case
    ip = cell(ndim,1);
    ipsz = sz;
  elseif size(ip,1) ~= ndim
    error('ndop:e1','Ambiguous dimension of input.');
  else
    ipsz = cellfun('prodofsize',ip)';
    ipsz(ipsz == 0) = sz(ipsz == 0); % 'colon'-notation
  end

  % build jj by indexing according to input ip, assembling all indices per row
  % at the same time (this is the last dimension in jj)
  jj = ones(size(imm,1),1);
  ii = 1;
  for j = 1:ndim
    if isempty(ip{j}), ip{j} = 1:sz(j); end
    jj = tsum(jj,mod(tsum(ip{j}(:)-1,imr(:,j),1,2),sz(j))*str(j),...
              [1:j-1 j+1],[j j+1]);
    ii = tsum(ii(:),(ip{j}(:)-1)*str(j),1,2);
  end

  % final assembly
  ii = ii(:);
  jj = reshape(jj,[],size(imr,1));
else
  error('ndop:e4','Unknown syntax.');
end

% scaling (if any)
if ~isempty(L)
  ss = L(ii)*ss;
  if ~isempty(R)
    ss = ss.*reshape(R(jj),size(ss));
  end
elseif ~isempty(R)
  ss = tprod(ss,reshape(R(jj),size(jj)),[3 2],[1 2]);
end

% output
if nargout <= 1
  varargout{1} = fsparse(ii,jj,ss,[str(end) str(end)]);
else
  % cannot use the 'row' format of fsparse for this output
  if size(ii,2) ~= size(jj,2)
    ii = frepmat(ii,[1 size(jj,2)]);
  end
  if size(ss,1) ~= size(ii,1)
    ss = frepmat(ss,size(ii,1));
  end
  varargout(1:3) = {[varargin{4}(:); ii(:)] ...
                    [varargin{5}(:); jj(:)] ...
                    [varargin{6}(:); ss(:)]};
end
