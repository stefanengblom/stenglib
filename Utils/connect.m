function S = connect(z,tol,exact)
%CONNECT Connectivity information for points.
%   S = CONNECT(Z,tol) determines closest neighbours for all the (complex) 
%   coordinates Z so that find(S(:,i)) contains a complete list of indices
%   to elements Z(j) that satisfy abs(Z(j)-Z(i)) <= tol.
%
%   S = CONNECT(Z,tol,0) returns a slightly less sharp result. The list is
%   still complete but generally contains more elements. All elements Z(j)
%   are guaranteed to satisfy abs(Z(j)-Z(i)) <= 2*sqrt(2)*tol ~ 2.82*tol.
%
%   The output S is an N-by-N sparse matrix, where N is the number of
%   coordinates.
%
%   Cautionary: if tol is small compared to the largest bounding box
%   containing all points in Z, a lot of memory is allocated.
%
%   Example:
%     % 10000 random points
%     np = 10000;
%     z = 2*complex(rand(np,1),rand(np,1))-(1+1i);
%
%     S = connect(z,0.5e-1);
%     S0 = connect(z,0.5e-1,0); % faster but less sharp
%
%     figure, plot(z,'b.'); hold on
%     i = find(S0(:,1)); plot(z(i),'r.');
%     i = find(S(:,1)); plot(z(i),'c.'); plot(z(1),'k*');
%     j = find(S0(:,2)); plot(z(j),'r.');
%     j = find(S(:,2)); plot(z(j),'c.'); plot(z(2),'k*');
%     axis equal, axis([-1 1 -1 1]);

% S. Engblom 2009-09-24 (Minor revision)
% S. Engblom 2008-04-29

% check input
if ~isscalar(tol) || tol <= 0
  error('Invalid input.');
end

% column vector
z = z(:); N = size(z,1);

% extract coordinates
x = real(z); y = imag(z);

% bounding box
xmin = min(x); hx = max(x)-xmin;
ymin = min(y); hy = max(y)-ymin;

% size of each box
heff = 2*tol;

% number of boxes
nx = max(1,ceil(hx/heff)); ny = max(1,ceil(hy/heff));
NB = nx*ny;

% a possibly slightly larger bounding box
hx = nx*heff; hy = ny*heff;

% change to local coordinates in [0,1]
x = (x-xmin)/hx; y = (y-ymin)/hy;

% index into primary mesh
i1 = 1+min(floor(x*nx),nx-1);
j1 = min(floor(y*ny),ny-1);

% index into shifted mesh
i2 = 1+min(floor(x*nx+0.5),nx);
j2 = min(floor(y*ny+0.5),ny);

% assemble the result
ix = repmat(1:N,1,4);
jx = [i1 i1+NB i2+2*NB+nx i2+3*NB+nx+ny]+ ...
     [j1*nx j2*nx j1*(nx+1) j2*(nx+1)];
S = sparse(ix,jx,1);
S = S*S';

% optionally improve the result
if nargin < 3 || exact
  [ix,jx] = find(S);
  mask = find(abs(z(ix)-z(jx)) <= tol);
  S = sparse(ix(mask),jx(mask),1,N,N);
end

% this code resembles the multipole mesh more closely but is
% slightly less effective:

% % size of each box
% heff = tol/2;
% 
% % number of boxes
% nx = max(1,ceil(hx/heff)); ny = max(1,ceil(hy/heff));
% NB = nx*ny;
% 
% % a possibly slightly larger bounding box
% hx = nx*heff; hy = ny*heff;
% 
% % change to local coordinates in [0,1]
% x = (x-xmin)/hx; y = (y-ymin)/hy;
% 
% % index into an imagined nx-by-ny matrix
% ii = 1+min(floor(x*nx),nx-1);
% jj = min(floor(y*ny),ny-1);
% ii = [max(ii-1,1); ii; min(ii+1,nx)];
% jj = [max(jj-1,0) jj min(jj+1,ny-1)];
% ii = repmat(ii,1,3);
% jj = repmat(jj,3,1);
% 
% % assemble the result
% ix = repmat(1:N,1,9);
% jx = ii(:)+jj(:)*nx;
% S = sparse(ix,jx,1,N,NB); % N-by-NB (coordinate-to-box)
% S = S*S'; % N-by-N (coordinate-to-coordinate in the same box)
