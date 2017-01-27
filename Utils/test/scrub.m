function scrub(ix)
%SCRUB Tests for UTILS.

% S. Engblom 2004-10-29 (revision)
% S. Engblom 2003-03-27

if nargin == 0, ix = 1:5; end
ftests = {@l_scrub1 @l_scrub2 @l_scrub3 @l_scrub4 @l_scrub5};
stests = {'Test #1 (ndop)' 'Test #2 (ndop)' 'Test #3 (ndop)' ...
          'Test #4 (ndop)' 'Test #5 (ndop)'};
runtest('SCRUB (utils)',ftests(ix),stests(ix));

%--------------------------------------------------------------------------
function ok = l_scrub1
%L_SCRUB1 Variable coefficients with periodic BCs.

ok = 1;

M = rand(3);
ix = [2 2];
sz = [40 40];
for i = 1:4
  % various 'empty' syntaxes
  switch i
   case 1, l = rand(sz); r = rand(sz);
   case 2, l = []; r = rand(sz);
   case 3, l = rand(sz); r = [];
   case 4, l = []; r = [];
  end

  S1 = ndop(M,ix,sz,l,r);
  S1 = S1+ndop(M,ix,sz,{[1 sz(2)],[]},l,r);
  S1 = S1+ndop(M,ix,sz,{[2:sz(1)-1],[1 sz(2)]},l,r);

  % the same thing
  S2 = ndop(M,ix,sz,[],l,r);
  S = ndop(M,ix,sz,[]);
  if ~isempty(l)
    S = sparse(1:prod(sz),1:prod(sz),l)*S;
  end
  if ~isempty(r)
    S = S*sparse(1:prod(sz),1:prod(sz),r);
  end
  ok = ok & nnz(S-S1) == 0 & nnz(S1-S2) == 0;
end

% 1-D example with variable coefficients according to P. Sjöberg
sz = [100];       % size of the discretization
h = rand(sz-1,1); % something like a varying stepsize...
M = [1 -1]';

% zeros where operator is undefined (boundaries)
gam = [0; 1./(h(1:end-1).*h(2:end).*(h(1:end-1)+h(2:end))); 0];
S1 = ndop(M,2,sz,[0; h.^2].*gam,[]);
S2 = ndop(M,1,sz,[h.^2; 0].*gam,[]);
S1D = (S1+S2);

% same thing, but 2-D instead
sz = [sz 10];
h = repmat(h,1,10);
M = [1 -1]';
gam = [zeros(1,10); ...
       1./(h(1:end-1,:).*h(2:end,:).*...
	   (h(1:end-1,:)+h(2:end,:))); ...
       zeros(1,10)];
hh = h.^2;
S1 = ndop(M,[2 1],sz,[zeros(1,10); hh].*gam,[]);
S2 = ndop(M,[1 1],sz,[hh; zeros(1,10)].*gam,[]);
S2D = (S1+S2);

x1 = rand(sz);
u1 = S1D*x1;
u2 = S2D*x1(:);

ok = ok & norm(u1(:)-u2) == 0;

%--------------------------------------------------------------------------
function ok = l_scrub2
%L_SCRUB2 Larger tests: comparision with true operators.

ok = 1;

sz = 10000;

x = linspace(0,1,sz(1))';
h = x(2)-x(1);
u = 1+2*x-3*x.^2+0.25*x.^3;

M = [-1 0 1]'/2/h;
ix = 2;
S1 = ndop(M,ix,sz);

M = [-3 4 -1]'/2/h;
ix = 1; ip = {1};
S1 = S1+ndop(M,ix,sz,ip);

M = [1 -4 3]'/2/h;
ix = 3; ip = {sz(end)};
S1 = S1+ndop(M,ix,sz,ip);

u1 = S1*u;
d1 = u1-(2-6*x+3/4*x.^2);

% periodic BCs
x(end) = [];
sz = sz-1;
u = sin(2*pi*x);
M = [1 -2 1]'/h^2;
ix = 2;
S2 = ndop(M,ix,sz);

ip = {[1 sz(end)]};
S2 = S2+ndop(M,ix,sz,ip);

u2 = S2*u;
d2 = u2-(-4*pi^2*sin(2*pi*x));

% 2-D example
sz = [400 500];
x = linspace(-1,1,sz(2)+1); x(end) = [];
y = linspace(0,1,sz(1))'; y = y(end:-1:1);
X = repmat(x,sz(1),1);
Y = repmat(y,1,sz(2));
hx = x(2)-x(1); hy = y(1)-y(2);
u = sin(pi*X)+Y.*cos(pi*X)+(1-Y.^2).*Y;
du = pi*cos(pi*X)-Y.*pi.*sin(pi*X)+2*(cos(pi*X)+1-3*Y.^2);

% operator is u_x+2*u_y
Mx = [0 0 0; -1/2/hx 0 1/2/hx; 0 0 0];
My = [0 2/2/hy 0; 0 0 0; 0 -2/2/hy 0];
M = Mx+My;
ix = [2 2];
S3 = ndop(M,ix,sz);

% include single boundaries
S3b = ndop(Mx,ix,sz,{[1 sz(1)] [2:sz(2)-1]});
S3b = S3b+ndop(My,ix,sz,{[2:sz(1)-1] [1 sz(2)]});

% periodic molecule in x
S3b = S3b+ndop(Mx,ix,sz,{[] [1 sz(2)]});

% modified (skew) molecule in y
Myb1 = [3 -4 1]'/hy;
ix = [1 1];
S3b = S3b+ndop(Myb1,ix,sz,{1 []});
Myb2 = [-1 4 -3]'/hy;
ix = [3 1];
S3b = S3b+ndop(Myb2,ix,sz,{sz(1) []});

S3 = S3+S3b; 
u3 = reshape(S3*u(:),size(u));

d3 = u3-du;

% (Padé-) operator is u_x+2*u_y
sz = [100 80];
x = linspace(-1,1,sz(2)+1); x(end) = [];
y = linspace(0,1,sz(1))'; y = y(end:-1:1);
X = repmat(x,sz(1),1);
Y = repmat(y,1,sz(2));
hx = x(2)-x(1); hy = y(1)-y(2);
u = sin(pi*X)+Y.*cos(pi*X)+(1-Y.^2).*Y;
du = pi*cos(pi*X)-Y.*pi.*sin(pi*X)+2*(cos(pi*X)+1-3*Y.^2);

Mxl = [1 4 1]/6;
Mxr = [-1 0 1]/2/hx;
Myl = [1 4 1]'/6;
Myr = [1 0 -1]'/hy;
% periodic molecule in x
S4xl = ndop(Mxl,[1 2],sz,[]);
S4xr = ndop(Mxr,[1 2],sz,[]);
S4yl = ndop(Myl,[2 1],sz);
S4yr = ndop(Myr,[2 1],sz);

% modified (skew) molecule in y
Myb1l = [1 2]'/12;
Myb1r = [5 -4 -1]'/12/hy;
S4byl = ndop(Myb1l,[1 1],sz,{1 []});
S4byr = ndop(Myb1r,[1 1],sz,{1 []});

Myb2l = [2 1]'/12;
Myb2r = [1 4 -5]'/12/hy;
S4byl = S4byl+ndop(Myb2l,[2 1],sz,{sz(1) []});
S4byr = S4byr+ndop(Myb2r,[3 1],sz,{sz(1) []});

S4yl = S4yl+S4byl;
S4yr = S4yr+S4byr;

u4 = reshape(S4xl\(S4xr*u(:))+S4yl\(S4yr*u(:)),size(u));

d4 = u4-du;

ok = ok & all([norm(d1) < 3e-7 ...
	       norm(d2) < 1e-4 ...
	       norm(d3(:)) < 0.04 ...
	       norm(d4(:)) < 5e-5]) & ...
     nnz(S1) == 20002 & nnz(S1) == nzmax(S1) & ...
     nnz(S2) == 29997 & nnz(S2) == nzmax(S2) & ...
     nnz(S3) == 801000 & nnz(S3) == nzmax(S3) & ...
     nnz(S3b) == 8184 & nnz(S3b) == nzmax(S3b) & ...
     nnz(S4xl) == 24000 & nnz(S4xl) == nzmax(S4xl) & ...
     nnz(S4xr) == 16000 & nnz(S4xr) == nzmax(S4xr) & ...
     nnz(S4yl) == 23840 & nnz(S4yl) == nzmax(S4yl) & ...
     nnz(S4yr) == 16160 & nnz(S4yr) == nzmax(S4yr);

%---------------------------------------------------------------------------
function ok = l_scrub3
%L_SCRUB3 Performance tests.

ok = 1;

% 3-D Laplace -- 10^6 DOFs
M = [0 0 0; 0 1 0; 0 0 0];
M(:,:,3) = M;
M(:,:,2) = [0 1 0; 1 -8 1; 0 1 0];
sz = [100 100 100];
S = ndop(M,[2 2 2],sz);
x = rand(10^6,1); y = S*x;
x = rand(10^6,1); y = S*x;
[ii,jj,ss] = ndop(1,[1 1 1],sz,{[1 100] [] []},[],[],[],[],[]);
[ii,jj,ss] = ndop(1,[1 1 1],sz,{[2:99] [1 100] []},[],[],ii,jj,ss);
[ii,jj,ss] = ndop(1,[1 1 1],sz,{[2:99] [2:99] [1 100]},[],[],ii,jj,ss);
S = S+fsparse(ii,jj,ss,[prod(sz) prod(sz)]);

%---------------------------------------------------------------------------
function ok = l_scrub4
%L_SCRUB4 Performance tests.

ok = 1;

% 2-D Biharmonic operator -- 10^6 DOFs
M = [0 0 1 0 0; 0 2 -8 2 0; 1 -8 20 -8 1; 0 2 -8 2 0; 0 0 1 0 0];
sz = [1000 1000];
S = ndop(M,[3 3],sz);
x = rand(10^6,1); y = S*x;
x = rand(10^6,1); y = S*x;
[ii,jj,ss] = ndop(1,[1 1],sz,{[1 2 999 1000] []},[],[],[],[],[]);
[ii,jj,ss] = ndop(1,[1 1],sz,{[3:998] [1 2 999 1000]},[],[],ii,jj,ss);
S = S+fsparse(ii,jj,ss,[prod(sz) prod(sz)]);

%---------------------------------------------------------------------------
function ok = l_scrub5
%L_SCRUB5 Performance tests.

ok = 1;

% 6-D Laplace -- 10**6 DOFs
M = zeros([3 3 3 3 3 3]);
M(2,2,2,2,2,2) = -2^6;
M([1 3],2,2,2,2,2) = 1;
M(2,[1 3],2,2,2,2) = 1;
M(2,2,[1 3],2,2,2) = 1;
M(2,2,2,[1 3],2,2) = 1;
M(2,2,2,2,[1 3],2) = 1;
M(2,2,2,2,2,[1 3]) = 1;
sz = [10 10 10 10 10 10];
S = ndop(M,[2 2 2 2 2 2],sz);
x = rand(10^6,1); y = S*x;
x = rand(10^6,1); y = S*x;

[ii,jj,ss] = ndop(1,[1 1 1 1 1 1],sz,...
		  {[1 10] [] [] [] [] []},[],[],[],[],[]);
[ii,jj,ss] = ndop(1,[1 1 1 1 1 1],sz,...
		  {[2:9] [1 10] [] [] [] []},[],[],ii,jj,ss);
[ii,jj,ss] = ndop(1,[1 1 1 1 1 1],sz,...
		  {[2:9] [2:9] [1 10] [] [] []},[],[],ii,jj,ss);
[ii,jj,ss] = ndop(1,[1 1 1 1 1 1],sz,...
		  {[2:9] [2:9] [2:9] [1 10] [] []},[],[],ii,jj,ss);
[ii,jj,ss] = ndop(1,[1 1 1 1 1 1],sz,...
		  {[2:9] [2:9] [2:9] [2:9] [1 10] []},[],[],ii,jj,ss);
[ii,jj,ss] = ndop(1,[1 1 1 1 1 1],sz,...
		  {[2:9] [2:9] [2:9] [2:9] [2:9] [1 10]},[],[],ii,jj,ss);
S = S+fsparse(ii,jj,ss,[prod(sz) prod(sz)]);
%---------------------------------------------------------------------------
