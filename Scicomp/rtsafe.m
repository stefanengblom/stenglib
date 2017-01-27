function x = rtsafe(x1,f1,x2,f2,fun,fargs,tol,maxits)
%RTSAFE Scalar nonlinear solver.
%   X = RTSAFE(X1,F1,X2,F2,FUN,FARGS,TOL,MAXITS) computes a solution X
%   to the equation FUN(X) = 0. FUN has the signature [F,DF] =
%   FUN(X,FARGS{:}) and should compute both function values and
%   derivatives. On entry, X1 and X2 must bracket a root and the
%   function values F1 and F2 must therefore satisfy F1.*F2 <= 0. Only
%   the signs of F1 and F2 need to be correct in order for the
%   algorithm to work properly.
%
%   If there are more than one root, RTSAFE will find one of them
%   (including possible singularities). If FUN is identical to zero on
%   the whole interval [X1 X2], RTSAFE will converge to the midpoint.
%
%   FARGS is a cell-array and should contain additional inputs to
%   FUN. An empty cell is the default.
%  
%   TOL = [RTOL ATOL] is the absolute and relative tolerance. If TOL
%   is a scalar, then RTOL = TOL and ATOL = 1e-2*TOL are used. The
%   default is TOL = [10*eps 0.1*eps].
%
%   MAXITS is the maximum number of iterations allowed. The default is
%   100.
%
%   The method used is a fail-safe combination of Newtons method and
%   bisection. Only limited error-checking is performed.

% S. Engblom 2004-03-03

if any(f1.*f2 > 0), error('Roots must be bracketed.'); end

% defaults
if nargin < 8
  maxits = 100;
  if nargin < 7
    tol = 10*eps;
    if nargin < 6
      fargs = {};
    end
  end
end

% absolute/relative tolerance
rtol = tol(1);
if prod(size(tol)) == 1
  atol = 1e-2*rtol;
else
  atol = tol(2);
end

% form initial brackets according to the following rules:
%   xl(f2 < 0 f1 < 0 f2 == 0 f1 == 0) = [x2 x1 x2 x1]
%   xh(f1 > 0 f2 > 0 f1 == 0 f2 == 0) = [x1 x2 x1 x2]
% (the brackets has to be formed exactly as written)
xl = zeros(size(x1));
xh = zeros(size(x2));
ix1 = find(f1 < 0);
ix2 = find(f2 < 0);
xl(ix2) = x2(ix2); xl(ix1) = x1(ix1);
ix1 = find(f1 > 0);
ix2 = find(f2 > 0);
xh(ix1) = x1(ix1); xh(ix2) = x2(ix2);

ix1 = find(f1 == 0);
ix2 = find(f2 == 0);
xl(ix2) = x2(ix2); xl(ix1) = x1(ix1);
xh(ix1) = x1(ix1); xh(ix2) = x2(ix2);

% evaluation in midpoint
x = 0.5*(xl+xh);
[f,df] = feval(fun,x,fargs{:});

% step and 'step before last'
dxold = xh-xl;
dx = dxold;

for i = 1:maxits
  % prefer bisection whenever the Newton step jumps out of brackets or
  % when the size of the interval is not decreasing fast enough
  ix1 = (((x-xh).*df-f).*((x-xl).*df-f) >= 0 | ...
         abs(2*f) > abs(dxold.*df));
  ix2 = find(~ix1);
  ix1 = find(ix1);

  dxold = dx;

  % bisection step
  dx(ix1) = 0.5*(xh(ix1)-xl(ix1));
  x(ix1) = xl(ix1)+dx(ix1);

  % Newton step
  dx(ix2) = f(ix2)./df(ix2);
  x(ix2) = x(ix2)-dx(ix2);

  if all(abs(dx) <= rtol*abs(x)+atol), return; end

  % new evaluation and brackets
  [f,df] = feval(fun,x,fargs{:});
  ix1 = find(f <= 0);
  ix2 = find(f >= 0);
  xl(ix1) = x(ix1);
  xh(ix2) = x(ix2);
end
error('Maximum number of iterations reached.');
