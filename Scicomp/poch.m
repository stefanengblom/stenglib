function y = poch(a,x,b)
%POCH Pochhammer's function.
%  Y = POCH(A,X), for a scalar A and a vector X, computes the value of
%  GAMMA(A+X)./GAMMA(A). The output Y is a vector the same size as X.
%
%  If A is a vector, the ouput is Y(i) = prod_i POCH(A(i),X).
%
%  Y = POCH(A,B,X) is the same as POCH(A,X)./POCH(B,X).
%
%  See also GAMMA, GAMMALN, PSI.

% S. Engblom 2006-01-12

if nargin == 3
  t = x(:);
  x = b(:);
  a = a(:); b = t;
  y = exp(sum(l_gammaln(tsum(a,x,[1],[2])),1)- ...
          sum(l_gammaln(tsum(b,x,[1],[2])),1)- ...
          sum(l_gammaln(a))+sum(l_gammaln(b)));
else
  a = a(:);
  x = x(:);
  y = exp(sum(l_gammaln(tsum(a,x,[1],[2])),1)- ...
          sum(l_gammaln(a)));
end

%--------------------------------------------------------------------------
function y = l_gammaln(x)
%L_GAMMALN Fix for bug in GAMMALN.
%   Annoyingly, GAMMALN is incorrectly imlemented. This function is a
%   patch for the case when there might be negative elements in
%   X. Still, complex values of X are not allowed.

isneg = find(x < 0);
x(isneg) = 1-x(isneg);
y = gammaln(x);
y(isneg) = -y(isneg)+log(pi./sin(pi*(x(isneg))));

%--------------------------------------------------------------------------
