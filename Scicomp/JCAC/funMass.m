function [fimpl,fexpl] = funMass(t,y,told,yold)

% Stiffness matrix A
n = size(y,1);

e = ones(n,1);
A = spdiags([e -2*e e], -1:1, n, n);
A(1,1) = -1;
A(n,n) = -1;
A = 10*A;

fimpl= (A*(y.*y));
if nargout > 1
  fexpl = zeros(size(fimpl));
end