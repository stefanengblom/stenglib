function [p,c] = consistency(rho,sigma,N)
%CONSISTENCY Local truncation error of linear multistep method.
%   [P,C] = CONSISTENCY(RHO,SIGMA) determines the local truncation
%   error for the linear multistep method defined by the polynomials
%   RHO and SIGMA. For the ODE Y'= F(t,Y) the method is given by
%   RHO(E)*Y(n) = h*SIGMA(E)*F(t(n),Y(n)), where E is the forward
%   shift operator.
%
%   The local error is returned as an integer P (the order) and a
%   length-2 vector C (error constants). The local error is then to
%   leading order given by C(1)*h^(P+1)*D^(P+1)*Y(n), where D =
%   d/dt. The error constant is given by C(2) = C(1)/sum(sigma).
%
%   P = CONSISTENCY(RHO,SIGMA) returns a residual polynomial P such
%   that the local error is given by POLYVAL(P,DELTA)*Y(n), where the
%   discrete difference DELTA = E-1 = exp(hD)-1 = hD+(hD)^2/2+...
%
%   P = CONSISTENCY(RHO,SIGMA,N) returns the first N terms, where N =
%   max(numel(RHO),numel(SIGMA))+2 is the default.
%
%   Examples:
%     % explicit Euler
%     [p,c] = consistency([1 -1],[0 1])
%
%     % BDF2
%     p1 = consistency([3/2 -2 1/2],[1 0 0])
%
%     % "explicit" BDF2
%     p2 = consistency([1 -4/3 1/3],2/3*[0 2 -1])
%
%     % 2nd order Adams-Bashforth
%     p3 = consistency([1 -1 0],[0 3/2 -1/2])
%
%     % (note that the actual error constants differ more than
%     % perhaps expected judging from p2 and p3:)
%     [p,c2] = consistency([1 -4/3 1/3],2/3*[0 2 -1])
%     [p,c3] = consistency([1 -1 0],[0 3/2 -1/2])
%
%     % trapezoidal rule
%     p4 = consistency([1 -1],1/2*[1 1],10)'
%
%   See also STABILITY, POLYVAL.

% S. Engblom 2010-01-31

% input
rlen = numel(rho);
slen = numel(sigma);
len = max(rlen,slen);
if nargin < 3, N = len+2; end

% code becomes much slicker if the constant term is first:
rho = [reshape(rho(end:-1:1),1,[]) zeros(1,len-rlen) ];
sigma = [reshape(sigma(end:-1:1),1,[]) zeros(1,len-slen)];

% series for log(1+x)
logpoly = [0 1./(1:N-1)]; % enough for terms 1..N
logpoly(1:2:end) = -logpoly(1:2:end);
% (understood with x = DELTA = exp(hD)-1, logpoly = hD)

% re-expand rho and sigma around x = 1
sigma1 = zeros(1,len);
rho1 = zeros(1,len);
fac = 1;
dfac = 1:len;
for d = 1:len
  % derivatives at x = 1
  sigma1(d) = sum(sigma)/fac;
  sigma = sigma(2:end).*dfac(1:end-d);
  rho1(d) = sum(rho)/fac;
  rho = rho(2:end).*dfac(1:end-d);
  fac = d*fac;
end

% residual polynomial
p = [rho1 zeros(1,N-1)]-conv(logpoly,sigma1);
if nargout > 1
  ip = find(abs(p) > eps(1000*(norm(rho1,1)+norm(sigma1,1))),1,'first');
  c = [p(ip) p(ip)/sigma1(1)];
  p = ip-2;
else
  p = p(1:N);

  % Stirling transformation (in case one prefers a residual polynomial
  % in hD rather than in DELTA = exph(hD)-1)
  %v = [1 zeros(1,N-2)];
  %q = [p(1) zeros(1,N-1)];
  %p = p(2:N);
  %for i = 2:N
  %  q(i) = v*p';
  %  v = (1:N-1)/i.*(v+[0 v(1:end-1)]);
  %end
  %p = q;

  % Matlab convention:
  p = p(end:-1:1);
end
