function h = stability(rho,sigma,varargin)
%STABILITY Plot of stability region for linear multistep method.
%   STABILITY(RHO,SIGMA) plots the unit root boundary locus for the
%   linear multistep method defined by the polynomials RHO and
%   SIGMA. For the ODE Y'= F(t,Y) the method is given by RHO(E)*Y(n) =
%   h*SIGMA(E)*F(t(n),Y(n)), where E is the forward shift operator.
%
%   STABILITY(RHO,SIGMA,N) uses N points on the unit circle, with N =
%   200 the default.
%
%   H = STABILITY(RHO,SIGMA,...) returns a handle H to the plot and
%   sends any additional arguments to PLOT.
%
%   Examples:
%     % explicit/implicit Euler
%     figure, hold on,
%     stability([1 -1],[0 1],'r');
%     stability([1 -1],[1 0],'r--');
%
%     % BDF2
%     stability([3/2 -2 1/2],[1 0 0],'b--');
%     axis equal, hold off,
%     legend('Euler Fwd','Euler Bwd','BDF2');
%
%   See also CONSISTENCY, PLOT, POLYVAL.

% S. Engblom 2010-01-30 (Minor modifications)
% G. Söderlind 2010-01-29

% straightforward
if nargin < 3 || ~isnumeric(varargin{1})
  N = 200;
else
  N = varargin{1};
  varargin = varargin(2:end);
end
z = exp(1i*linspace(0,2*pi,N+1));
y = polyval(rho,z)./polyval(sigma,z);
if nargout > 0
  h = plot(real(y),imag(y),varargin{:});
else
  plot(real(y),imag(y),varargin{:});
end
