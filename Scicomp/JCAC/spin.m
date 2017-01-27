function spin(ix)
%SPIN Tests for SCICOMP.

% J. C. Araujo Cabarcas 2011-06-13 (spin2 and 3)
% S. Engblom 2008-06-17

ftests = {@l_spin1 @l_spin2 @l_spin3}; % remove spin3
n = size(ftests,2);
stests = cellstr(reshape(sprintf('Test #%2d',1:n)',[],n)')';

if nargin == 0, ix = 1:n; end
runtest('SPIN (scicomp)',ftests(ix),stests(ix));

%--------------------------------------------------------------------------
function [f,df] = l_tfun(x)
%L_TFUN Small testfunction for L_SPIN1.

f = (0.5 <= x & x <= 1.5).*(1-x)+(x < 0.5)*0.5-(1.5 < x)*0.5;
if nargout > 1, df = (0.5 <= x & x <= 1.5).*(-1); end

%--------------------------------------------------------------------------
function ok = l_spin1
%L_SPIN1 Test of RTSAFE.

ok = 1;

x1 = [linspace(0,1,10) 1 linspace(1,2,10)]';
x2 = [linspace(1,2,10) 1 linspace(0,1,10)]';

x = rtsafe(x1,l_tfun(x1),x2,l_tfun(x2),@l_tfun);
ok = ok && all(x == 1);

%--------------------------------------------------------------------------
function [fimpl,fexpl] = l_funbruss(t,y,told,yold,flag)
%L_FUNBRUSS Testfunction (the Brusselator) for ODE1s.

A = 1; B = 3; a1 = 1; a2 = 1; a3 =1; a4 = 1;
if flag == 1
  % explicit
  f(1) = a1*A - a2*B*yold(1) + a3*yold(1)*yold(1)*yold(2) - a4*yold(1);
  f(2) = a2*B*yold(1) - a3*yold(1)*yold(1)*yold(2);
  fexpl=f';
  fimpl = zeros(size(fexpl));
elseif flag == 2
  % implicit
  f(1) = a1*A - a2*B*y(1) + a3*y(1)*y(1)*y(2) - a4*y(1);
  f(2) = a2*B*y(1) - a3*y(1)*y(1)*y(2);
  fimpl=f';
  if nargout > 1, fexpl = zeros(size(fimpl)); end
elseif flag == 3
  % trapezoidal
  f(1) = a1*A - a2*B*y(1) + a3*y(1)*y(1)*y(2) - a4*y(1);
  f(2) = a2*B*y(1) - a3*y(1)*y(1)*y(2);
  fimpl=0.5*f';
  f(1) = a1*A - a2*B*yold(1) + a3*y(1)*yold(1)*yold(2) - a4*yold(1);
  f(2) = a2*B*yold(1) - a3*yold(1)*yold(1)*yold(2);
  fexpl = 0.5*f';
end

%--------------------------------------------------------------------------
function [fimpl,fexpl] = l_funVDP(t,y,told,yold,flag)
%L_FUNVDP Testfunction (Van der Pol oscillator) for ODE1S.

o = 10;
if flag == 1
  % explicit
  f(1) = yold(2);
  f(2) = o*(1 - yold(1)^2)*yold(2) - yold(1);
  fexpl=f';
  fimpl = zeros(size(fexpl));
elseif flag == 2
  % implicit
  f(1) = y(2);
  f(2) = o*(1 - y(1)^2)*y(2) - y(1);
  fimpl=f';
  if nargout > 1, fexpl = zeros(size(fimpl)); end
elseif flag == 3
  % trapezoidal
  f(1) = y(2);
  f(2) = o*(1 - y(1)^2)*y(2) - y(1);
  fimpl = 0.5*f';
  f(1) = yold(2);
  f(2) = o*(1 - yold(1)^2)*yold(2) - yold(1);
  fexpl = 0.5*f';
end

%--------------------------------------------------------------------------
function ok = l_spin2
%L_SPIN2 Test of ODE1S.

ok = 1;

% the Brusselator
f = @l_funbruss;
xf = 10;
y0 = [0 0];

% compare to reference solution from ode15s
n = 100;
x = linspace(0,xf,n);
sol0 = ode15s(f,[0 xf],y0,odeset('reltol',1e-6,'abstol',1e-8),[],[],2);
sol1 = ode1s(f,[0 xf],y0,{'rtol',1e-4,'atol',1e-6},1);
sol2 = ode1s(f,[0 xf],y0,{'rtol',1e-4,'atol',1e-6},2);
sol3 = ode1s(f,[0 xf],y0,{'rtol',1e-4,'atol',1e-6},3);

y0 = deval(sol0,x);
y1 = interp1(sol1.x',sol1.y',x)';
y2 = interp1(sol2.x',sol2.y',x)';
y3 = interp1(sol3.x',sol3.y',x)';

% integrated L^2-error
e1 = sqrt(sum(sum((y1-y0).^2)*(x(2)-x(1))))/xf;
e2 = sqrt(sum(sum((y2-y0).^2)*(x(2)-x(1))))/xf;
e3 = sqrt(sum(sum((y3-y0).^2)*(x(2)-x(1))))/xf;

ok = ok && all([e1 e2 e3] < [0.02 0.01 0.06]);

%--------------------------------------------------------------------------
function ok = l_spin3
%L_SPIN2 Test of ODE1S.

ok = 1;

% Van der Pol
f = @l_funVDP;
xf = 20;
y0 = [2 0];

n = 100;
x = linspace(0,xf,n);
sol0 = ode15s(f,[0 xf],y0,odeset('reltol',1e-6,'abstol',1e-8),[],[],2);
sol1 = ode1s(f,[0 xf],y0,{'rtol',1e-4,'atol',1e-6},1);
sol2 = ode1s(f,[0 xf],y0,{'rtol',1e-4,'atol',1e-6},2);
sol3 = ode1s(f,[0 xf],y0,{'rtol',1e-4,'atol',1e-6},3);

y0 = deval(sol0,x);
y1 = interp1(sol1.x',sol1.y',x)';
y2 = interp1(sol2.x',sol2.y',x)';
y3 = interp1(sol3.x',sol3.y',x)';

e1 = sqrt(sum(sum((y1-y0).^2)*(x(2)-x(1))))/xf;
e2 = sqrt(sum(sum((y2-y0).^2)*(x(2)-x(1))))/xf;
e3 = sqrt(sum(sum((y3-y0).^2)*(x(2)-x(1))))/xf;

ok = ok && all([e1 e2 e3] < [0.5 0.5 0.05]);

%--------------------------------------------------------------------------

