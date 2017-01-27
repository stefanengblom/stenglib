function scrub(ix)
%SPIN Tests for SCICOMP.

% J.C Araujo 2011-06-14

ftests = {@l_scrub1 @l_scrub2 @l_scrub3}; % remove spin3
n = size(ftests,2);
stests = cellstr(reshape(sprintf('Test #%2d',1:n)',[],n)')';

if nargin == 0, ix = 1:n; end
runtest('SCRUB (scicomp)',ftests(ix),stests(ix));
%--------------------------------------------------------------------------
function [fimpl,fexpl] = fun_s(t,y,told,yold,flag)

if flag == 1
    %Brusselator explicit
    A = 1;  B = 3;  a1=1;  a2=1;  a3=1;  a4=1;
    f(1) = a1*A - a2*B*yold(1) + a3*yold(1)*yold(1)*yold(2) - a4*yold(1);
    f(2) = a2*B*yold(1) - a3*yold(1)*yold(1)*yold(2);
    fexpl=f';
    fimpl = zeros(size(fexpl));
elseif flag == 2
    %Brusselator implicit
    A = 1;  B = 3;  a1=1;  a2=1;  a3=1;  a4=1;
    f(1) = a1*A - a2*B*y(1) + a3*y(1)*y(1)*y(2) - a4*y(1);
    f(2) = a2*B*y(1) - a3*y(1)*y(1)*y(2);
    fimpl=f';
    if nargout > 1, fexpl = zeros(size(fimpl)); end
elseif flag == 3
    %Brusselator Trapezoidal
    A = 1;  B = 3;  a1=1;  a2=1;  a3=1;  a4=1;
    f(1) = a1*A - a2*B*y(1) + a3*y(1)*y(1)*y(2) - a4*y(1);
    f(2) = a2*B*y(1) - a3*y(1)*y(1)*y(2);
    fimpl=0.5*f';
    f(1) = a1*A - a2*B*yold(1) + a3*y(1)*yold(1)*yold(2) - a4*yold(1);
    f(2) = a2*B*yold(1) - a3*yold(1)*yold(1)*yold(2);
    fexpl = 0.5*f';
end
%--------------------------------------------------------------------------
function [fimpl,fexpl] = funVDP_s(t,y,told,yold,flag)

o = 10;

if flag == 1
    %Van der Pol explicit
    f(1) = yold(2);
    f(2) = o*(1 - yold(1)^2)*yold(2) - yold(1);
    fexpl=f';
    fimpl = zeros(size(fexpl));
elseif flag == 2
    %Van der Pol implicit
    f(1) = y(2);
    f(2) = o*(1 - y(1)^2)*y(2) - y(1);
    fimpl=f';
    if nargout > 1, fexpl = zeros(size(fimpl)); end
elseif flag == 3
    %Van der Pol Trapezoidal
    f(1) = y(2);
    f(2) = o*(1 - y(1)^2)*y(2) - y(1);
    fimpl = 0.5*f';
    f(1) = yold(2);
    f(2) = o*(1 - yold(1)^2)*yold(2) - yold(1);
    fexpl = 0.5*f';
end
function ok = l_spin1
%L_SPIN1 Test of RTSAFE.

ok = 1;

x1 = [linspace(0,1,10) 1 linspace(1,2,10)]';
x2 = [linspace(1,2,10) 1 linspace(0,1,10)]';

x = rtsafe(x1,l_tfun(x1),x2,l_tfun(x2),@l_tfun);
ok = ok && all(x == 1);

%--------------------------------------------------------------------------
function [fimpl,fexpl] = funMassScrub(t,y,told,yold)

% Stiffness matrix A
n = size(y,1);

e = ones(n,1);
A = spdiags([e -2*e e], -1:1, n, n);
A(1,1) = -1;
A(n,n) = -1;
A = 10*A;

fimpl= (A*y);
if nargout > 1
  fexpl = zeros(size(fimpl));
end
%--------------------------------------------------------------------------
function ok = l_scrub1
%L_SCRUB1 Test of ODE1S.
% J. C. Araujo Cabarcas 2011-06-13

ok = 1;

n = 100;

TOL = 0.1; % just guessing

% Brusselator
f = @fun_s;
xf = 10;
yi = [0 0];

rtols = [1e-4 3e-5 1e-5 3e-6 1e-6];

x = linspace(0,xf,n);
sol0 = ode15s(f,[0 xf],yi,odeset('reltol',1e-6,'abstol',1e-8),[],[],2);
y0 = interp1(sol0.x',sol0.y',x)';

i = 1;
for rtol = rtols
    atol = 1e-2*rtol;
    
    sol1=ode1s(f,[0 xf],yi,{'rtol',rtol,'atol',atol},1);
    sol2=ode1s(f,[0 xf],yi,{'rtol',rtol,'atol',atol},2);
    sol3=ode1s(f,[0 xf],yi,{'rtol',rtol,'atol',atol},3);

    y1 = interp1(sol1.x',sol1.y',x)';
    y2 = interp1(sol2.x',sol2.y',x)';
    y3 = interp1(sol3.x',sol3.y',x)';

    e1(i) = sqrt(sum(sum((y1-y0).^2)*(x(2)-x(1))))/xf;
    e2(i) = sqrt(sum(sum((y2-y0).^2)*(x(2)-x(1))))/xf;
    e3(i) = sqrt(sum(sum((y3-y0).^2)*(x(2)-x(1))))/xf;
    
    ok = ok && e1(i) < TOL && e2(i) < TOL && e3(i) < TOL;
    i = i+1;
end

%e1
%e2
%e3
%--------------------------------------------------------------------------
function ok = l_scrub2
%L_SCRUB2 Test of ODE1S.
% J. C. Araujo Cabarcas 2011-06-13

ok = 1;

n = 100;

TOL = 0.5; % just guessing

% Van de Pol
f = @funVDP_s;
xf = 20;
yi = [2 0];

rtols = [1e-4 3e-5 1e-5 3e-6 1e-6];

x = linspace(0,xf,n);
sol0 = ode15s(f,[0 xf],yi,odeset('reltol',1e-6,'abstol',1e-8),[],[],2);
y0 = interp1(sol0.x',sol0.y',x)';

i = 1;
for rtol = rtols
    atol = 1e-2*rtol;
    
    sol1=ode1s(f,[0 xf],yi,{'rtol',rtol,'atol',atol},1);
    sol2=ode1s(f,[0 xf],yi,{'rtol',rtol,'atol',atol},2);
    sol3=ode1s(f,[0 xf],yi,{'rtol',rtol,'atol',atol},3);

    y1 = interp1(sol1.x',sol1.y',x)';
    y2 = interp1(sol2.x',sol2.y',x)';
    y3 = interp1(sol3.x',sol3.y',x)';

    e1(i) = sqrt(sum(sum((y1-y0).^2)*(x(2)-x(1))))/xf;
    e2(i) = sqrt(sum(sum((y2-y0).^2)*(x(2)-x(1))))/xf;
    e3(i) = sqrt(sum(sum((y3-y0).^2)*(x(2)-x(1))))/xf;
    
    ok = ok && e1(i) < TOL && e2(i) < TOL && e3(i) < TOL;
    i = i+1;
end

e1
e2
e3
%--------------------------------------------------------------------------
function ok = l_scrub3
%L_SCRUB3 Test of ODE1S. Test based on the Diffusion equation, using a
% pseudo Laplacian operator A as RHS of u_t = Au.
% J. C. Araujo Cabarcas 2011-11-4
ok = 1;

tf = 10;
n = 11;

u0 = zeros(n,1);

for i = 1:n
    u0(i) = i-1;
end
   
for flag = 1:4
    if flag == 1
        M = 5*eye(n, n);
    elseif flag == 2
        e = ones(n,1);
        M = spdiags([e 5*e e], -1:1, n, n);
        M(1,1) = 1; M(n,n) = 1;
    elseif flag == 3
        e = ones(n,1);
        M = spdiags([3*e 5*e e], -1:1, n, n);
        M(1,1) = 1; M(n,n) = 1;
    elseif flag == 4
        e = ones(n,1);
        M = spdiags([3*e 5*e e], -1:1, n, n);
        M(1,1) = 1; M(n,n) = 1;

        M = 0.1*M;
    else
        M = eye( n,n );
    end
    rtol = 1e-4;

    % computing solutions
    options = odeset('Mass',M,'RelTol', rtol ); % options for ode23s
    sol1 =  ode1s(@funMassScrub,[0 tf],u0,{'mass',M,'rtol',rtol,'atol',rtol/100});
    x = sol1.x;
    y1 = sol1.y;

    y = zeros(size(u0),size(x) );
    y(:,1) = u0;

    %figure(1)
    %plot( sol1.x,sol1.y );

    %sol2 = ode23t(@funMass,[0 tf],u0,options);

    %figure(2)
    %plot( sol2.x,sol2.y );

    e(1)  = 0; % same initial condition
    ey(1) = 0; % same initial condition

    for i = 2:size(sol1.y,2)
        sol = ode23t(@funMassScrub,[0 x(i)],u0,options);
        y(:,i) = sol.y(:,size(sol.y,2));

        % comparison calculations
        ey(i)   = y(:,i)'*y(:,i);
        dy(:,i) = y(:,i) - y1(:,i);
        e(i)    = dy(:,i)'*dy(:,i); % squared difference (each time)
    end
    ey = sqrt( ey );

    e  = sqrt( e  );
    error = sqrt(e'*e)/norm(ey); % shall I take error/100?
    
    ok = error < 0.1;
    if ~ok
        break;
    end
end
