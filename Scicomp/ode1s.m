function [t,y] = ode1s(fun,tspan,y0,opts,varargin)
%ODE1S Solve stiff and strongly nonlinear ODE. Split-step Euler method.
%   [T,Y] = ODE1S(FUN,[TSTART TEND],Y0,OPTS,...) solves the ODE dY/dt
%   = FUN(t,Y,...) from TSTART to TEND with initial data Y0. OPTS is a
%   structure with options, see the table below.
%
%   For the right-hand side, the special syntax [dYimpl,dYexpl] =
%   FUN(tnew,Ynew,told,Yold,...) is used. Here Ynew(Yold) is the
%   latest iterate(solution) at the next(previous) step. This syntax
%   makes it easy to implement split-step semi-implicit methods. With
%   only one output, dYimpl = FUN(tnew,Ynew,told,Yold,...) is
%   understood.
%
%   The method is thus a forward/backward Euler combination. The error
%   is estimated by taking half-steps and a step-size control strategy
%   based on digital filters is used.
%
%   SOL = ODE1S(...) is also an allowed syntax, returning a solution
%   structure instead.
%
%   Field        Value/{Default}        Description
%   -----------------------------------------------------------------------
%   atol         Vector/scalar {1e-6}   Absolute tolerance.
%   rtol         Vector/scalar {1e-4}   Relative tolerance.
%
%   adaptivity   {'eps'} | 'epus' |     Error per step, per unit step, or
%                'eq'                   equilibrium control.
%   errnorm1     Function {@abs}        Inner and outer error-norms, see 
%   errnorm2     Function {@norm}       below.
%
%   initialstep  Scalar {0}             The initial step, zero means that
%                                       this is determined automatically.
%   maxstep      Scalar {inf}           Maximum step, 
%   minstep      Scalar {0}             minimum step.
%
%   nlinauto     {1} | 2 | 3            Tuning of the nonlinear solver. A
%                                       higher value means a more
%                                       suspicious solver, prepared to take 
%                                       more drastic measures when
%                                       poor convergence is detected.
%
%   reportfun    Function {[]}          Report function, see below.
%   info         {0} | 1 | 2            Level of information. 0:
%                                       only serious errors, 1:
%                                       include warnings, 2:
%                                       include final statistics
%                                       (not yet implemented).
%
%   In EPS-adaptivity, the control objective is to keep
%     errnorm2(TOL./ERR) = 1,
%   where in terms of the tentative solution y,
%     TOL = max(rtol.*errnorm1(y),atol),
%     ERR = errnorm1(error estimate).
%   The default choice is errnorm1(z) = abs(z) and errnorm2(z) =
%   norm(z), and is the usual component-wise control. An alternative
%   is the less stringent "norm control" defined by errnorm1(z) =
%   norm(z) and errnorm2(z) = abs(z). Other variants are easily
%   constructed. Note that the norms are automatically normalized to
%   unit measure so that errnorm2(errnorm1(ones(size(solution
%   vector)))) = 1.
%
%   In EPUS-adaptivity, one rather uses
%     TOL = h*max(rtol.*errnorm1(FUN(y)),atol),
%   while in EQ-adaptivity,
%     TOL = max(rtol.*errnorm1(y-y0),atol),
%   ensuring small steps near an equilibrium.
%
%   The report function is called according to the three
%   self-explaining syntaxes STATUS = reportfun([TSTART
%   TEND],suggested title,'init'), STATUS = reportfun(t,y,[],...),
%   and, STATUS = reportfun(t,y,'done',...). An occasional useful
%   feature is that the integration is automatically aborted whenever
%   the returned value STATUS is non-zero.
%
%   See also ODE15S, ODE23S.

% S. Engblom 2012-09-18 (info added)
% S. Engblom 2010-08-25 (Port from earlier code)
% S. Engblom 2007-03-13

yy = y0(:);
ndof = size(yy,1);

% default options
optdef.atol = 1e-6;
optdef.rtol = 1e-4;
optdef.adaptivity = 'eps';
optdef.errnorm1 = @abs;
optdef.errnorm2 = @norm;
optdef.initialstep = 0;
optdef.maxstep = inf;
optdef.minstep = 0;
optdef.nlinauto = 1;
optdef.reportfun = [];
optdef.info = 0;

% merge defaults with actual inputs
if nargin > 3
  opts = parseopts(optdef,opts);
else
  opts = optdef;
end
assignopts(opts);
info = opts.info; % (silly, but there is a Matlab-function 'info'...)

% EPS/EPUS/EQ control
[foo,adaptivity] = fsetop('ismember',{adaptivity},{'eps' 'epus' 'eq'});
if adaptivity == 0
  error('Unknown type of adaptivity.');
end

% normalize the norms...
NM = 1/feval(errnorm2,(feval(errnorm1,(ones(ndof,1)))));

% setup reporter
if ischar(reportfun)
  reportfun = str2func(reportfun);
end
if ~isempty(reportfun)
  status = reportfun([tspan(1) tspan(end)],'Solution progress...','init');
else
  reportfun = @report;
  status = reportfun([],[],'none');
end
% *** stat counter of some kind?

% step-size controller
Const = [0 0 1];      % constant step-size
Ho110 = [1 0 0];      % elementary controller (Deadbeat I)
H110 = [1/3 0 0];     % I controller (convolution filter)
PI34 = [7/10 -2/5 0]; % PI controller (Gustafsson)
PI42 = [3/5 -1/5 0];  % PI controller for nonstiff methods
H211PI = [1/6 1/6 0]; % 1st order LP filter of PI structure
gain = 1/4;           % gain in [1/6 1/2] (gain = 1/2 yields Ho211)
H211b = gain*[1 1 1]; % general purpose 1st order LP filter
% controller used:
ctrl = H211b;
% *** ctrl and Mfac (below) could be options

% (inverse of the) order of the scheme
switch adaptivity
 case 1 % EPS
  order = 1/2; 
 case {2 3} % EPUS/EQ
  order = 1;
end

% rejection logic is based on step-ratios
Mfac = 2; % rejection whenever err > Mfac*TOL

% reject when FILTER(err) > Mfac*TOL, where FILTER is the product of
% the controller and the limiters
rhomin = l_limiter(l_lpfilter(ctrl,l_limiter(1/Mfac^order,2),1,1),1);

% start-up
cerr = 1;
rho = 1;
% *** direct manipulation of the controller's state is done in three
% places in the code below when the nonlinear solver fails; it is not
% completely clear how this should be done in the best way

% build options for the nonlinear solver
if nlinauto == 1
  nlinopts.kappa = 10^-1.5;       % required extra precision
  nlinopts.thetaREF = 0.4;        % target convergence factor
  nlinopts.thetaJAC = 0.4;        % recompute Jacobian
  nlinopts.thetaLU = 0.4;         % refactorize Jacobian
  nlinopts.etainit = 0.9/(1-0.9); % (total) initial convergence factor
  nlinopts.etapow = 0.8;          % power law for reuse of estimate
  nlinopts.maxiter = 10;          % maximum number of iterations
  nlinopts.hfac = 0.5;            % for irregular convergence
elseif nlinauto == 2
  nlinopts.kappa = 10^-2;
  nlinopts.thetaREF = 0.3;
  nlinopts.thetaJAC = 0.3;
  nlinopts.thetaLU = 0.3;
  nlinopts.etainit = 0.95/(1-0.95);
  nlinopts.etapow = 0.7;
  nlinopts.maxiter = 9;
  nlinopts.hfac = 0.4;
  % (nlinauto >= 2 also uses a more conservative initial guess)
elseif nlinauto == 3
  nlinopts.kappa = 10^-3;
  nlinopts.thetaREF = 0.2;
  nlinopts.thetaJAC = 0.2;
  nlinopts.thetaLU = 0.2;
  nlinopts.etainit = 0.95/(1-0.95);
  nlinopts.etapow = 0.6;
  nlinopts.maxiter = 8;
  nlinopts.hfac = 0.3;
else
  error('Property ''nlinauto'' must be one of {1 2 3}.');
end

% solution and allocation logic
tspan = tspan(:);
nchunks = 10;
nchunksmax = 10000;
t = zeros(1,nchunks);
t(1) = tspan(1);
y = zeros(ndof,nchunks);
y(:,1) = yy;
i = 1;

% clear persistent data
l_implsolve(nlinopts);

% initial step
tt = tspan(1);
tend = tspan(end);
[h,TOL] = l_initialstep(rtol,atol,adaptivity,order, ...
                        NM,errnorm1,errnorm2,tt,yy,fun,varargin{:});
if initialstep > 0
  h = initialstep;
end
h = max(min(h,maxstep),minstep);
h2 = 0.5*h;

while 1
  if h2 < eps(16*(1+abs(tt)))
    warning(sprintf(['Unable to meet integration tolerances without ' ...
                     'reducing the step-size below the smallest ' ...
                     'value allowed. Failure at time t = %f.'],tt));
    t = t(1:i);
    y = y(:,1:i);
    break; % return solution obtained so far
  end

  % First take the full step y3 = T_1(h) using zero as the initial
  % approximation for the implicit part.
  [y3,rhon] = l_step(true,tt,yy,h,[], ...
                     TOL,NM,errnorm1,errnorm2,nlinopts,fun,varargin{:});
  if isempty(y3)
    h = l_limiter(rhon,1)*h; % retry with smaller h
    h2 = 0.5*h;
    cerr = 1;
    rho = 1;
    if info > 0
      disp('Nonlinear iteration converged poorly.');
    end
    continue;
    % *** the continue-statement implies that minstep+maxstep are
    % disregarded!
  end

  % Then take two half-steps T_1(h/2) and T_2(h/2) using the full step
  % T_1(h) as a predictor. These steps are usually easier than the
  % previous full step, so something strange is going on if this does
  % not converge.

  % y1 = T_1(h/2)
  if nlinauto < 3
    % interpolate using the full step
    [y1,rhon] = l_step(false,tt,yy,h2,0.5*(yy+y3), ...
                       TOL,NM,errnorm1,errnorm2,nlinopts,fun,varargin{:});
  else
    % be more suspicious...
    [y1,rhon] = l_step(false,tt,yy,h2,[], ...
                       TOL,NM,errnorm1,errnorm2,nlinopts,fun,varargin{:});
  end
  if isempty(y1)
    h2 = l_limiter(rhon,1)*h2;
    h = 2*h2; % (step just attempted was h2)
    cerr = 1;
    rho = 1;
    if info > 0
      warning('Nonlinear iteration did not converge!');
    end
    continue;
  end

  % y2 = T_2(h/2)
  if nlinauto == 1
    % predictor by combining the two previous steps
    [y2,rhon] = l_step(false,tt+h2,y1,h2,y1+0.5*(y3-yy), ...
                       TOL,NM,errnorm1,errnorm2,nlinopts,fun,varargin{:});
  elseif nlinauto == 2
    % using only the full step
    [y2,rhon] = l_step(false,tt+h2,y1,h2,y3, ...
                       TOL,NM,errnorm1,errnorm2,nlinopts,fun,varargin{:});
  else
    % suspicious...
    [y2,rhon] = l_step(false,tt+h2,y1,h2,[], ...
                       TOL,NM,errnorm1,errnorm2,nlinopts,fun,varargin{:});
  end
  if isempty(y2)
    h2 = l_limiter(rhon,1)*h2;
    h = 2*h2;
    cerr = 1;
    rho = 1;
    if info > 0
      warning('Nonlinear iteration did not converge!');
    end
    continue;
  end

  % Richardson's estimate...
  err = errnorm1((y3-y2)/3);

  % control objective
  switch adaptivity
   case 1 % EPS
    TOL = max(rtol.*errnorm1(y2),atol);
   case 2 % EPUS
    [dydt,dydtexpl] = fun(tt+h,y2,tt+h,y2,varargin{:});
    dydt = dydt+dydtexpl;
    TOL = h*max(rtol.*errnorm1(dydt),atol);
   case 3 % EQ
    TOL = max(rtol.*errnorm1(y2-yy),atol);
  end

  % digital control of step-size
  cerrold = cerr;
  cerr = l_limiter(1/(NM*errnorm2(err./TOL)^order),2);
  rho = l_limiter(l_lpfilter(ctrl,cerr,cerrold,rho),1);

  % acceptance/rejection based on smooth rho
  if rho >= rhomin
    % step accepted: record the value of T_2(h/2)
    tt = tt+h;
    yy = y2;
    i = i+1;
    if i > size(t,2)
      t = [t zeros(1,nchunks)];
      y = [y zeros(ndof,nchunks)];
      nchunks = min(2*nchunks,nchunksmax);
    end
    t(i) = tt;
    y(:,i) = yy;
    status = reportfun(tt,yy,[],varargin{:});

    % normal end of algorithm is here
    if tt-tend > -eps(32*abs(tt))
      status = reportfun(tt,yy,'done',varargin{:});
      t = t(1:i);
      y = y(:,1:i);
      break;
    elseif status ~= 0
      fprintf('Reporter returned non-zero status = %d. Bailing out.\n',status);
      status = reportfun(tt,yy,'done',varargin{:});
      t = t(1:i);
      y = y(:,1:i);
      break;
    end
  else
    if info > 0
      disp('Failed integration tolerance.');
      fprintf('   [h = %g, err = %f, cerr = %f, rho = %f]\n', ...
              h,1/(NM*errnorm2(err./TOL)^order),cerr,rho);
    end
  end
  h = rho*h;

  % in these two cases (not the preferred regime!) we must ensure that
  % the controller gets at least some kind of 'reasonable' information
  if h > maxstep
    rho = rho*maxstep/h;
    h = maxstep;
    cerr = 1;
  elseif h < minstep
    rho = rho*minstep/h;
    h = minstep;
    cerr = 1;
  end
  h2 = 0.5*h;
end

% output solution structure
if nargout < 2
  t = struct('x',t,'y',y,'solver','ode1s');
end

%-------------------------------------------------------------------------
function [h,TOL] = l_initialstep(rtol,atol,adaptivity,order, ...
                                 NM,errnorm1,errnorm2,t0,y0,fun,varargin)
%L_INITIALSTEP Starting step-size and tolerance.
%   [H,TOL] = L_INITIALSTEP(RTOL,ATOL,...) provides a starting
%   step-size whenever none is supplied. A value of the effective
%   tolerance TOL is also returned. It "usually gives a good guess for
%   the initial step-size... (or at least avoids a very bad choice)"
%   [1, p. 169]. The algorithm is described in [1, Chap. II.4,
%   p. 169].
%
% Reference:
%   [1] E. Hairer, S. P. Norsett and G. Wanner: "Solving Ordinary
%   Differential Equations I, Nonstiff Problems", 2nd ed, Springer
%   (1993).

% h0 by ensuring that an explicit Euler step is small
[dydt,dydtexpl] = fun(t0,y0,t0,y0,varargin{:});
dydt0 = dydt+dydtexpl;
% *** evaluation can be reused later

d0 = errnorm1(y0);
TOL = max(rtol.*d0,atol);
d0 = NM*errnorm2(d0./TOL);
d1 = NM*errnorm2(errnorm1(dydt0)./TOL);
if d0 >= 1e-5 && d1 >= 1e-5
  h0 = 0.01*d0/d1;
else
  h0 = 1e-6;
end

% take step and estimate the 2nd derivative
y1 = y0+h0*dydt0;
[dydt,dydtexpl] = fun(t0+h0,y1,t0+h0,y1,varargin{:});
dydt1 = dydt+dydtexpl;
d2 = NM*errnorm2(errnorm1(dydt1-dydt0)./(TOL*h0));

% update TOL (used in the implicit step)
switch adaptivity
 case 2 % EPUS
  TOL = h0*max(rtol.*errnorm1(dydt0),atol);
 case 3 % EQ
  TOL = max(rtol.*errnorm1(y1-y0),atol);
end

% h1 from simple error estimate
d3 = max(d1,d2);
if d3 > 1e-15
  h1 = (0.01/d3)^order;
else
  h1 = max(1e-6,h0*1e-3);
end

h = min(100*h0,h1);
% or otherwise integration might be aborted immediately:
h = max(h,eps(1024*(1+abs(t0))));

%-------------------------------------------------------------------------
function [y,rhon] = l_step(fullstep,t0,y0,h,yest,TOL, ...
                           NM,errnorm1,errnorm2, ...
                           nlinopts,fun,varargin)
%L_STEP One step with the split-step Euler method.
%   [Y,RHON] = L_STEP(FULLSTEP,T0,Y0,H,YEST,...) returns the step from
%   [T0,Y0] to [T0+H,Y] with a split-step Euler method using YEST as
%   an initial guess for the implicit solver.
%
%   See L_IMPLSOLVE.

if ~isempty(yest)
  ydelta = yest-y0;
else
  ydelta = zeros(size(y0));
end

% step: implicit step by solving a nonlinear system
[ydelta,rhon] = l_implsolve(fullstep,ydelta,t0,h,y0,TOL, ...
                            NM,errnorm1,errnorm2, ...
                            nlinopts,fun,varargin{:});
if isempty(ydelta)
  % nonlinear iteration did not converge, empty return
  y = [];
  return;
end

% assemble step
y = y0+ydelta;

%-------------------------------------------------------------------------
function [ydelta,rhon] = l_implsolve(fullstep,ydelta,t0,h,y0,TOL, ...
                                     NM,errnorm1,errnorm2, ...
                                     nlinopts,fun,varargin)
%L_IMPLSOLVE Solver for split-step Euler method.
%   [YDELTA,RHON] = L_IMPLSOLVE(FULLSTEP,YDELTA,T0,H,Y0,[TOL,...],
%   NLINOPTS,FUN,...) solves the equation YDELTA =
%   H*F(T0,Y0)+H*G(T0+H,Y0+YDELTA) for YDELTA. FULLSTEP is true if H
%   is a full step, false for half-steps. The split in
%   implicit/explicit parts is input as [G,F] =
%   FUN(T0+H,T0+YDELTA,T0,Y0,...). Empty YDELTA is returned whenever
%   poor convergence is detected, in which case RHON < 1 is a
%   suggested step-size reduction factor.
%
%   The implementation follows rather closely the prescription in [2,
%   Chap. IV.8, pp. 118--121], but has also been augmented with ideas
%   from [1].
%
% Reference:
%   [1] K. Gustafsson and G. Söderlind: "Control Strategies for the
%   Iterative Solution of Nonlinear Equations in ODE Solvers", SIAM
%   J. Sci. Comput. 18(1):23--40 (1997).
%   [2] E. Hairer and G. Wanner: "Solving Ordinary Differential
%   Equations II, Stiff and Differential-Algebraic Problems", 2nd ed,
%   Springer (1996).

% control of Jacobian evaluations and refactorizations
persistent J wksp tJAC L U P hLU eta theta L2 U2 P2;

if nargin == 1
  % clear persistent data
  nlinopts = fullstep;
  J = [];
  wksp = [];
  eta = nlinopts.etainit;
  theta = eta/(1+eta);
  return;
elseif isempty(eta)
  % previous attempt did not converge well
  eta = nlinopts.etainit;
  theta = eta/(1+eta);
end

% parameters of the solver
rhon = 1;                           % only changed in case of problems
eta = max(eta,eps)^nlinopts.etapow; % convergence factor (reused on success)
kappa = nlinopts.kappa;             % extra accuracy
maxiter = nlinopts.maxiter;         % max number of iterations

% The equation solved is now
%
%   ydelta = h*F(t0,y0)+h*G(t0+h,y0+ydelta,[t0 y0]),
%
% where [G,F] = fun(t0+h,y0+ydelta,t0,y0,varargin{:}) is the
% [implicit,explicit] part.

% first residual
[dydt,dydtexpl] = fun(t0+h,y0+ydelta,t0,y0,varargin{:});

% full/half-step
if fullstep
  % new Jacobian and/or factorization
  if isempty(J) || theta*hLU-abs(h-hLU) > nlinopts.thetaJAC*hLU
    % convergence was slow previously despite a not very large difference
    % in step-size (or there is no saved Jacobian)
    [J,wksp] = numjac(fun,t0+h,y0+ydelta,dydt,{0 1},wksp,0,[],[], ...
                      t0,y0,varargin{:});
    % *** user-defined Jacobian (exact and numerical)
    tJAC = t0; % (formed with data available up to t = t0)
    [L,U,P] = lu(eye(size(J))-h*J);
    hLU = h;
    P2 = []; % force refactorization at the half-step
  elseif isempty(P) || abs(h-hLU) > nlinopts.thetaLU*hLU
    % anticipated convergence failure due to a step-size change:
    % refactorize the Jacobian (or a new factorization was explicitly
    % called for)
    [L,U,P] = lu(eye(size(J))-h*J);
    hLU = h;
    P2 = [];
  end
else
  % when taking a half-step, the previous full step always
  % converged and hence we only care about the factorization
  if isempty(P2)
    [L2,U2,P2] = lu(eye(size(J))-h*J);
  end
  % During the iterations below we let the half-step change the
  % convergence indicators theta and eta since they should normally be
  % about the same. If this is not the case and the half-step is
  % rejected, then the full step is rejected automatically.
end
 
% first Newton-step
res = ydelta-h*dydtexpl-h*dydt;
if fullstep
  delta = U\(L\(P*res));
else
  delta = U2\(L2\(P2*res));
end
ydelta = ydelta-delta;
ndelta = NM*errnorm2(errnorm1(delta)./TOL);

% early return
if eta*ndelta <= kappa
  return;
end

for i = 2:maxiter
  % new update and residual
  res =  ydelta-h*dydtexpl-h*fun(t0+h,y0+ydelta,t0,y0,varargin{:});
  if fullstep
    delta = U\(L\(P*res));
  else
    delta = U2\(L2\(P2*res));
  end
  ydelta = ydelta-delta;

  % check convergence
  ndeltaold = ndelta;
  ndelta = NM*errnorm2(errnorm1(delta)./TOL);
  theta = ndelta/ndeltaold;
  eta = theta/(1-theta);
  
  if theta >= 1 || theta^(maxiter-i)*ndelta > kappa*(1-theta)
    % failure: diverging or not converging fast enough
    break;
  elseif eta*ndelta <= kappa
    % converged, return solution
    return;
  end
end

% iteration not terminated successfully, take action
if theta >= 1
  % diverging, suggest a smaller step
  rhon = nlinopts.thetaREF/theta;
elseif tJAC == t0
  % slow convergence despite a fresh Jacobian
  if theta > nlinopts.thetaREF
    % suggest a smaller step
    rhon = nlinopts.thetaREF/theta;
  else
    % rare: irregular convergence, shrink step by an arbitrary factor
    rhon = nlinopts.hfac;
  end
elseif ~fullstep
  % rare: half-step, irregular convergence despite a successful full
  % step?
  rhon = nlinopts.hfac;
end

% note: rejecting a half-step implies also a new full step and hence
% we postpone handling the Jacobian to after a smaller step has been
% tried
if fullstep
  % handle the Jacobian
  if tJAC ~= t0
    J = []; % force a new evaluation
  else
    P = []; % force refactorization
  end
end

ydelta = []; % step failed
eta = [];    % don't reuse convergence indicators 

%-------------------------------------------------------------------------
function u = l_limiter(u,kappa)
%L_LIMITER Step-size ratio limiter.
%   U = L_LIMITER(U,KAPPA) applies an arctangent limiter parametrized
%   by KAPPA to U.
%
% Reference:
%   [1] G. Söderlind, L. Wang: "Adaptive time-stepping and
%   computational stability", J. Comput. Appl. Math. 185:225--243 (2006).

u = 1+kappa*atan((u-1)/kappa);

%-------------------------------------------------------------------------
function rho = l_lpfilter(fir,cerr,cerrold,rho)
%L_LPFILTER Low-pass 2nd order digital filter.
%   RHO = L_LPFILTER(FIR,CERR,CERROLD,RHO) evaluates the digital
%   filter FIR = [k*beta_1 k*beta_2 alpha_2] using scaled controls
%   (CERR,CERROLD) and previous ratio RHO.
% 
%   The proposed next ratio is returned.
%
% Reference:
%   [1] G. Söderlind: "Digital Filters in Adaptive Time-Stepping",
%   ACM Trans. Math. Software, 29:1--26 (2003).

% immediate
rho = prod([cerr cerrold 1/rho].^fir);

%-------------------------------------------------------------------------
