function [x,fx,flag,out,jac] = ainsolve(fun,x0,opts,varargin)
%AINSOLVE Solver for large sets of nonlinear equations.
%   AINSOLVE solves nonlinear equations of the form F(X) = 0, where F
%   and X are vectors, by using either a damped Newton method (DN) or
%   an accelerated Newton method (AN) where old search directions as
%   well as a preconditioned Krylov space is used together with the
%   Newton direction in order to accelerate the convergence. The
%   inversion of the Jacobian may be avoided using an inexact variant
%   (AIN) where the Newton direction is not used. Furthermore, there
%   are variants (XAN/XAIN) of both accelerated algorithms in which
%   the Krylov space is recomputed for every inner iteration. This is
%   more expensive but may also produce a faster convergence.
%
%   X = AINSOLVE(FUN,X0) starts at the initial guess X0 and tries to
%   solve the equations in FUN. FUN has the signature [Y,J] =
%   FUN(X,...) and returns a vector of values Y together with a
%   Jacobian J evaluated at the point X. Note that FUN will be called
%   with both one and two output arguments which can be exploited in
%   order to improve on the efficiency.
%
%   X = AINSOLVE(FUN,X0,OPTS) solves with the default parameters
%   replaced by the values in the structure OPTS. The syntax OPTDEF =
%   AINSOLVE returns the default options and the syntax OPTHELP =
%   AINSOLVE('opts') produces a short description of available
%   options.
%
%   X = AINSOLVE(FUN,X0,OPTS,P1,P2,...) passes the parameters P1, P2,
%   ... directly to the function FUN; FUN(X,P1,P2,...).
%
%   An important option is the use of the field OPTS.prefun which
%   specify the preconditioner for the Krylov space. The syntax works
%   as follows; initially let M = [] and for every assembled Jacobian
%   J, initialize the preonditioner by calling [PFUN,M] =
%   feval(OPTS.prefun,J,M,OPTS.prepar,P1,P2,...). Then the
%   preconditioner is defined by the call X = feval(PFUN,R,M), where X
%   should be an approximative solution to the equation J*X = R.
%
%   [X,FX] = AINSOLVE(FUN,X0,...) returns the value of the objective
%   function at X.
%
%   [X,FX,FLAG] = AINSOLVE(FUN,X0,...) returns a number FLAG that
%   describes the exit condition of AINSOLVE.
%   If FLAG is:
%      0 then AINSOLVE converged to a solution X,
%      1 then AINSOLVE converged to a solution X that does not satisfy
%        the tolerance on the norm of the residual,
%      2 then the maximum number of nonlinear iterations was reached,
%      3 then the maximum number of function evaluations was reached.
%
%   The error criteria used works as follows:
%     if all(abs(r) <= opts.atolf+opts.rtolf*abs(r0))
%       then exit with FLAG = 0,
%     elseif all(abs(x-xold) <= opts.atolx+opts.rtolx*abs(x))
%       then exit with FLAG = 1.
%   Here r0 is the initial residual and xold is the previous
%   approximation.
%
%   [X,FX,FLAG,OUT] = AINSOLVE(FUN,X0,...) returns a structure OUT
%   with the number of iterations taken in OUT.iterations, the number
%   of function (Jacobian) evaluations in OUT.funcCount (OUT.jacCount)
%   and the algorithm used in OUT.algorithm.
%
%   [X,FX,FLAG,OUT,JAC] = AINSOLVE(FUN,X0,...) also returns the
%   Jacobian JAC of FUN at X.
%
%   Cautionary: this is an experimental version with very little error
%   checking.
%
%   See also FSOLVE.

% S. Engblom 2010-09-03 (Minor revision)
% S. Engblom 2004-10-11 (Revision)
% S. Engblom 2004-01-26

% *** Use parseopts!
% *** Option jacmul
% *** Use of numjac
% *** For some combinations, only the action of J to a vector v is
% needed. This is just (F(u+h*v)-F(u))/h for a suitable small h.
% *** Check the angle between the latest step and the one just taken
% before adding it to the trace space? Alternative: update a
% Gram-Schmidt orthogonalization? How should the LS-problem really be
% solved? GMRES?

% special call: existing options
if nargin == 1 && ischar(fun) && strcmpi(fun,'opts')
  optdef.maxnlinit = 'integer > 0 {500}';
  optdef.maxfunevals = 'integer > 0 {1000}';
  optder.report = '{''on''} | ''off''';
  optdef.atolx = 'scalar {1e-8}';
  optdef.rtolx = 'scalar {1e-6}';
  optdef.atolf = 'scalar {1e-8}';
  optdef.rtolf = 'scalar {1e-6}';
  optdef.algorithm = ['''xain'' | {''ain''} | ' ...
                      '''xan'' | ''an'' | ''dn'''];
  optdef.linesearch = '{''bisect''} | ''poly''';
  optdef.tracedim = 'integer >= 0 {6}';
  optdef.krylovdim = 'integer >= 0 {6}';
  optdef.accmaxits = 'integer > 0 {5}';
  optdef.accminfac = 'scalar in (0,1) {0.95}';
  optdef.prefun = 'function {[]}';
  optdef.prepar = 'argument {[]}';
  x = optdef;
  return;
end

% default options
optdef.maxnlinit = 500;
optdef.maxfunevals = 1000;
optdef.report = 'on';
optdef.atolx = 1e-8;
optdef.rtolx = 1e-6;
optdef.atolf = 1e-8;
optdef.rtolf = 1e-6;
optdef.algorithm = 'ain';
optdef.linesearch = 'bisect';
optdef.tracedim = 6;
optdef.krylovdim = 6;
optdef.accmaxits = 5;
optdef.accminfac = 0.95;
optdef.prefun = [];
optdef.prepar = [];

% special call: return defaults
if nargin == 0, x = optdef; return; end

% merge defaults with actual inputs
if nargin > 2
  if iscell(opts), opts = struct(opts{:}); end
  fn = fieldnames(opts);
  for i = 1:length(fn)
    optdef = setfield(optdef,fn{i},getfield(opts,fn{i}));
  end
end
opts = optdef;

% size of the problem
ndof = size(x0,1);

% algorithm
switch opts.algorithm
 case 'dn', opts.alg = 1;
 case 'an', opts.alg = 2;
 case 'xan', opts.alg = 3;
 case 'ain', opts.alg = 4;
 case 'xain', opts.alg = 5;
 otherwise, error('Unknown algorithm.');
end
if opts.alg ~= 1, Vt = zeros(ndof,0); end

% linesearch
switch opts.linesearch
 case 'bisect', lines = @l_linesearch_bisect;
 case 'poly', lines = @l_linesearch_poly;
 otherwise, error('Unknown linesearch.');
end

% preconditioner
hasPrecond = ~isempty(opts.prefun);
if hasPrecond, M = []; end

% evaluate the first residual and Jacobian
[r,J] = feval(fun,x0,varargin{:});
nfun = 1; njac = 1;
rcmp = opts.atolf+opts.rtolf*abs(r);
normr = norm(r);

% diagnostics
opts.report = strcmpi(opts.report,'on');
if opts.report
  disp(sprintf([' Iteration  Func-count  Jac-count  ' ...
                ' |residual|     |step|\n' ...
                ' %5.0f      %5.0f      %5.0f       ' ...
                ' %0.6e   -'], ...
               0,nfun,njac,normr));
end

% early return
if all(abs(r) <= rcmp)
  x = x0;
  nliniter = 0;
else
  % nonlinear (outer) iterations
  for nliniter = 1:opts.maxnlinit

    if opts.alg == 1 % damped Newton (DN)
      en = J\(-r);
      [x,nfun] = l_linesearch(lines,x0,J,r,normr,en, ...
                              nfun,fun,varargin);

      % next residual and Jacobian
      [r,J] = feval(fun,x,varargin{:});
      nfun = nfun+1; njac = njac+1;
      normr = norm(r);
      step = x-x0;
    else % accelerated variants (AN/XAN/AIN/XAIN)
      % needed for progress report
      xold = x0;

      % compute the Newton direction (AN/XAN)
      if opts.alg == 2 || opts.alg == 3, en = J\(-r); end

      % inner iterations
      for j = 1:opts.accmaxits
        if j == 1 || opts.alg == 3 || opts.alg == 5
          % create Krylov space
          Vk = zeros(ndof,opts.krylovdim);
          if hasPrecond
            [prefun,M] = feval(opts.prefun,J,M,opts.prepar, ...
                               varargin{:});
            Vk(:,1) = feval(prefun,-r,M);
            for j = 2:opts.krylovdim
              Vk(:,j) = feval(prefun,J*Vk(:,j-1),M);
            end
          else
            Vk(:,1) = -r;
            for j = 2:opts.krylovdim
              Vk(:,j) = J*Vk(:,j-1);
            end
          end

          if opts.alg == 2 || opts.alg == 3 % (AN/XAN)
            W = orth([Vt en Vk]); % *** sloppy
          else % (AIN/XAIN)
            W = orth([Vt Vk]); % *** sloppy
          end
        end

        % restricted (local) Jacobian
        U = J*W;

        % Newton direction for local problem
        ea = W*(U\(-r)); % *** sloppy(?)
        [x,nfun] = l_linesearch(lines,x0,J,r,normr,ea, ...
                                nfun,fun,varargin);

        % next residual and Jacobian
        [r,J] = feval(fun,x,varargin{:});
        nfun = nfun+1; njac = njac+1;
        normr0 = normr;
        normr = norm(r);
 
        % should we stop?
        if all(abs(r) <= rcmp) || ...
              normr > opts.accminfac*normr0 || ...
              nfun >= opts.maxfunevals
          break;
        end
        x0 = x;
      end

      % restore old solution
      x0 = xold;

      % add the step thus produced to the trace space
      step = x-x0;
      if size(Vt,2) == opts.tracedim
        Vt = [Vt(:,2:end) step];
      else
        Vt = [Vt step];
      end
    end

    % diagnostics
    if opts.report
      disp(sprintf([' %5.0f      %5.0f      %5.0f       ' ...
                    ' %0.6e   %0.6e'], ...
                   nliniter,nfun,njac,normr,norm(step)));
    end

    % stopping criteria
    if  all(abs(r) <= rcmp) || ...
          all(abs(step) <= opts.atolx+opts.rtolx*abs(x)) || ...
          nfun >= opts.maxfunevals
      break;
    end

    % cycle solution for next round
    x0 = x;
  end
end

% various outputs
if nargout > 1
  fx = r;
  if nargout > 2
    flag = 0;
    if nfun >= opts.maxfunevals
      flag = 3;
    elseif nliniter >= opts.maxnlinit
      flag = 2;
    elseif any(abs(r) > rcmp)
      flag = 1;
    end
    if nargout > 3
      out = struct('algorithm',opts.algorithm, ...
                   'iterations',nliniter, ...
                   'funcCount',nfun, ...
                   'jacCount',njac);
      if nargout > 4
        jac = J;
      end
    end
  end
end

%--------------------------------------------------------------------------
function [x,fc] = l_linesearch(lines,xold,J,r,normr,eo,fcold,fun,fargs)
%L_LINESEARCH Linesearch algorithm.
%   [X,FC] = L_LINSEARCH(LINES,XOLD,J,R,NORMR,EO,FCOLD,FUN,FARGS)
%   performs an optimal scaling of the search direction EO and then
%   determines X = X0+LAMBDA*EO for which the norm of the function FUN
%   (with additional arguments FARGS) has decreased
%   sufficiently. LINES is a linesearch algorithm, J is the Jacobian
%   at XOLD, R and NORMR the residual and norm of the residual. FCOLD
%   is the function count before the call; it is increased for every
%   call to FUN and the accumulated number of calls is returned in FC.

% *** scaling really needed for AN/XAN/AIN/XAIN?
% *** factor 1/2 really needed for the norm? (kept for the sake of poly)

% straightforward
Je = J*eo;
mrtJe = -r'*Je;
normJe = norm(Je);
eta = mrtJe/(normr*normJe);
e = (mrtJe/(normJe*normJe))*eo;
grad = -(mrtJe/normJe)^2;

[x,fc] = feval(lines,xold,normr^2/2,grad,e, ...
               fcold,fun,fargs,0.5*eta*eta);
%--------------------------------------------------------------------------
function [x,fc] = l_linesearch_bisect(xold,fold,dfold,e, ...
                                      fcold,fun,fargs,eta2)
%L_LINESEARCH_BISECT Linesearch by bisection.

fc = fcold; % count the number of residuals evaluated
lambda = 1; % start by trying the full step

% backtrack by bisection
while lambda > eps
  x = xold+lambda*e;
  f = norm(feval(fun,x,fargs{:}))^2/2;
  fc = fc+1;
  if f <= (1-lambda*eta2)*fold, return; end
  lambda = 0.5*lambda;
end

% otherwise we are probably on a local minimum -- this should be
% handled by the caller by monitoring the residual
warning('Linesearch direction does not appear to be a descent.');
%--------------------------------------------------------------------------
function [x,fc] = l_linesearch_poly(xold,fold,dfold,e, ...
                                    fcold,fun,fargs,eta2)
%L_LINESEARCH_POLY Linesearch by polynomial interpolation.

% first try the full Newton step...
lambda1 = 1;
x = xold+e;
f1 = norm(feval(fun,x,fargs{:}))^2/2;
fc = fcold+1;

% this much decrease is sufficient
if f1 <= (1-eta2)*fold, return; end

% ...then backtrack using a quadratic form...
lambda2 = -dfold/(2*(f1-fold-dfold));
lambda2 = min(lambda2,0.5);
lambda2 = max(lambda2,0.1);
x = xold+lambda2*e;
f2 = norm(feval(fun,x,fargs{:}))^2/2;
fc = fc+1;
if f2 <= (1-lambda2*eta2)*fold, return; end

% ...otherwise backtrack using a cubic form
while lambda2 > eps
  % solve for certain cubic coefficients
  ab = ([1 -1; -lambda1 lambda2]* ...
        [(f2-dfold*lambda2-fold)/lambda2^2; ...
         (f1-dfold*lambda1-fold)/lambda1^2])./(lambda2-lambda1);

  % cycle the old values
  lambda1 = lambda2;
  f1 = f2;

  % compute next lambda safely
  if ab(1) ~= 0
    lambda2 = ab(2)^2-3*ab(1)*dfold;
    if lambda2 <= 0
      lambda2 = 1;
    elseif ab(2) <= 0
      lambda2 = (sqrt(lambda2)-ab(2))/(3*ab(1));
    else
      lambda2 = -dfold/(sqrt(lambda2)+ab(2));
    end
  elseif ab(2) ~= 0
    lambda2 = -dfold/(2*ab(2));
  else
    lambda2 = -fold/dfold;
  end

  % enforce bounds
  lambda2 = min(lambda2,0.5*lambda1);
  lambda2 = max(lambda2,0.1*lambda1);

  x = xold+lambda2*e;
  f2 = norm(feval(fun,x,fargs{:}))^2/2;
  fc = fc+1;
  if f2 <= (1-lambda2*eta2)*fold, return; end
end

% otherwise we are probably on a local minimum -- this should be
% handled by the caller by monitoring the residual
warning('Linesearch direction does not appear to be a descent.');
%--------------------------------------------------------------------------
