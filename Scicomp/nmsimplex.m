function [x,fx,flag] = nmsimplex(fun,x0,opts,varargin)
%NMSIMPLEX Nelder-Meads simplex-algorithm.
%   X = NMSIMPLEX(FUN,X0) starts at the initial guess X0 and tries to
%   minimize the objective function FUN. FUN has the signature Y =
%   FUN(X,...) and returns a scalar value.
%
%   X = NMSIMPLEX(FUN,X0,OPTS) solves with the default parameters replaced
%   by the values in the structure OPTS. The syntax OPTDEF = NMSIMPLEX
%   returns the default options and the syntax OPTHELP = NMSIMPLEX('opts')
%   produces a short description of available options (see also the table
%   below).
%
%   X = NMSIMPLEX(FUN,X0,OPTS,P1,P2,...) passes the parameters P1, P2,
%   ... directly to the function FUN; FUN(X,P1,P2,...).
%
%   [X,FX] = NMSIMPLEX(FUN,X0,...) returns the value of the objective
%   function at X.
%
%   [X,FX,FLAG] = NMSIMPLEX(FUN,X0,...) returns a number FLAG that
%   describes the exit condition of NMSIMPLEX.
%   If FLAG is:
%      0 then NMSIMPLEX converged to a solution X,
%      1 then NMSIMPLEX halted at a flat local minimum of the objective
%        function,
%      2 then the maximum number of function evaluations was
%        reached,
%      3 then the report function returned non-zero (see below).
%
%   The convergence criterion used works as follows:
%     if all(abs(xstep) <= opts.atolx+opts.rtolx*abs(xmin))
%       then exit with FLAG = 0,
%     elseif ystep <= opts.atolf+opts.rtolf*abs(ymin)
%       then exit with FLAG = 1.
%   Here (xmin,ymin) is the current best approximation and (xstep,ystep)
%   is the step just taken.
%
%   Property    Value/{Default}           Description
%   -----------------------------------------------------------------------
%   init        Scalar {0.25} | Simplex   A scalar indicating an
%                                         estimate of the error in the 
%                                         initial guess used to produce the
%                                         initial simplex. Alternatively, 
%                                         an ND-by-ND matrix where each
%                                         column defines a corner in the 
%                                         initial simplex. Together with 
%                                         the initial guess X0 this must 
%                                         define a feasible simplex.
%
%   maxfunevals integer > 0 {1000}        Maximum number of function 
%                                         evaluations.
%
%   report      {'on'} | 'off'            Print diagnostics.
%   reportfun   Function {[]}             Report function, see below.
%
%   atolx       Scalar {1e-8}             Absolute tolerance in solution.
%   rtolx       Scalar {1e-6}             Relative tolerance in solution.
%
%   atolf       Scalar {1e-14}            Absolute tolerance in function.
%   rtolf       Scalar {1e-12}            Relative tolerance in function.
%
%   lb          Scalar | Vector {-Inf}    Lower bound in each dimension.
%   ub          Scalar | Vector {Inf}     Upper bound in each dimension.
%
%   The report function is called according to the three syntaxes STATUS =
%   reportfun([t(0) -1],suggested title,'init'), STATUS =
%   reportfun(t,xsimplex,[],...), and, STATUS =
%   reportfun(t,xsimplex,'done',...). The "time" is defined by t =
%   -max(abs(xstep)./(opts.atolx+opts.rtolx*abs(xmin))). An occasional
%   useful feature is that the minimization is automatically aborted
%   whenever the returned value STATUS is non-zero.
%
%   Example:
%     % Rosenbrock function
%     f = @(x)((1-x(1,:)).^2+100*(x(2,:)-x(1,:).^2).^2);
%     [x,y] = meshgrid(linspace(-2,2,40),linspace(-1,3));
%     figure, contour(x,y,reshape(f([x(:)'; y(:)']),size(x)),30),
%     hold on
%
%     % reporter to plot convergence (could also be an .m-file function!)
%     iif = @(varargin)(varargin{2*find([varargin{1:2:end}],1,'first')}());
%     simplex_report = @(t,x,s)(iif(~ischar(x), ...
%       @()0*plot(x(1,[1:end 1]),x(2,[1:end 1]),'b.-'), ...
%       true,@()0));
%
%     % solve
%     x = nmsimplex(f,[-1 2],{'reportfun' simplex_report})
%
%   See also FMINSEARCH, FSOLVE, REPORT.

% S. Engblom 2014-08-28 (Revision, reportfun and plot example)
% S. Engblom 2006-09-19 (Revision)
% S. Engblom 2006-02-03

% *** stat counter of some kind?

% special call: existing options
if nargin == 1 & ischar(fun) & strcmpi(fun,'opts')
  optdef.init = 'scalar | simplex {0.25}'
  optdef.maxfunevals = 'integer > 0 {1000}';
  optdef.report = '{''on''} | ''off''';
  optdef.reportfun = 'Function {[]}';
  optdef.atolx = 'scalar {1e-8}';
  optdef.rtolx = 'scalar {1e-6}';
  optdef.atolf = 'scalar {1e-14}';
  optdef.rtolf = 'scalar {1e-12}';
  optdef.lb = 'scalar | vector {-inf}';
  optdef.ub = 'scalar | vector {inf}';
  x = optdef;
  return;
end

% default options
optdef.init = 0.25;
optdef.maxfunevals = 1000;
optdef.report = 'on';
optdef.reportfun = [];
optdef.atolx = 1e-8;
optdef.rtolx = 1e-6;
optdef.atolf = 1e-14;
optdef.rtolf = 1e-12;
optdef.lb = -inf;
optdef.ub = inf;

% special call: return defaults
if nargin == 0, x = optdef; return; end

% parse options
if nargin < 3
  opts = optdef;
else
  opts = parseopts(optdef,opts);
end

% initial guess
x0 = x0(:);
nd = size(x0,1);
opts.lb = opts.lb(:);
if isscalar(opts.lb)
  opts.lb = frepmat(opts.lb,nd);
end
opts.ub = opts.ub(:);
if isscalar(opts.ub)
  opts.ub = frepmat(opts.ub,nd);
end
if any(x0 < opts.lb) || any(opts.ub < x0)
  error('Initial guess not feasible.');
end

nfun = 0;
if isscalar(opts.init)
  % try creating a non-degenerate simplex
  r = opts.init*abs(x0)+opts.atolx;
  x1 = tsum(x0,diag(r),[1],[1 2]);
  y1 = zeros(1,size(x1,2));
  for i = 1:size(x1,2)
    x1(:,i) = min(max(x1(:,i),opts.lb),opts.ub);
    y1(i) = feval(fun,x1(:,i),varargin{:});
  end
  nfun = nfun+size(x1,2);

  x2 = tsum(x0,diag(-r),[1],[1 2]);
  y2 = zeros(1,size(x2,2));
  for i = 1:size(x2,2)
    x2(:,i) = min(max(x2(:,i),opts.lb),opts.ub);
    y2(i) = feval(fun,x2(:,i),varargin{:});
  end
  nfun = nfun+size(x2,2);

  % use best corners for the initial guess
  [foo,i] = min([y1; y2]);
  x = [x1(:,i == 1) x2(:,i == 2)];
  y = [y1(i == 1) y2(i == 2)];

  % add user's corner...
  x = [x x0];
  y = [y feval(fun,x0,varargin{:})];
  nfun = nfun+1;
else
  % intial simplex supplied by user
  x = [x0 opts.init];
  y = zeros(1,size(x,2));
  for i = 1:size(x,2)
    x(:,i) = min(max(x(:,i),opts.lb),opts.ub);
    y(i) = feval(fun,x(:,i),varargin{:});
  end
  nfun = nfun+size(x,2);
end

if size(x,2) ~= nd+1 || rank(tsum(x(:,2:end),-x(:,1),[1 2],[1])) < nd
  error(['Degenerate initial simplex. Either provide a feasible initial ' ...
         'simplex or try a smaller scalar option .init.']);
end

opts.report = strcmpi(opts.report,'on');
if opts.report
  disp(' Iteration  Func-count   Func-value      |step|');
end

% find best and worst values
[ymin,imin] = min(y);
[ymax,imax] = max(y);
xmin = x(:,imin);
xmax = x(:,imax);
ymax2 = max(y([1:imax-1 imax+1:end]));

% setup reporter
if ischar(opts.reportfun)
  reportfun = str2func(opts.reportfun);
else
  reportfun = opts.reportfun;
end
if ~isempty(reportfun)
  status = reportfun([-max(abs(xmax-xmin)./(opts.atolx+opts.rtolx*abs(xmin))) -1], ... 
                     'Solution progress...','init');
else
  % void function
  reportfun = @report;
  status = reportfun([],[],'none');
end

nit = 0;
while 1
  if status ~= 0
    fprintf('Reporter returned non-zero status = %d. Bailing out.\n',status);
    status = reportfun(-max(abs(xmax-xmin)./(opts.atolx+opts.rtolx*abs(xmin))), ...
                       x,'done',varargin{:});
    x = xmin;
    if nargout > 1
      fx = ymin;
      if nargout > 2
        flag = 3;
      end
    end
    return;
  end
  nit = nit+1;  

  xcen = (-xmax+sum(x,2))/nd;
  xref = min(max(2*xcen-xmax,opts.lb),opts.ub);
  yref = feval(fun,xref,varargin{:});
  nfun = nfun+1;

  % form new simplex
  if ymin > yref
    % xref seems promising; try to extrapolate
    xexp = min(max(2*xref-xcen,opts.lb),opts.ub);
    yexp = feval(fun,xexp,varargin{:});
    nfun = nfun+1;
    if yexp < yref
      % successful, good decrease
      x(:,imax) = xexp;
      y(imax) = yexp;
    else
      % a fair decrease anyway...
      x(:,imax) = xref;
      y(imax) = yref;
    end
  elseif ymax2 > yref
    % intermediate cost; slight improvement
    x(:,imax) = xref;
    y(imax) = yref;
  else
    % high cost; local contraction
    if ymax <= yref
      x(:,imax) = 0.5*(xmax+xcen);
    else
      x(:,imax) = 0.5*(xref+xcen);
    end

    % new evaluation
    y(imax) = feval(fun,x(:,imax),varargin{:});
    nfun = nfun+1;

    % still not an improvement: try multiple contraction around the
    % best corner
    if y(imax) > ymax
      for i = [1:imin-1 imin+1:size(x,2)]
        x(:,i) = 0.5*(x(:,i)+xmin);
        y(i) = feval(fun,x(:,i),varargin{:});
      end
      nfun = nfun+nd;
    end
  end

  % update best and worst values
  [ymin,imin] = min(y);
  [ymax,imax] = max(y);
  xmin = x(:,imin);
  xmax = x(:,imax);
  ymax2 = max(y([1:imax-1 imax+1:end]));

  % report progress
  xstep = xmax-xmin;
  ystep = ymax-ymin;
  if opts.report
    fprintf(1,' %5.0f       %5.0f       %0.6e    %0.6e\n', ...
            nit,nfun,ymin,norm(xstep));
  end
  status = reportfun(-max(abs(xstep)./(opts.atolx+opts.rtolx*abs(xmin))),x,'');

  % convergence criteria
  ok = all(abs(xstep) <= opts.atolx+opts.rtolx*abs(x(:,imin)));
  if ok || abs(ystep) <= opts.atolf+opts.rtolf*abs(ymin) || ...
        nfun >= opts.maxfunevals
    break;
  end
end

% various outputs
x = xmin;
if nargout > 1
  fx = ymin;
  if nargout > 2
    flag = 0;
    if nfun >= opts.maxfunevals
      flag = 2;
    elseif ~ok
      flag = 1;
    end
  end
elseif nargout < 3
  if nfun >= opts.maxfunevals
    warning('The maximum number of function evaluations was reached.');
  elseif ~ok
    warning('Halted at a flat local minimum of the objective function.');
  end
end

status = reportfun(-max(abs(xstep)./(opts.atolx+opts.rtolx*abs(xmin))),x,'done');
