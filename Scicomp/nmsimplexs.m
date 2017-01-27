function [x,fx,flag] = nmsimplexs(fun,x0,opts,varargin)
%NMSIMPLEXS Nelder-Meads simplex-algorithm, stochastic version.
%   X = NMSIMPLEXS(FUN,X0) starts at the initial guess X0 and tries to
%   minimize the expected value of the objective function FUN. FUN has the
%   signature Y = FUN(X,...) and returns a scalar value.
%
%   NMSIMPLEXS is designed with stochastic functions in mind where the value
%   of FUN can vary between calls. The algorithm tries to find a solution X
%   up to a given tolerance such that the interval MEAN(FUN(X,...))  +/-
%   alpha*STD(FUN(X,...))/SQRT(N) is minimized for some given precision
%   alpha, with N the number of evaluations of FUN in X.
%
%   X = NMSIMPLEXS(FUN,X0,OPTS) solves with the default parameters replaced
%   by the values in the structure OPTS. The syntax OPTDEF = NMSIMPLEXS
%   returns the default options and the syntax OPTHELP = NMSIMPLEXS('opts')
%   produces a short description of available options (see also the table
%   below).
%
%   X = NMSIMPLEXS(FUN,X0,OPTS,P1,P2,...) passes the parameters P1, P2,
%   ... directly to the function FUN; FUN(X,P1,P2,...).
%
%   [X,FX] = NMSIMPLEXS(FUN,X0,...) returns the value of the objective
%   function at X.
%
%   [X,FX,FLAG] = NMSIMPLEXS(FUN,X0,...) returns a number FLAG that
%   describes the exit condition of NMSIMPLEXS.
%   If FLAG is:
%      0 then NMSIMPLEXS converged to a solution X,
%      1 then NMSIMPLEXS halted at a flat local minimum of the objective
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
%   alpha       Scalar {2}                Stochastic precision in solution.
%   nmin        Scalar {6}                Minimum number of function
%                                         evaluations in any simplex corner.
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
%     % a stochastic Rosenbrock function
%     f = @(x)((1-x(1,:)).^2+100*(x(2,:)-x(1,:).^2).^2+randn(1,size(x,2)));
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
%     x = nmsimplexs(f,[-1 2],{'reportfun' simplex_report 'alpha' 2 ...
%       'maxfunevals' 1e5 'nmin' 5})
%
%   See also NMSIMPLEX, FMINSEARCH, FSOLVE, REPORT.

% S. Engblom 2014-08-30 (Initial version, built from NMSIMPLEX)

% *** stat counter of some kind?

% special call: existing options
if nargin == 1 & ischar(fun) & strcmpi(fun,'opts')
  optdef.init = 'scalar | simplex {0.25}'
  optdef.maxfunevals = 'integer > 0 {1000}';
  optdef.report = '{''on''} | ''off''';
  optdef.reportfun = 'Function {[]}';
  optdef.alpha = 'scalar {2}';
  optdef.nmin = 'scalar {6}';
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
optdef.nmin = 6;
optdef.alpha = 2;
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

% *** assignopts?
nmin = opts.nmin;

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
  % ***
  warning('Syntax implemented as a hack.');
  
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
  y = cell(1,size(x,2));
  for i = 1:size(x,2), y{i} = zeros(1,0); end
  my = zeros(1,size(x,2));
  sy = zeros(1,size(x,2));
else
  % intial simplex supplied by user
  x = [x0 opts.init];
  y = cell(1,size(x,2));
  for i = 1:size(x,2)
    x(:,i) = min(max(x(:,i),opts.lb),opts.ub);
    y{i} = zeros(1,0);
  end
  my = zeros(1,size(x,2));
  sy = zeros(1,size(x,2));
end

if size(x,2) ~= nd+1 || rank(tsum(x(:,2:end),-x(:,1),[1 2],[1])) < nd
  error(['Degenerate initial simplex. Either provide a feasible initial ' ...
         'simplex or try a smaller scalar option .init.']);
end

opts.report = strcmpi(opts.report,'on');
if opts.report
  disp(' Iteration  Func-count   Func-value      |step|');
end

% initial evaluation
for i = 1:size(x,2)
  for j = 1:nmin, y{i} = [y{i} feval(fun,x(:,i),varargin{:})]; end
  my(i) = mean(y{i});
  sy(i) = std(y{i})/sqrt(nmin);
end
nfun = nfun+nmin*size(x,2);

% determine the topology [y(imin) <= <bulk> <= y(imax2) <= y(imax)] up to
% precision alpha, adding function evaluations as needed
while 1 % (really a repeat...until)
  % current topology
  [foo,is] = sort(my);
  imax = is(end);
  imax2 = is(end-1);
  imin = is(1);
  alpha0 = [(my(imax2)-my(imin))/(sy(imax2)+sy(imin)) ...
            (my(imax)-my(imax2))/(sy(imax)+sy(imax2)) ...
            (my(is(2:end-2))-my(imin))./(sy(is(2:end-2))+sy(imin)) ...
            (my(imax2)-my(is(2:end-2)))./(sy(imax2)+sy(is(2:end-2)))];

  if all(alpha0 >= opts.alpha)
    break; % done when all alpha0's are large enough
  else
    % greedily improve alpha by adding evaluations, assuming for convenience
    % that the topology does not change during the current round
    if alpha0(1) < opts.alpha
      % (imax2,imin) cannot be compared-alpha
      if sy(imax2) >= sy(imin), i = imax2; else, i = imin; end
      y{i} = [y{i} feval(fun,x(:,i),varargin{:})];
      my(i) = mean(y{i});
      sy(i) = std(y{i})/sqrt(size(y{i},2));
      nfun = nfun+1;
    end
    if alpha0(2) < opts.alpha
      % (imax,imax2) cannot be compared-alpha
      if sy(imax) >= sy(imax2), i = imax; else, i = imax2; end
      y{i} = [y{i} feval(fun,x(:,i),varargin{:})];
      my(i) = mean(y{i});
      sy(i) = std(y{i})/sqrt(size(y{i},2));
      nfun = nfun+1;
    end
    for ia = find(alpha0(3:nd) > opts.alpha)
      % (imin,is(ia+1)) cannot be compared-alpha
      if sy(imin) >= sy(is(ia+1)), i = imin; else, i = is(ia+1); end
      y{i} = [y{i} feval(fun,x(:,i),varargin{:})];
      my(i) = mean(y{i});
      sy(i) = std(y{i})/sqrt(size(y{i},2));
      nfun = nfun+1;
    end
    for ia = find(alpha0(nd+1:end) > opts.alpha)
      % (imax2,is(ia+1)) cannot be compared-alpha
      if sy(imax2) >= sy(is(ia+1)), i = imax2; else, i = is(ia+1); end
      y{i} = [y{i} feval(fun,x(:,i),varargin{:})];
      my(i) = mean(y{i});
      sy(i) = std(y{i})/sqrt(size(y{i},2));
      nfun = nfun+1;
    end
  end
end

xmin = x(:,imin);
xmax = x(:,imax);

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
      fx = my(imin);
      if nargout > 2
        flag = 3;
      end
    end
    return;
  end
  nit = nit+1;  

  % center point of nd best corners
  xcen = (-xmax+sum(x,2))/nd;

  % reflection away from worst corner
  xref = min(max(2*xcen-xmax,opts.lb),opts.ub);
  yref = zeros(1,0);
  for j = 1:nmin, yref = [yref feval(fun,xref,varargin{:})]; end
  myref = mean(yref);
  syref = std(yref)/sqrt(nmin);
  nfun = nfun+nmin;

  % ensure the reflection point can be compared with the best corner to within
  % precision alpha
  alpha1 = abs(my(imin)-myref)/(sy(imin)+syref);
  while alpha1 < opts.alpha
    yref = [yref feval(fun,xref,varargin{:})];
    myref = mean(yref);
    syref = std(yref)/sqrt(size(yref,2));
    nfun = nfun+1;
    if syref <= sy(imin)
      i = imin;
      y{i} = [y{i} feval(fun,x(:,i),varargin{:})];
      my(i) = mean(y{i});
      sy(i) = std(y{i})/sqrt(size(y{i},2));
      nfun = nfun+1;
    end
    alpha1 = abs(my(imin)-myref)/(sy(imin)+syref);
  end

  % form new simplex
  if my(imin) > myref
    % xref seems promising; try to extrapolate
    xexp = min(max(2*xref-xcen,opts.lb),opts.ub);
    yexp = zeros(1,0);
    for j = 1:nmin, yexp = [yexp feval(fun,xexp,varargin{:})]; end
    myexp = mean(yexp);
    syexp = std(yexp)/sqrt(nmin);
    nfun = nfun+nmin;
  
    % ensure the extrapolation point can be compared with the reflection point
    % to within precision alpha
    alpha2 = abs(myref-myexp)/(syref+syexp);
    while alpha2 < opts.alpha
      yexp = [yexp feval(fun,xexp,varargin{:})];
      myexp = mean(yexp);
      syexp = std(yexp)/sqrt(size(yexp,2));
      nfun = nfun+1;
      if syexp <= syref
        yref = [yref feval(fun,xref,varargin{:})];
        myref = mean(yref);
        syref = std(yref)/sqrt(size(yref,2));
        nfun = nfun+1;
      end
      alpha2 = abs(myref-myexp)/(syref+syexp);
    end
    
    if myexp < myref
      % successful, good decrease
      x(:,imax) = xexp;
      y{imax} = yexp;
      my(imax) = myexp;
      sy(imax) = syexp;
    else
      % a fair decrease anyway...
      x(:,imax) = xref;
      y{imax} = yref;
      my(imax) = myref;
      sy(imax) = syref;
    end
  else
    % ensure the reflection point can be compared with the second worst point to
    % within precision alpha
    alpha3 = abs(myref-my(imax2))/(syref+sy(imax2));
    while alpha3 < opts.alpha
      yref = [yref feval(fun,xref,varargin{:})];
      myref = mean(yref);
      syref = std(yref)/sqrt(size(yref,2));
      nfun = nfun+1;
      if syref <= sy(imax2)
        i = imax2;
        y{i} = [y{i} feval(fun,x(:,i),varargin{:})];
        my(i) = mean(y{i});
        sy(i) = std(y{i})/sqrt(size(y{i},2));
        nfun = nfun+1;
      end
      alpha3 = abs(myref-my(imax2))/(syref+sy(imax2));
    end
    
    if my(imax2) > myref
      % intermediate cost; slight improvement
      x(:,imax) = xref;
      y{imax} = yref;
      my(imax) = myref;
      sy(imax) = syref;
    else
      % ensure the reflection point can be compared with the worst point to within
      % precision alpha
      alpha4 = abs(myref-my(imax))/(syref+sy(imax));
      while alpha4 < opts.alpha
        yref = [yref feval(fun,xref,varargin{:})];
        myref = mean(yref);
        syref = std(yref)/sqrt(size(yref,2));
        nfun = nfun+1;
        if syref <= sy(imax)
          i = imax;
          y{i} = [y{i} feval(fun,x(:,i),varargin{:})];
          my(i) = mean(y{i});
          sy(i) = std(y{i})/sqrt(size(y{i},2));
          nfun = nfun+1;
        end
        alpha4 = abs(myref-my(imax))/(syref+sy(imax));
      end

      % high cost; local contraction
      if my(imax) <= myref
        xcontr = 0.5*(xmax+xcen);
      else
        xcontr = 0.5*(xref+xcen);
      end

      % new evaluation
      ycontr = zeros(1,0);
      for j = 1:nmin, ycontr = [ycontr feval(fun,xcontr,varargin{:})]; end
      mycontr = mean(ycontr);
      sycontr = std(ycontr)/sqrt(nmin);
      nfun = nfun+nmin;

      % ensure the contraction point can be compared with the worst point to
      % within precision alpha
      alpha5 = abs(mycontr-my(imax))/(sycontr+sy(imax));
      while alpha5 < opts.alpha
        ycontr = [ycontr feval(fun,xcontr,varargin{:})];
        mycontr = mean(ycontr);
        sycontr = std(ycontr)/sqrt(size(ycontr,2));
        nfun = nfun+1;
        if sycontr <= sy(imax)
          i = imax;
          y{i} = [y{i} feval(fun,x(:,i),varargin{:})];
          my(i) = mean(y{i});
          sy(i) = std(y{i})/sqrt(size(y{i},2));
          nfun = nfun+1;
        end
        alpha5 = abs(mycontr-my(imax))/(sycontr+sy(imax));
      end

      if mycontr <= my(imax)
        % an improvement now?
        x(:,imax) = xcontr;
        y{imax} = ycontr;
        my(imax) = mycontr;
        sy(imax) = sycontr;
      else
        % still not an improvement: try multiple contractions around the best corner
        for i = [1:imin-1 imin+1:size(x,2)]
          x(:,i) = 0.9*x(:,i)+0.1*xmin; % (0.9 instead of 0.5 is the RS9 modification)
          y{i} = zeros(1,0);
          for j = 1:nmin, y{i} = [y{i} feval(fun,x(:,i),varargin{:})]; end
          my(i) = mean(y{i});
          sy(i) = std(y{i})/sqrt(nmin);
        end
        nfun = nfun+nmin*nd;
      end
    end
  end

  % update best and worst values

  % determine the topology [y(imin) <= <bulk> <= y(imax2) <= y(imax)] up to
  % precision alpha, adding function evaluations as needed
  while 1 % (really a repeat...until)
    % current topology
    [foo,is] = sort(my);
    imax = is(end);
    imax2 = is(end-1);
    imin = is(1);
    alpha0 = [(my(imax2)-my(imin))/(sy(imax2)+sy(imin)) ...
              (my(imax)-my(imax2))/(sy(imax)+sy(imax2)) ...
              (my(is(2:end-2))-my(imin))./(sy(is(2:end-2))+sy(imin)) ...
              (my(imax2)-my(is(2:end-2)))./(sy(imax2)+sy(is(2:end-2)))];
    
    if all(alpha0 >= opts.alpha)
      break; % done when all alpha0's are large enough
    else
      % greedily improve alpha by adding evaluations, assuming for convenience
      % that the topology does not change during the current round
      if alpha0(1) < opts.alpha
        % (imax2,imin) cannot be compared-alpha
        if sy(imax2) >= sy(imin), i = imax2; else, i = imin; end
        y{i} = [y{i} feval(fun,x(:,i),varargin{:})];
        my(i) = mean(y{i});
        sy(i) = std(y{i})/sqrt(size(y{i},2));
        nfun = nfun+1;
      end
      if alpha0(2) < opts.alpha
        % (imax,imax2) cannot be compared-alpha
        if sy(imax) >= sy(imax2), i = imax; else, i = imax2; end
        y{i} = [y{i} feval(fun,x(:,i),varargin{:})];
        my(i) = mean(y{i});
        sy(i) = std(y{i})/sqrt(size(y{i},2));
        nfun = nfun+1;
      end
      for ia = find(alpha0(3:nd) > opts.alpha)
        % (imin,is(ia+1)) cannot be compared-alpha
        if sy(imin) >= sy(is(ia+1)), i = imin; else, i = is(ia+1); end
        y{i} = [y{i} feval(fun,x(:,i),varargin{:})];
        my(i) = mean(y{i});
        sy(i) = std(y{i})/sqrt(size(y{i},2));
        nfun = nfun+1;
      end
      for ia = find(alpha0(nd+1:end) > opts.alpha)
        % (imax2,is(ia+1)) cannot be compared-alpha
        if sy(imax2) >= sy(is(ia+1)), i = imax2; else, i = is(ia+1); end
        y{i} = [y{i} feval(fun,x(:,i),varargin{:})];
        my(i) = mean(y{i});
        sy(i) = std(y{i})/sqrt(size(y{i},2));
        nfun = nfun+1;
      end
    end
  end

  xmin = x(:,imin);
  xmax = x(:,imax);

  % report progress
  xstep = xmax-xmin;
  ystep = my(imax)-my(imin);
  if opts.report
    fprintf(1,' %5.0f       %5.0f       %0.6e    %0.6e\n', ...
            nit,nfun,my(imin),norm(xstep));
  end
  status = reportfun(-max(abs(xstep)./(opts.atolx+opts.rtolx*abs(xmin))),x,'');

  % convergence criteria
  ok = all(abs(xstep) <= opts.atolx+opts.rtolx*abs(x(:,imin)));
  if ok || abs(ystep) <= opts.atolf+opts.rtolf*abs(my(imin)) || ...
        nfun >= opts.maxfunevals
    break;
  end
end

% various outputs
x = xmin;
if nargout > 1
  fx = my(imin);
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
