function status = report(t,title,s,varargin)
%REPORT Report progress of solver.
%   The function REPORT is designed to report progress of any type of solver
%   and for Matlab's ODE-solvers in particular.
%
%   Only one reporter may be active at any given time. Use WAITBAR
%   directly to achieve multiple reporters.
%
%   STATUS = REPORT([T0 Tend],TITLE,'init',...) sets up the progress
%   reporter. T0 is the initial time/initial accuracy and Tend the final
%   time/requested accuracy. The character array TITLE is the title of the 
%   reporter; if TITLE is a non-character array or empty, the default used 
%   is 'Solution progress'. This somewhat singular behavior is due to 
%   considerations of compatibility with Matlab's ODE-solvers.
%
%   STATUS = REPORT([T0 Tend],'timeleft','init',...) continuously
%   measures time and reports an estimate of the time until
%   completion.
%
%   STATUS = REPORT([],[],'none') causes the reporter to be inactive until
%   reinitialized again. This is useful for avoiding if-clauses before
%   calling REPORT.
%
%   STATUS = REPORT(T,[],'',...) reports the progress T, where T0 <= T <=
%   Tend.
%
%   STATUS = REPORT([],[],'',...) relies on the 'timeleft'-syntax and
%   updates the estimated time without computing a new estimate. This
%   is useful when the number of tasks are limited, but each task
%   takes a considerable amount of time.
%
%   STATUS = REPORT([],[],'done',...) finishes the reporter.
%
%   STATUS = 0 is returned in all syntaxes.
%
%   See also WAITBAR.

% S. Engblom 2009-03-06 (Minor revision)
% S. Engblom 2007-04-05

persistent T0 Tend percent hwait t0 est lastupdate;

if isempty(s) && ~isempty(hwait)
  if isempty(t)
    if ~isempty(est) && ~isempty(lastupdate)
      % blind update (relying on est and lastupdate being properly initiated)
      est = est-toc(lastupdate);
      hrs = floor(est/3600);
      mins = floor(est/60)-60*hrs;
      waitbar(percent/100,hwait, ...
              sprintf('Estimated time left: %d hour(s), %d minute(s).', ...
                      hrs,mins));
      lastupdate = tic;
    end
    % call ignored otherwise
  else
    % main use
    now = round((t(end)-T0)/(Tend-T0)*100);
    if now > percent
      percent = now;
      if isempty(t0) || percent < eps
        waitbar(percent/100,hwait);
      else
        t1 = toc(t0);
        est = (100-percent)/percent*t1;
        hrs = floor(est/3600);
        mins = floor(est/60)-60*hrs;
        waitbar(percent/100,hwait, ...
                sprintf('Estimated time left: %d hour(s), %d minute(s).', ...
                        hrs,mins));
      end
    end
    lastupdate = tic;
  end
elseif strcmp(s,'init')
  T0 = t(1);
  Tend = t(end);
  t0 = [];
  est = [];
  lastupdate = [];
  percent = 0;
  try, close(hwait); catch, end
  if isempty(title) || ~isa(title,'char')
    title = 'Solution progress';
  elseif strcmp(title,'timeleft')
    title = 'Estimated time left: ';
    t0 = tic;
  end
  hwait = waitbar(0.0,title);
else % 'none' and 'done' empties hwait
  try, close(hwait); catch, end
  hwait = [];
  t0 = [];
  est = [];
  lastupdate = [];
end

status = 0;
