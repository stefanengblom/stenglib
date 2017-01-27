function allok = runtest(idstr,funs,funnames,varargin)
%RUNTEST General test facility.
%   OK = RUNTEST(IDSTR,FUNS,FUNNAMES,...) calls the functions in the
%   cell-vector FUNS, one by one, and displays appropriate
%   messages. FUNNAMES is a cell-vector with corresponding names of
%   the functions and IDSTR is an identifying string which appears
%   first in all messages produces by RUNTEST (use empty to supress
%   all such messages). Additional arguments, if any, are passed to
%   the functions.
%
%   Each function in FUNS should be written so that it returns 1 for
%   'passed' and zero for 'failed'.
%
%   Example:
%     runtest('Medium tests (myfun)', ...
%       {@test1 @test2},{'Test #1' 'Test #2'});

% S. Engblom 2013-12-03 (empty IDSTR, output OK)
% S. Engblom 2012-05-25 (Additional arguments)
% S. Engblom 2004-10-29

% input
ntest = size(funs,2);
if size(funnames,2) ~= ntest, error('Input does not match.'); end

% run tests
nf = 0;
allok = 1;
for i = 1:size(funs,2)
  try
    ok = feval(funs{i},varargin{:});
    allok = allok && ok;
    if ~ok
      if ~isempty(idstr)
        fprintf(1,'%s: %s failed.\n',idstr,funnames{i});
      end
      nf = nf+1;
    else
      if ~isempty(idstr)
        fprintf(1,'%s: %s passed.\n',idstr,funnames{i});
      end
    end
  catch
    allok = 0;
    if ~isempty(idstr)
      fprintf(1,'%s: Error caught in %s.\n%s\n', ...
              idstr,funnames{i},lasterr);
    end
    nf = nf+1;
  end
end

% concluding diagnostics
if ~isempty(idstr)
  if nf == 0
    fprintf(1,'%s passed: %d test(s).\n',idstr,ntest);
   else
     fprintf(1,'%s failed: %d test(s) passed, %d failed.\n',idstr,ntest-nf,nf);
  end
end
