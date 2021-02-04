function S = arr2latex(T,fmt,varargin)
%ARR2LATEX LaTeX-table from matrix.
%   S = ARR2LATEX(T,FMT) creates a LaTeX-table S from a matrix T and a
%   format string FMT.
%
%   T must be a matrix containing scalars that SPRINTF
%   accepts. Currently, complex matrices are not supported.
%
%   FMT must be a cell-matrix of the same size as T. Singleton
%   dimensions are, however, automatically expanded. ARR2LATEX
%   supports all format-flags as supported by SPRINTF but the
%   percent-character should not be included as it is added
%   automatically.
%
%   Additionally, flags ending with a dollar-sign (see the second
%   example below) expands the 'scientific' notation so that
%   '6.022e23' becomes '$6.022 \cdot 10^{23}$' instead. Also, 'n_'
%   reflects n whitespaces in a row and is useful for typesetting
%   empty entries.
%
%   S = ARR2LATEX(T,FMT,...) allows for options to be passed; the
%   table below explains the available options.
%
%   Property    Value/{Default}           Description
%   -----------------------------------------------------------------------
%   collabel,   Cell-vector with          Label of columns and rows.
%   rowlabel    strings {''}
%
%   breakrow    Integer {74}              Where to insert a newline.
%
%   hline       {'on'} | 'off'            Insert a \hline between
%                                         each row.
%
%   posinf      String {'$\infty$'}       Symbols of plus/minus infinity
%   neginf      String {'$-\infty$'}      and NaN's. Use empty to type
%   nan         String {'-'}              them using SPRINTFs defaults.
%
%   caption     String {''}               The caption of the table.
%
%   centering   {'on'} | 'off'            Use \centering.
%
%   colspec     String {'|c| ... |c|'}    The format of each column.
%
%   label       String {''}               The label of the table.
%
%   pos         String {'htp'}            The position of the table.
%
%   Examples:
%     % a small table of the BesselJ-function
%     A = [0 1:5; ...
%          [1:4]' besselj([1:5],[1:4]')];
%     fmt = [{'1_'} frepmat({'7d'},[1 5]); ...
%            frepmat([{'d'} frepmat({'7.4f'},[1 5])],[4])];
%     S = arr2latex(A,fmt,'colspec','|c||r|r|r|r|r|','hline','off')
%
%     % the Gamma-function
%     A = [10:15; gamma(10:15)]';
%     fmt = {'3d' '6.4e$'};
%     T = arr2latex(A,fmt,'colspec','|c|l|', ...
%                         'collabel',{'$x$' '$\Gamma(x)$'}, ...
%                         'hline','off')
%
%   See also SPRINTF.

% S. Engblom 2005-09-22

% checks...
if ~isreal(T) || ndims(T) > 2
  error('Real, two-dimensional matrix expected.');
end
if ndims(fmt) > 2 || ~isa(fmt,'cell') || ...
      ~all(all(cellfun('isclass',fmt,'char') & ...
               cellfun('ndims',fmt) == 2 & ...
               cellfun('size',fmt,1) == 1))
  error('Cell-matrix containing row strings expected.');
end

% expand fmt to match T
rep = size(T)./size(fmt);
if any(rep ~= ceil(rep))
  error('The sizes of the table and the format strings must match.');
end
fmt = frepmat(fmt,rep);

% parse options
opts = struct('collabel','', ...
              'rowlabel','', ...
              'breakrow',74, ...
              'hline','on', ...
              'posinf','$\infty$', ...
              'neginf','$-\infty', ...
              'nan','  -', ...
              'caption','', ...
              'centering','on', ...
              'colspec',['|' frepmat('c|',[1 size(T,2)])], ...
              'label','', ...
              'pos','htp');
[opts,got] = parseopts(opts,varargin);
if strcmpi(opts.hline,'on')
  opts.hline = '\t\\hline\n';
else
  opts.hline = '';
end
if ~isempty(opts.label)
  opts.label = ['\\label{' opts.label '}\n'];
end
if strcmpi(opts.centering,'on')
  opts.centering = '\\centering\n';
else
  opts.centering = '';
end
% rowlabels? add one column!
if ~isempty(opts.rowlabel) && ~got.colspec
  opts.colspec = ['|c|' opts.colspec];
end
if ~isempty(opts.caption)
  opts.caption = ['\\caption{' l_string(opts.caption) '}\n'];
end

% the beginning...
S = ['\\begin{table}[' opts.pos ']\n' ...
     opts.centering ...
     '\\begin{tabular}{' opts.colspec '}\n' ...
     '\t\\hline\n'];

% column labels
if ~isempty(opts.collabel)
  s = '\t';
  for j = 1:size(T,2)+(~isempty(opts.rowlabel))
    s = [s l_string(opts.collabel{j}) ' & '];
  end
  s(end-1:end) = '';
  s = l_break(s,opts.breakrow);
  s = [s '\t\\\\\n\t\\hline\n'];
  S = [S s];
end

% ...the body...
for i = 1:size(T,1)
  s = '\t';
  if ~isempty(opts.rowlabel)
    s = [s l_string(opts.rowlabel{i}) ' & '];
  end
  for j = 1:size(T,2)
    ss = l_format(fmt{i,j},T(i,j));

    % signs for infinity and NaN's, if any -- try to keep the tabs 
    % of the table whenever possible
    if ~isempty(opts.posinf) && T(i,j) == Inf
      ss = [frepmat(' ',[1 max(size(ss,2)-size(opts.posinf,2),0)]) ...
            l_string(opts.posinf)];
    elseif ~isempty(opts.neginf) && T(i,j) == -Inf
      ss = [frepmat(' ',[1 max(size(ss,2)-size(opts.neginf,2),0)]) ...
            l_string(opts.neginf)];
    elseif ~isempty(opts.nan) && isnan(T(i,j))
      ss = [frepmat(' ',[1 max(size(ss,2)-size(opts.nan,2),0)]) ...
            l_string(opts.nan)];
    end
    s = [s ss ' & '];
  end
  s(end-1:end) = '';

  % possibly insert breaks
  s = l_break(s,opts.breakrow);

  % new row
  s = [s '\t\\\\\n' opts.hline];
  S = [S s];
end

% always produce a 'bounding box'
if isempty(opts.hline), S = [S '\t\\hline\n']; end

% ...and the end!
S = [S ...
     '\\end{tabular}\n' ...
     opts.caption ...
     opts.label ...
     '\\end{table}'];

% let sprintf do the rest
S = sprintf(S);

%--------------------------------------------------------------------------
function s = l_string(t)
%L_STRING Formatting one string.
%   S = L_FORMAT(T) formats the string T. Care is taken to insert a
%   double '\\' whenever a single slash is encountered.

s = t(sort([find(t == '\') 1:size(t,2)]));

%--------------------------------------------------------------------------
function s = l_break(t,n)
%L_BREAK Break rows.
%   S = L_BREAK(T,N) tries to insert new lines in the string T at
%   every Nth character position. The new lines are always inserted
%   after a '&'-character and a double tab is used so as to keep the
%   text readable.

ii = find(t == '&');
jj = find(diff(mod(ii,n)) < 0);
if ~isempty(jj)
  jj = ii(jj);
  s = t(1:jj(1));
  for k = [jj+1; jj(2:end) size(t,2)]
    s = [s '\n\t\t' t(k(1):k(2))];
  end
else
  s = t;
end

%--------------------------------------------------------------------------
function s = l_format(f,t)
%L_FORMAT Formatting one scalar.
%   S = L_FORMAT(F,T) formats the scalar T using the flag F and
%   returns the result in the string S.

if f(end) == '_'
  s = frepmat(' ',[1 str2num(f(1:end-1))]);
elseif f(end) == '$'
  s = sprintf(['%' f(1:end-1)],t);
  i = find(s == 'e',1,'first');
  if ~isempty(i)
    s1 = s(1:i-1); s2 = s(i+1:end);
    % masks for 'e-08' ('-8'), 'e+01' ('1'), 'e+00' ('0')
    s2 = num2str(str2num(s2));
    s = ['$' s1 ' \\cdot 10^{' s2 '}$'];
  end
else
  s = sprintf(['%' f],t);
end
%--------------------------------------------------------------------------
