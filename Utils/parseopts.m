function  [opts,got] = parseopts(optdef,optin)
%PARSEOPTS Parse options.
%   [OPTS,GOT] = PARSEOPTS(OPTDEF,OPTIN) merges the options OPTIN with
%   the default options OPTDEF and produces the options structure
%   OPTS. OPTDEF and OPTIN can either be structures or cell-vectors
%   with property/value pairs.
%
%   Fields in OPTDEF that are structures are recursively parsed, see
%   the example below.
%
%   The optional output GOT is a structure containing the same fields
%   as OPTS. The values of each field are boolean variables indicating
%   whether the corresponging field was received in OPTIN or not.
%
%   It is considered an error for OPTIN to contain options not in
%   OPTDEF.
%
%   Example:
%     f1 = struct('delta',0.5,'rho',0.8);
%     f2 = struct('delta',0.2,'rho',0.1);
%     od = struct('nlin','yes','theta',0.1,'relax','on', ...
%                 'filter1',f1,'filter2',f2);
%     oi = {'theta' 0.25 'relax' 'off' 'filter1' {'rho' 0.2}};
%     [opts,got] = parseopts(od,oi)
%
%   See also ASSIGNOPTS.

% S. Engblom 2005-07-21

% extract fields and values
if isstruct(optdef)
  if any(size(optdef) ~= 1)
    error('Expecting a scalar structure.');
  end
  dfield = fieldnames(optdef);
  dval = struct2cell(optdef);
elseif iscell(optdef)
  dfield = reshape(optdef(1:2:end),[],1);
  dval = reshape(optdef(2:2:end),[],1);
  if size(dfield,1) ~= size(dval,1)
    error('Cell-vector must contain property/value pairs.');
  end
else
  error('Default options must either be a struct or a cell-vector.');
end

if isstruct(optin)
  if any(size(optin) ~= 1)
    error('Expecting a scalar structure.');
  end
  ifield = fieldnames(optin);
  ival = struct2cell(optin);
elseif iscell(optin)
  ifield = reshape(optin(1:2:end),[],1);
  ival = reshape(optin(2:2:end),[],1);
  if size(ifield,1) ~= size(ival,1)
    error('Cell-vector must contain property/value pairs.');
  end
else
  error('Input options must either be a struct or a cell-vector.');
end

% check for illegal fields
d = fsetop('setdiff',ifield,dfield);
if ~isempty(d)
  error(sprintf('Unsupported option "%s".',d{1}));
end

% merge and create structure of options
[sfield,ii,id] = fsetop('union',ifield,dfield);
opts = cell2struct([ival(ii); dval(id)],sfield(:));
got = cell2struct([frepmat({true},size(ii,2)); ...
                   frepmat({false},size(id,2))],sfield(:));

% recursive call for fields containing structures
is = find(cellfun('isclass',dval,'struct'))';
for i = is
  f = dfield{i};
  if got.(f)
    [opts.(f),got.(f)] = parseopts(dval{i},opts.(f));
  else
    [opts.(f),got.(f)] = parseopts(dval{i},{});
  end
end
