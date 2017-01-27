function assignopts(opts)
%ASSIGNOPTS Assign options.
%   ASSIGNOPTS(OPTS) assigns the options in OPTS. This is useful
%   whenever you prefer to use 'name' rather than 'OPTS.name'
%   throughout the code. OPTS can either be a structure or a
%   cell-vector with property/value pairs.
%
%   Example:
%     % caller's site:
%     opts.order = 2;
%
%     % inside fun(opts, ...)
%     def.rho = 0.01;
%     def.method = 'lagrange';
%     def.order = 4;
%     opts = parseopts(def,opts);
%     assignopts(opts);
%
%     % use the variables directly
%     order % note: 2
%     method
%
%   See also PARSEOPTS, ASSIGNIN.

% S. Engblom 2010-08-25

if isstruct(opts)
  if any(size(opts) ~= 1)
    error('Expecting a scalar structure.');
  end
  field = fieldnames(opts);
  val = struct2cell(opts);
elseif iscell(opts)
  field = reshape(opts(1:2:end),[],1);
  val = reshape(opts(2:2:end),[],1);
  if size(field,1) ~= size(val,1)
    error('Cell-vector must contain property/value pairs.');
  end
else
  error('Options must either be a struct or a cell-vector.');
end

% immediate
for i = 1:size(field,1)
  assignin('caller',field{i},val{i});
end
