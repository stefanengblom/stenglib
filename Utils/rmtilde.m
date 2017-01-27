function nrem = rmtilde(wd)
%RMTILDE Remove files ending with a tilde ('~').
%   NREM = RMTILDE Recursively removes all files with names ending
%   with a tilde, starting from the current directory. The total
%   number of removed files is returned.
%
%   NREM = RMTILDE(WD) does the same thing, but starts from the
%   directory WD instead.
%
%   In any of the syntaxes, filenames or directories starting with a
%   period are left unaffected.
%
%   Cautionary: the function is recursive and may run slowly in
%   certain situations.
%
%   See also DELETE, RECYCLE, RMDIR.

% S. Engblom 2005-06-16

% remember callers working directory
cwd = pwd;
nrem = 0;

% input
if nargin == 0, wd = cwd; end

% change working directory
cd(wd);

% fetch all files/directories
d = dir;

% loop over all of them
for i = 1:size(d,1)
  % avoid all names starting with a period
  if d(i).name(1) ~= '.'
    if d(i).isdir
      % directory: recursive call
      fprintf(1,'Searching directory %s...\n',d(i).name);
      nrem = nrem+rmtilde(d(i).name);
    elseif d(i).name(end) == '~'
      % file ending with '~': remove it
      delete(d(i).name);
      nrem = nrem+1;
      fprintf(1,'Removed file: %s.\n',d(i).name);
    end
  end
end

% switch back to callers directory
cd(cwd);
