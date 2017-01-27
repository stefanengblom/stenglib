%STARTUP Add paths to SCICOMP stuff.

% S. Engblom 2008-06-17

s = pwd;
if exist('test','dir')
  addpath([s '/test']);
end
