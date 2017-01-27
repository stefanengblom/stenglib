%STARTUP Add paths to FAST stuff.

% S. Engblom 2005-03-22

s = pwd;
addpath([s '/source']);
if exist('test','dir')
  addpath([s '/test']);
end
