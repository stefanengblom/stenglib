%STARTUP Add paths to TENSOR stuff.

% S. Engblom 2005-04-10

s = pwd;
if exist('test','dir')
  addpath([s '/test']);
end
addpath([s '/source']);
