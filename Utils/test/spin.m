function spin(ix)
%SPIN Tests for UTILS.

% S. Engblom 2012-03-31 (Revision)
% S. Engblom 2007-01-29 (Revision)
% S. Engblom 2004-10-29 (Revision)
% S. Engblom 2003-03-27

ftests = {@l_spin1 @l_spin2 @l_spin3 @l_spin4 @l_spin5 @l_spin6 @l_spin7};
stests = {'Test #1 (ndop)' 'Test #2 (ndop)' 'Test #3 (ndop)' ...
          'Test #4 (ndop)' 'Test #5 (ndop)' 'Test #6 (ndop)' ...
          'Test #7 (matmerge)'};
if nargin == 0, ix = 1:size(ftests,2); end

runtest('SPIN (utils)',ftests(ix),stests(ix));

%--------------------------------------------------------------------------
function ok = l_spin1
%L_SPIN1 Test of illegal input.

ok = 1;

M = rand(4,4);
try
  S = ndop(M,1,10);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('ndop:e2',li);
end

try
  S = ndop(M,[1 10],2);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('ndop:e1',li);
end

try
  S = ndop(M,2,[1 4]);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('ndop:e1',li);
end

try
  S = ndop(rand(4,4,4),[2 3],[2 3]);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('ndop:e2',li);
end

try
  S = ndop([1 2],[1 1],[10 10],{3 4},[1 3 4]);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('ndop:e3',li);
end

try
  S = ndop(1,1,10,[],1:20);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('ndop:e3',li);
end

try
  [S,T] = ndop([1 2],[1 1],[10 10]);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('ndop:e4',li);
end

% not caught by ndop
try
  S = ndop([1 2],[1 1],[10 10],[3 4]);
  ok = 0;
catch
  ok = ok;
end

%--------------------------------------------------------------------------
function ok = l_spin2
%L_SPIN2 Special cases when target point is on one boundary.

ok = 1;

M = rand(4,4);
S = ndop(M,[1 3],[40 30]);
S = ndop(M,[2 1],[40 30]);
S = ndop(M,[1 1],[40 30]);
S = ndop(M,[4 1],[40 30]);
S = ndop(M,[4 4],[40 30]);
S = ndop(M,[1 4],[40 30]);

% same in 3-D
M = rand(5,5,5);
S = ndop(M,[1 3 4],[8 7 6]);
S = ndop(M,[2 1 5],[8 7 6]);
S = ndop(M,[1 1 1],[8 7 6]);
S = ndop(M,[5 1 5],[8 7 6]);
S = ndop(M,[5 5 1],[8 7 6]);
S = ndop(M,[1 5 3],[8 7 6]);

% from gnats
S = full(ndop([1; 0],[2 1],[3 3]));
ok = ok && all(size(S) == [9 9]) && norm(S-diag(diag(S,-1),-1),'fro') == 0;

% similar case
S = full(ndop([1 0; 0 0],[2 2],[3 3]));
ok = ok && all(size(S) == [9 9]) && norm(S-diag(diag(S,-4),-4),'fro') == 0;

%--------------------------------------------------------------------------
function ok = l_spin3
%L_SPIN3 Test of correct boundary handling.

ok = 1;

sz = [20 33];
M = [0 1/3 0; 1 -4 2; 0 10 0];
S1 = ndop(M,[2 2],sz);
M = [0 1/3 0 0; 1 -4 2 0; 0 10 0 0];
S2 = ndop(M,[2 2],sz);
M = [0 0 1/3 0; 0 1 -4 2; 0 0 10 0];
S3 = ndop(M,[2 3],sz);
M = [0 0 0 0; 0 0 1/3 0; 0 1 -4 2; 0 0 10 0; 0 0 0 0];
S4 = ndop(M,[3 3],sz);

ok = ok && all([norm(S1(:)-S2(:)) ...
                norm(S2(:)-S3(:)) ...
                norm(S3(:)-S4(:))] < 1e-16);

%--------------------------------------------------------------------------
function ok = l_spin4
%L_SPIN4 "Collapsed 1-D" case.

ok = 1;

sz = [1 20];
M1 = rand(1,4);
S1 = ndop(M1,[1 3],sz);
M2 = [zeros(1,6); [0 M1 0]; zeros(1,6)];
S2 = ndop(M2,[2 4],sz);
S3 = ndop(M1',[3 1],sz(end:-1:1));
S4 = ndop(M2',[4 2]',sz(end:-1:1));

ok = ok && all([norm(S1(:)-S2(:)) ...
                norm(S2(:)-S3(:)) ...
                norm(S3(:)-S4(:))] < 1e-16);

%--------------------------------------------------------------------------
function ok = l_spin5
%L_SPIN5 Empty molecules.

ok = 1;

S = ndop(0,1,3);
ok = ok && all(size(S) == [3 3]) && nnz(S) == 0; 
S = ndop([0; 0],1,4);
ok = ok && all(size(S) == [4 4]) && nnz(S) == 0;
S = ndop([0 0; 0 0],[1 1],[4 4]);
ok = ok && all(size(S) == [16 16]) && nnz(S) == 0;

%--------------------------------------------------------------------------
function ok = l_spin6
%L_SPIN6 Complex cases.

ok = 1;

% from help
N = 10; h = 1/(N-1);
M = [0 1 0; 1 -4 1; 0 1 0]/h^2;
S = ndop(M,[2 2],[N N]);
S = (1+2i)*S;

S1 = ndop((1+2i)*M,[2 2],[N N]);
S2 = ndop(M,[2 2],[N N],frepmat(1+2i,[N N]),[]);
S3 = ndop(M,[2 2],[N N],[],frepmat(1+2i,[N N]));

ok = ok && norm(S-S1,'fro') == 0 && ...
     norm(S-S2,'fro') == 0 && norm(S-S3,'fro') == 0;

%--------------------------------------------------------------------------
function ok = l_spin7
%L_SPIN 7 Test of MATMERGE.

ok = 1;

% ensure that file exists
arr = nan;
s = nan;
t = nan;
save test/spin7c arr s t;
try
  matmerge('test/spin7c','test/spin7a','test/spin7b');
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp(li,'matmerge:e1');
end

% remove it
delete('test/spin7c.mat');
try
  matmerge('test/spin7c','test/spin7a');
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp(li,'matmerge:e2');
end

% standard behavior
matmerge('test/spin7c','test/spin7a','test/spin7b');

load test/spin7c
ok = ok && all(size(arr) == [2 7]) && ...
     all(size(s) == [2 1]) && ...
     all(size(t) == [1 10 2]);

matmerge('test/spin7c','test/spin7a');
load test/spin7c
ok = ok && all(size(arr) == [2 10]) && ...
     all(size(s) == [3 1]) && ...
     all(size(t) == [1 10 3]);

%--------------------------------------------------------------------------
