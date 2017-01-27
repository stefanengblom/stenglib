function ok = spin(ix)
%SPIN Tests for FAST.

% S. Engblom 2015-01-19 (removed empty rep in repmat)
% S. Engblom 2013-12-03 (silent syntax)
% S. Engblom 2010-02-09 (extended syntax for 'union')
% S. Engblom 2010-01-13 ('nosort'-tests modyfied)
% S. Engblom 2005-08-03 (clenshaw added)
% S. Engblom 2005-06-22 (fsetop added)
% S. Engblom 2005-04-27 (merging frepmat/fsparse)
% S. Engblom 2004-10-22

ftests = {@l_spin1 @l_spin2 @l_spin3 @l_spin4 ...
          @l_spin5 @l_spin6 @l_spin7 @l_spin8 ...
          @l_spin9 @l_spin10 @l_spin11 @l_spin12 ...
          @l_spin13 @l_spin14 @l_spin15 @l_spin16 @l_spin17 ...
          @l_spin18 @l_spin19 @l_spin20 @l_spin21 @l_spin22};
stests = {'Test #1 (frepmat)' 'Test #2 (frepmat)' 'Test #3 (frepmat)' ...
          'Test #4 (frepmat)' 'Test #5 (frepmat)' 'Test #6 (frepmat)' ...
          'Test #7 (frepmat)' 'Test #8 (fsparse)' 'Test #9 (fsparse)' ...
          'Test #10 (fsparse)' 'Test #11 (fsparse)' 'Test #12 (fsparse)' ...
          'Test #13 (fsetop)' 'Test #14 (fsetop)' ...
          'Test #15 (fsetop)' 'Test #16 (fsetop)' 'Test #17 (fsetop)' ...
          'Test #18 (clenshaw)' 'Test #19 (clenshaw)' ...
          'Test #20 (clenshaw)' 'Test #21 (clenshaw)' ...
          'Test #22 (clenshaw)'};
if nargin == 0, ix = 1:size(ftests,2); end
if nargout > 0
  % for the case that spin is called from scrub
  ok = runtest('',ftests(ix),stests(ix));
else
  runtest('SPIN (fast)',ftests(ix),stests(ix));
end

%--------------------------------------------------------------------------
function ok = l_spin1
%L_SPIN1 (FREPMAT) Test of illegal input.

ok = 1;

try
  a = mexfrepmat(1,2,3);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('frepmat:e1',li);
end

try
  a = mexfrepmat({1,2,3},[1 2 3]);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('frepmat:e2',li);
end

try
  a = mexfrepmat(struct('foo','bar'),[2 3]);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('frepmat:e2',li);
end

try
  a = mexfrepmat(@sin,[2 3]);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('frepmat:e2',li);
end

try
  a = mexfrepmat(inline('sin(x)'),[2 3]);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('frepmat:e2',li);
end

try
  a = frepmat([1 2 3],[1 2 3]+i);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('frepmat:e3',li);
end

try
  a = frepmat([1 2 3],sparse([1 2 3]));
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('frepmat:e3',li);
end

try
  a = frepmat([1 2 3],[1 2.3 3]);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('frepmat:e4',li);
end

try
  a = frepmat(sparse([1 0 3; 0 0 1]),[2 2 3]);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('frepmat:e5',li);
end

try
  a = frepmat([1 0 3; 0 0 1],sparse([2 1]));
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('frepmat:e3',li);
end

% one positive case...
try
  a = frepmat([1 2 3],[1 2 3]);
  ok = ok & all(size(a) == [1 6 3]);
catch
  ok = 0;
end

% from GNATS
try
  frepmat(1,[1 NaN]);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('frepmat:e4',li);
end

try
  frepmat({1},[1 NaN]);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('frepmat:e4',li);
end

%--------------------------------------------------------------------------
function ok = l_spin2
%L_SPIN2 (FREPMAT) Various small cases.

ok = 1;

for bb = {1:5 [1:7]' reshape(1:6,2,3) reshape(1:5,1,1,5) ...
          reshape(1:24,[2 1 3 1 1 4]) reshape(1:16,2,4,2)}
  b = bb{1};
  a1 = frepmat(b,2);
  a2 = repmat(b,2,1);
  ok = ok & isequal(a1,a2);

  a1 = frepmat(b,[1 2]);
  a2 = repmat(b,[1 2]);
  ok = ok & isequal(a1,a2);

  a1 = frepmat(b,[2 1]);
  a2 = repmat(b,[2 1]);
  ok = ok & isequal(a1,a2);

  a1 = frepmat(b,[2 3]);
  a2 = repmat(b,[2 3]);
  ok = ok & isequal(a1,a2);
end

b = reshape(1:24,[2 1 3 1 1 4]);
a1 = frepmat(b,[1 2 3]);
a2 = repmat(b,[1 2 3]);
ok = ok & isequal(a1,a2);

a1 = frepmat(b,[2 1 1]);
a2 = repmat(b,[2 1 1]);
ok = ok & isequal(a1,a2);

a1 = frepmat(b,[1 1 1 4]);
a2 = repmat(b,[1 1 1 4]);
ok = ok & isequal(a1,a2);

a1 = frepmat(b,[1 4 1 4]);
a2 = repmat(b,[1 4 1 4]);
ok = ok & isequal(a1,a2);

a1 = frepmat(b,[1 3 1 1]);
a2 = repmat(b,[1 3 1 1]);
ok = ok & isequal(a1,a2);

b = reshape(1:16,2,4,2);
a1 = frepmat(b,[1 4 1 4]);
a2 = repmat(b,[1 4 1 4]);
ok = ok & isequal(a1,a2);

a1 = frepmat(b,[2 3 1 1 4]);
a2 = repmat(b,[2 3 1 1 4]);
ok = ok & isequal(a1,a2);

%--------------------------------------------------------------------------
function ok = l_spin3
%L_SPIN3 (FREPMAT) Empty cases.

ok = 1;

for bb = {ones(2,3,4) ones(2,0,4) ones(0,3,4) ones(2,3,0)}
  b = bb{1};
  a1 = frepmat(b,[]);
  a2 = repmat(b,1);
  ok = ok & isequal(a1,a2);

  a1 = frepmat(b,ones(2,0));
  a2 = b; % produces a warning: repmat(b,ones(2,0));
  ok = ok & isequal(a1,a2);

  a1 = frepmat(b,[1 0]);
  a2 = repmat(b,[1 0]);
  ok = ok & isequal(a1,a2);

  a1 = frepmat(b,[0 1]);
  a2 = repmat(b,[0 1]);
  ok = ok & isequal(a1,a2);

  a1 = frepmat(b,[2 0 1]);
  a2 = repmat(b,[2 0 1]);
  ok = ok & isequal(a1,a2);

  a1 = frepmat(b,[0 3 1]);
  a2 = repmat(b,[0 3 1]);
  ok = ok & isequal(a1,a2);

  a1 = frepmat(b,[0 0 1]);
  a2 = repmat(b,[0 0 1]);
  ok = ok & isequal(a1,a2);

  a1 = frepmat(b,[0 0 2 1]);
  a2 = repmat(b,[0 0 2 1]);
  ok = ok & isequal(a1,a2);

  a1 = frepmat(b,[0 0 2 3 0]);
  a2 = repmat(b,[0 0 2 3 0]);
  ok = ok & isequal(a1,a2);
end

%--------------------------------------------------------------------------
function ok = l_spin4
%L_SPIN4 (FREPMAT) Non-doubles.

ok = 1;

for ttype = {'single' 'uint8' 'uint16' 'uint32' 'uint64' ...
             'int8' 'int16' 'int32' 'int64'}
  type = ttype{1};
  b = feval(type,1:5);
  a1 = frepmat(b,[1 2]);
  a2 = frepmat(b,[1 2]);
  ok = ok & isequal(a1,a2);

  a1 = frepmat(b,[3 1]);
  a2 = frepmat(b,[3 1]);
  ok = ok & isequal(a1,a2);

  a1 = frepmat(b,[2 1 3]);
  a2 = frepmat(b,[2 1 3]);
  ok = ok & isequal(a1,a2);

  b = feval(type,reshape(1:8,1,1,8));
  a1 = frepmat(b,[1 2 3]);
  a2 = frepmat(b,[1 2 3]);
  ok = ok & isequal(a1,a2);
end

% characters and logicals
for bb = {'goo, went, gone!' 'MatLab' == 'MATLAB'}
  b = bb{1};
  a1 = frepmat(b,[1 2]);
  a2 = frepmat(b,[1 2]);
  ok = ok & isequal(a1,a2);

  a1 = frepmat(b,[3 1]);
  a2 = frepmat(b,[3 1]);
  ok = ok & isequal(a1,a2);

  a1 = frepmat(b,[2 1 3]);
  a2 = frepmat(b,[2 1 3]);
  ok = ok & isequal(a1,a2);
end

b = reshape('bye, bye',1,1,8);
a1 = frepmat(b,[1 2 3]);
a2 = frepmat(b,[1 2 3]);
ok = ok & isequal(a1,a2);

b = reshape('Bye, Bye' < 'BYE, BYE',1,1,8);;
a1 = frepmat(b,[1 2 3]);
a2 = frepmat(b,[1 2 3]);
ok = ok & isequal(a1,a2);

%--------------------------------------------------------------------------
function ok = l_spin5
%L_SPIN5 (FREPMAT) Complex cases.

ok = 1;

c = 3+i;
for ttype = {'double' 'single' 'uint16' 'int32' 'int64'}
  type = ttype{1};
  b = feval(type,c*[1:5]);
  a1 = frepmat(b,[1 2]);
  a2 = frepmat(b,[1 2]);
  ok = ok & isequal(a1,a2);

  a1 = frepmat(b,[3 1]);
  a2 = frepmat(b,[3 1]);
  ok = ok & isequal(a1,a2);

  a1 = frepmat(b,[2 1 3]);
  a2 = frepmat(b,[2 1 3]);
  ok = ok & isequal(a1,a2);

  b = feval(type,c*reshape(1:8,1,1,8));
  a1 = frepmat(b,[1 2 3]);
  a2 = frepmat(b,[1 2 3]);
  ok = ok & isequal(a1,a2);
  
  b = feval(type,c*reshape(1:16,1,1,8,2));
  a1 = frepmat(b,[1 1 3 2]);
  a2 = frepmat(b,[1 1 3 2]);
  ok = ok & isequal(a1,a2);
end

%--------------------------------------------------------------------------
function ok = l_spin6
%L_SPIN6 (FREPMAT) Sparse matrices.

ok = 1;

for c = [1 3+2*i]
  for bb = {sparse([1 0 0; 0 0 3; 2 4 0; 1 0 7]) ...
            sparse([1 2 0; 0 0 0; 3 4 0; 0 0 0]) ...
            sparse([0 0 0; 0 2 0; 0 4 0; 0 7 0])}
    b = c*bb{1};
    a1 = frepmat(b,2);
    a2 = repmat(b,[2 1]);
    ok = ok & isequal(a1,a2) & nzmax(a1) == 2*nnz(b);

    a1 = frepmat(b,[1 3]);
    a2 = repmat(b,[1 3]);
    ok = ok & isequal(a1,a2) & nzmax(a1) == 3*nnz(b);

    a1 = frepmat(b,[4 3]);
    a2 = repmat(b,[4 3]);
    ok = ok & isequal(a1,a2) & nzmax(a1) == 4*3*nnz(b);
  end
end

% empty sparse cases
for bb = {sparse([2 3 0; 0 0 1]) sparse(0,2) sparse(2,0) sparse(2,2)}
  b = bb{1};
  a1 = frepmat(b,[0 2]);
  a2 = repmat(b,[0 2]);
  ok = ok & isequal(a1,a2);

  a1 = frepmat(b,[3 0]);
  a2 = repmat(b,[3 0]);
  ok = ok & isequal(a1,a2);

  a1 = frepmat(b,[0 0]);
  a2 = repmat(b,[0 0]);
  ok = ok & isequal(a1,a2);

  a1 = frepmat(b,[3 2]);
  a2 = repmat(b,[3 2]);
  ok = ok & isequal(a1,a2);
end

% sparse logicals
for bb = {sparse([1 0 0; 0 0 3; 2 4 0; 1 0 7]) ...
          sparse([1 2 0; 0 0 0; 3 4 0; 0 0 0]) ...
          sparse([0 0 0; 0 2 0; 0 4 0; 0 7 0])}
  b = bb{1} > 0;
  a1 = frepmat(b,2);
  a2 = repmat(b,[2 1]);
  ok = ok & isequal(a1,a2) & nzmax(a1) == 2*nnz(b);

  a1 = frepmat(b,[1 3]);
  a2 = repmat(b,[1 3]);
  ok = ok & isequal(a1,a2) & nzmax(a1) == 3*nnz(b);

  a1 = frepmat(b,[4 3]);
  a2 = repmat(b,[4 3]);
  ok = ok & isequal(a1,a2) & nzmax(a1) == 4*3*nnz(b);
end

%--------------------------------------------------------------------------
function ok = l_spin7
%L_SPIN7 (FREPMAT) Cell- and structure-arrays, function-handles and
%objects.

ok = 1;

bb1 = {{} [1:3] [1:4]' pi {2*pi [1:5]}};
bb2 = cell(2,3); bb2{1,2} = 1:4; bb2{1,3} = {pi 2*pi};
bb3 = cell(2,0);
bb4 = cell(0,2);
bb5 = cell(0,0);
bb6 = struct('foo','bar','boo',1:4);
bb7 = {@sin @cos};
bb8 = {inline('sin(x)')};
for bb = {bb1 bb2 bb3 bb4 bb5 bb6 bb7 bb8}
  b = bb{1};
  a1 = frepmat(b,[1 2 3]);
  a2 = repmat(b,[1 2 3]);
  ok = ok & isequal(a1,a2);

  a1 = frepmat(b,[1 2 0]);
  a2 = repmat(b,[1 2 0]);
  ok = ok & isequal(a1,a2);

  a1 = frepmat(b,[0 2 0]);
  a2 = repmat(b,[0 2 0]);
  ok = ok & isequal(a1,a2);

  a1 = frepmat(b,[3 0]);
  a2 = repmat(b,[3 0]);
  ok = ok & isequal(a1,a2);
end
%--------------------------------------------------------------------------
function ok = l_spin8
%L_SPIN8 (FSPARSE) Test of illegal input.

ok = 1;

try
  s = fsparse(sparse(1:4) < 2);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e2',li);
end

try
  s = fsparse(ones(2,3,4));
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e3',li);
end

try
  s = fsparse(1,2,3,4,5,6,7);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e1',li);
end

try
  [s,t,u] = fsparse(1,2,3);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e1',li);
end

try
  s = fsparse(1,2+3i,3);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e4',li);
end

try
  s = fsparse(uint32(1),2,3);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e4',li);
end

try
  s = fsparse(1,2,sparse(3));
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e5',li);
end

try
  s = fsparse(ones(2,2,2),ones(2,2,2),zeros(2,2,2));
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e6',li);
end

try
  s = fsparse(1,2,3,int32([1 2]));
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e7',li);
end

try
  s = fsparse(1,2,3,[1 2],3);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e8',li);
end

try
  s = fsparse(1:3,1:2,3);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e9',li);
end

try
  s = fsparse(1:3,[1:2]',[3 3]);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e9',li);
end

try
  s = fsparse(1:3,[1:2]',[3 3 3]');
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e9',li);
end

try
  s = fsparse(int32(0:3)',1,2);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e10',li);
end

try
  s = fsparse(int32(0:3),int32(1:4),2);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e10',li);
end

try
  s = fsparse(0:3,int32(1:4),2);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e10',li);
end

try
  s = fsparse(int32(0:3),int32(1:4),2);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e10',li);
end

try
  s = fsparse(1:3,[1 1.2 2.3],2);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e10',li);
end

try
  s = fsparse(1:3,[1 NaN 2.3],2);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e10',li);
end

try
  s = fsparse(1:3,1:3,2,[2.2 3]);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e11',li);
end

try
  s = fsparse(1:3,1:3,2,[2 3]);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e12',li);
end

try
  s = fsparse(1:3,1:3,2,[3 3 -4]);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e13',li);
end

try
  s = fsparse(1:3,1:3,2,[3 3 4 5]);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e14',li);
end

try
  s = fsparse(1:10,1:10,1:10,[15 15 9]);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e16',li);
end

try
  s = fsparse(1:10,1:10,1:10,[15 15 9],'nosort');
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e16',li);
end

try
  s = fsparse([1:3]',1:2,1,[15 15 4]);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e16',li);
end

try
  s = fsparse([1:3]',1:2,1,[15 15 4],'nosort');
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e16',li);
end

try
  s = fsparse(1:3,1:3,2,[3 3],'~sort');
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e15',li);
end

try
  s = fsparse(1:3,1:3,2,[3 3],'alfkjfk94esgwrgewrgh3ertherhewrghewrhrhr');
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e15',li);
end

% from GNATS
try
  fsparse(1:3,[1 NaN 2],1);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e10',li);
end

try
  fsparse(1:3,[1 3 2],1,[6 NaN]);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e11',li);
end

try
  fsparse(1:3,[1 3 2],1,[6 4 NaN]);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e13',li);
end

try
  fsparse(1:3,[1 3 2],1,[],'sort',[1 2]);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e17',li);
end

try
  fsparse(1:3,[1 3 2],1,[],'nosort',-2);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsparse:e18',li);
end

% should work, though...
ix = 1:3;
ss = ones(1,3);
s = fsparse(ix,ix,ss);
s = fsparse(ix,ix,ss,[]);
s = fsparse(ix,ix,ss,[4]);
s = fsparse(ix,ix,ss,[4 4]);
s = fsparse(ix,ix,ss,[4 4],'nosort');
s = fsparse(ix,ix,ss,[4 4 4]);
s = fsparse(ix,ix,ss,[4 4 4],'nosort');
s = fsparse(ix,ix,ss,[4 4],'sort',2);

ix = int32(ix);
s = fsparse(ix,ix,ss);
s = fsparse(ix,ix,ss,[]);
s = fsparse(ix,ix,ss,[4]);
s = fsparse(ix,ix,ss,[4 4]);
s = fsparse(ix,ix,ss,[4 4],'nosort');
s = fsparse(ix,ix,ss,[4 4 4]);
s = fsparse(ix,ix,ss,[4 4 4],'nosort');
s = fsparse(ix,ix,ss,[4 4 4],'nosort',4);

s = fsparse(ones(3,4));
s = fsparse(sparse(ones(3,4)));

%--------------------------------------------------------------------------
function ok = l_spin9
%L_SPIN9 (FSPARSE) Small tests.

ok = 1;

for c = [1 1+1i]
  for typec = {@int32 @double}
    for flagc = {'sort' 'nosort'}
      type = typec{1};
      flag = flagc{1};

      ii = 1:3;
      jj = [1 1 1];
      ss = c*[1:3];
      S1 = fsparse(feval(type,ii),feval(type,jj),ss,[6 6],flag);
      S2 = sparse(ii,jj,ss,6,6);
      ok = ok & isequal(S1,S2) & nnz(S1) == nzmax(S1);

      S1 = fsparse(feval(type,ii),feval(type,jj),ss,[12 10],flag);
      S2 = sparse(ii,jj,ss,12,10);
      ok = ok & isequal(S1,S2) & nnz(S1) == nzmax(S1);

      ii = 4:6;
      jj = [6 6 6];
      ss = c*[1:3];
      S1 = fsparse(feval(type,ii),feval(type,jj),ss,[],flag);
      S2 = sparse(ii,jj,ss);
      ok = ok & isequal(S1,S2) & nnz(S1) == nzmax(S1);

      S1 = fsparse(feval(type,ii),feval(type,jj),ss,[12 10],flag);
      S2 = sparse(ii,jj,ss,12,10);
      ok = ok & isequal(S1,S2) & nnz(S1) == nzmax(S1);

      iic = {1:3 [1 1 1] 1:3 [1 1 1]};
      jjc = {1:3 1:3 [1 1 1] [1 1 1]};
      ss = c*[1:3];
      for i = 1:4
        ii = iic{i};
        jj = jjc{i};
        S1 = fsparse(feval(type,ii),feval(type,jj),ss,[3 3],flag);
        S2 = sparse(ii,jj,ss,3,3);
        ok = ok & isequal(S1,S2) & nnz(S1) == nzmax(S1);

        S1 = fsparse(feval(type,ii),feval(type,jj),ss,[8 8],flag);
        S2 = sparse(ii,jj,ss,8,8);
        ok = ok & isequal(S1,S2) & nnz(S1) == nzmax(S1);
      end

      iic = {[1:3 1 1 1] [1:3 1:3 1:3] [1:3 1:3 1:3 1 1 1 2 2 2 3 3 3]};
      jjc = {[1 1 1 1:3] [1 1 1 2 2 2 3 3 3] [1 1 1 2 2 2 3 3 3 1:3 1:3 1:3]};
      ssc = {c*[1:6] c*[1:9] c*[1:18]};
      for i = 1:3
        ii = iic{i};
        jj = jjc{i};
        ss = ssc{i};
        S1 = fsparse(feval(type,ii),feval(type,jj),ss,[],flag);
        S2 = sparse(ii,jj,ss);
        ok = ok & isequal(S1,S2) & nnz(S1) == nzmax(S1);

        S1 = fsparse(feval(type,ii),feval(type,jj),ss,[8 8],flag);
        S2 = sparse(ii,jj,ss,8,8);
        ok = ok & isequal(S1,S2) & nnz(S1) == nzmax(S1);
      end
    end
  end
end

%--------------------------------------------------------------------------
function ok = l_spin10
%L_SPIN10 (FSPARSE) Test of 'assembly' syntax.

ok = 1;

% the cases are 2x2 3x5 6x4
siz = [8 9];
for ic = 1:3
  switch ic
   case 1
    iic = {[2 8; 8 4] [1 8; 8 8] [2; 8] [8; 8]};
    jjc = {[9 3; 1 2] [1 1; 9 8] [1 9] [9 5]};
    ssc = {[1 2; 3 4] [1; 2] [1 3] 1};
   case 2
    iic = {[2 1 2 3 5; 1 1 8 3 4; 1 8 8 1 2] ...
           [2 3 2 2 1; 4 5 4 6 5; 2 2 8 3 8] ...
           [2; 8; 8] [1; 7; 8]};
    jjc = {[9 3 2 3 4; 1 2 7 8 7; 3 4 2 1 9] ...
           [1 1 2 3 2; 1 1 1 1 1; 7 8 8 9 8] ...
           [1 9 4 4 3] [1 2 8 9 5]};
    ssc = {[1 2 3 4 5; 6 7 8 9 1; 2 3 4 5 6] ...
           [1; 2; 3] [1 3 5 7 9] 1};
   case 3
    iic = {[2 1 2 3; 1 1 8 3; 1 8 8 2; 2 3 2 3; 3 4 5 6; 4 8 7 6] ...
           [2 3 2 2; 4 5 4 6; 2 2 8 3; 8 8 8 8; 7 7 6 6; 4 5 6 4] ...
           [1; 2; 6 ; 7; 8; 8;] [7; 7; 7; 1; 2; 8]};
    jjc = {[9 3 2 3; 1 2 7 8; 3 4 2 9; 3 3 9 9; 9 1 2 3; 2 1 2 3] ...
           [1 1 2 3; 1 1 1 1; 7 8 9 8; 5 6 6 4; 1 2 1 1; 3 3 3 2] ...
           [1 9 5 4] [5 8 9 5]};
    ssc = {[1 2 3 4; 6 7 8 9; 2 3 4 5; 6 7 8 9; 1 2 3 4; 5 6 7 8] ...
           [1; 2; 3; 4; 5; 6] [1 3 5 7] 1};
  end

  for i1 = 1:length(iic)
    for i2 = 1:length(jjc)
      for i3 = 1:length(ssc)
        for c = [1 1+2i]
          S = l_sparse(iic{i1},jjc{i2},c*ssc{i3});
          for s = 0:2
            S1 = fsparse(iic{i1},jjc{i2},c*ssc{i3},siz(1:s));
            S2 = fsparse(iic{i1},jjc{i2},c*ssc{i3},siz(1:s),'nosort');
            S2 = (S2')'; % 100113: newer version of operator minus does not work
                         % properly for unsorted matrix
            ok = ok & isequal(S1,S) & norm(S2-S,inf) < 1e-14 & ...
                 nnz(S1) == nzmax(S1) & nnz(S2) == nzmax(S2);
            if ~ok, [ic i1 i2 i3 c s], return; end
            % nzmax-cases
            S3 = fsparse(iic{i1},jjc{i2},c*ssc{i3},[siz nzmax(S1)+s]);
            S4 = fsparse(iic{i1},jjc{i2},c*ssc{i3},[siz nzmax(S2)+s]);
            ok = ok & isequal(S1,S3) & isequal(S1,S3) & ...
            nzmax(S3) == nzmax(S1)+s & nzmax(S4) == nzmax(S2)+s;
            if ~ok, [ic i1 i2 i3 c s], return; end
          end
        end
      end
    end
  end
end

%--------------------------------------------------------------------------
function S = l_sparse(i,j,s);
%L_SPARSE Simulator of the 'assembly' syntax in FSPARSE.

% no error-checking whatsoever
isiz = size(i);
jsiz = size(j);
ssiz = size(s);

j = repmat(j,[isiz(1)/jsiz(1) 1]);
i = repmat(i,[1 jsiz(2)/isiz(2)]);
s = repmat(s,[isiz(1) jsiz(2)]./ssiz);

S = sparse(i(:),j(:),s(:));

%--------------------------------------------------------------------------
function ok = l_spin11
%L_SPIN11 (FSPARSE) Empty input.

ok = 1;

siz = [17 7 4];
for s = 0:3
  S{1} = fsparse(ones(1,0),ones(1,0),ones(1,0),siz(1:s));
  S{2} = fsparse(ones(0,2),ones(0,2),ones(0,2),siz(1:s));

  S{3} = fsparse(ones(1,1),ones(1,0),ones(1,1),siz(1:s),'nosort');
  S{4} = fsparse(ones(0,2),ones(1,2),ones(0,2),siz(1:s),'sort');

  S{5} = fsparse(ones(0,0),ones(0,0),1,siz(1:s));
  S{6} = fsparse(ones(0,0),ones(0,0),ones(0,1),siz(1:s),'nosort');
  S{7} = fsparse(ones(0,0),ones(0,0),ones(1,0),siz(1:s),'sort');
  
  for i = 1:7
    Ssiz = [size(S{i}) nzmax(S{i})];
    ok = ok & all(siz(1:s) == Ssiz(1:s));
  end
end

% complex cases
S{1} = fsparse(ones(1,0),ones(1,0),1+1i);
S{2} = fsparse(ones(0,0),ones(1,0),1+1i);
S{3} = fsparse(ones(0,1),ones(0,0),1+1i);
S{4} = fsparse(ones(0,0),ones(0,0),1+1i);
S{5} = fsparse(ones(1,1),ones(1,0),1+1i);
S{6} = fsparse(ones(0,1),ones(1,1),1+1i);
for i = 1:6
  ok = ok & ~isreal(S{i});
end

S = fsparse([],[],[]);
ok = ok & all(size(S) == 0);

% from GNATS
S = fsparse(ones(1,0),ones(1,0),ones(1,0),[1 0]);
ok = ok && all(size(S) == [1 0]);
S = fsparse(ones(1,0),ones(1,0),ones(1,0),[0 0]);
ok = ok && all(size(S) == 0);
S = fsparse(ones(1,0),ones(1,0),ones(1,0),[1 1]);
ok = ok && all(size(S) == [1 1]);

% full2sparse/sparse2sparse
S = fsparse([]);
ok = ok & all(size(S) == 0);
S = fsparse(S);
ok = ok & all(size(S) == 0);
S = fsparse(ones(2,0));
ok = ok & all(size(S) == [2 0]);
S = fsparse(S);
ok = ok & all(size(S) == [2 0]);
S = fsparse(ones(0,3));
ok = ok & all(size(S) == [0 3]);
S = fsparse(S);
ok = ok & all(size(S) == [0 3]);

%--------------------------------------------------------------------------
function ok = l_spin12
%L_SPIN12 (FSPARSE) Test of 'squeezing' of zero elements, test of
%full2sparse/sparse2sparse at the same time.

ok = 1;

for c = [1 1+2i]
  % full case
  i = [1:4]'; j = 1:5;

  ss = c*reshape(1:20,4,5);
  ss(2,[1 3]) = 0;
  S = fsparse(i,j,ss);
  ok = ok & nnz(S) == sum(ss(:) ~= 0) & isequal(full(S),ss);
  S1 = fsparse(S);
  S2 = fsparse(ss);
  ok = ok & isequal(S,S1) & isequal(S1,S2) & nnz(S1) == nzmax(S1);

  ss(3:4,3) = 0;
  S = fsparse(i,j,ss);
  ok = ok & nnz(S) == sum(ss(:) ~= 0) & isequal(full(S),ss);
  S1 = fsparse(S);
  S2 = fsparse(ss);
  ok = ok & isequal(S,S1) & isequal(S1,S2) & nnz(S1) == nzmax(S1);

  ss(4,1) = 0; ss(4,5) = 0;
  S = fsparse(i,j,ss);
  ok = ok & nnz(S) == sum(ss(:) ~= 0) & isequal(full(S),ss);
  S1 = fsparse(S);
  S2 = fsparse(ss);
  ok = ok & isequal(S,S1) & isequal(S1,S2) & nnz(S1) == nzmax(S1);

  ss(1,4:5) = 0; ss(2,4) = 0;
  S = fsparse(i,j,ss);
  ok = ok & nnz(S) == sum(ss(:) ~= 0) & isequal(full(S),ss);
  S1 = fsparse(S);
  S2 = fsparse(ss);
  ok = ok & isequal(S,S1) & isequal(S1,S2) & nnz(S1) == nzmax(S1);

  ss(4,4) = 0; ss(2:3,5) = 0;
  S = fsparse(i,j,ss);
  ok = ok & nnz(S) == sum(ss(:) ~= 0) & isequal(full(S),ss);
  S1 = fsparse(S);
  S2 = fsparse(ss);
  ok = ok & isequal(S,S1) & isequal(S1,S2) & nnz(S1) == nzmax(S1);

  ss(1,1) = 0; ss(1:2,2) = 0;
  S = fsparse(i,j,ss);
  ok = ok & nnz(S) == sum(ss(:) ~= 0) & isequal(full(S),ss);
  S1 = fsparse(S);
  S2 = fsparse(ss);
  ok = ok & isequal(S,S1) & isequal(S1,S2) & nnz(S1) == nzmax(S1);

  ss(3:4,2) = 0;
  S = fsparse(i,j,ss);
  ok = ok & nnz(S) == sum(ss(:) ~= 0) & isequal(full(S),ss);
  S1 = fsparse(S);
  S2 = fsparse(ss);
  ok = ok & isequal(S,S1) & isequal(S1,S2) & nnz(S1) == nzmax(S1);

  ss(1,3) = 0;
  S = fsparse(i,j,ss);
  ok = ok & nnz(S) == sum(ss(:) ~= 0) & isequal(full(S),ss);
  S1 = fsparse(S);
  S2 = fsparse(ss);
  ok = ok & isequal(S,S1) & isequal(S1,S2) & nnz(S1) == nzmax(S1);

  ss(3,4) = 0;
  S = fsparse(i,j,ss);
  ok = ok & nnz(S) == sum(ss(:) ~= 0) & isequal(full(S),ss);
  S1 = fsparse(S);
  S2 = fsparse(ss);
  ok = ok & isequal(S,S1) & isequal(S1,S2) & nnz(S1) == nzmax(S1);

  ss(:) = 0;
  S = fsparse(i,j,ss);
  ok = ok & nnz(S) == sum(ss(:) ~= 0) & isequal(full(S),ss);
  S1 = fsparse(S);
  S2 = fsparse(ss);
  ok = ok & isequal(S,S1) & isequal(S1,S2) & nzmax(S1) == 1;

  ss = c*reshape(1:20,4,5);
  ss(1:2:end) = 0;
  S = fsparse(i,j,ss);
  ok = ok & nnz(S) == sum(ss(:) ~= 0) & isequal(full(S),ss);
  S1 = fsparse(S);
  S2 = fsparse(ss);
  ok = ok & isequal(S,S1) & isequal(S1,S2) & nnz(S1) == nzmax(S1);

  ss = c*reshape(1:20,4,5);
  ss(1:2:end) = 0;
  S = fsparse(i,j,ss);
  ok = ok & nnz(S) == sum(ss(:) ~= 0) & isequal(full(S),ss);
  S1 = fsparse(S);
  S2 = fsparse(ss);
  ok = ok & isequal(S,S1) & isequal(S1,S2) & nnz(S1) == nzmax(S1);

  ss = c*reshape(1:20,4,5);
  ss(1:3:end) = 0;
  S = fsparse(i,j,ss);
  ok = ok & nnz(S) == sum(ss(:) ~= 0) & isequal(full(S),ss);
  S1 = fsparse(S);
  S2 = fsparse(ss);
  ok = ok & isequal(S,S1) & isequal(S1,S2) & nnz(S1) == nzmax(S1);

  % sparser case
  siz = [9 7];
  i = [1 2 6 8 3 3 7 8 1 2 4 5 4 4 5 9];
  j = [2 3 4 7 2 3 5 6 7 1 4 2 2 3 4 6];
  ss = c*[1:16];

  ixout = {[] [10 9] [11 15 3] [8 16 4] [2 6 7] [13 12] [1 5 14]};
  for k = 1: length(ixout)
    ss(ixout{k}) = 0;
    S1 = fsparse(i,j,ss,siz);
    S2 = sparse(i,j,ss,siz(1),siz(2));
    ok = ok & nnz(S1) == sum(ss(:) ~= 0) & isequal(full(S1),full(S2));
    S3 = fsparse(S1);
    ok = ok & nzmax(S3) == max(nnz(S3),1);
  end
end

% mixed complex/real case
rand('state',7891);
i = [1:10]'; j = 1:17;
n = 10*17; nn = 1:n;
s = reshape(complex(mod(nn,2).*rand(1,n),mod(nn,3).*rand(1,n)),10,17);
S = fsparse(i,j,s);
ok = ok & nnz(S) == sum(s(:) ~= 0);
S1 = fsparse(S);
ok = ok & nzmax(S1) == nnz(S1) & isequal(S1,S);

%--------------------------------------------------------------------------
function ok = l_spin13
%L_SPIN13 (FSETOP) Illegal input.

ok = 1;

try
  c = fsetop('unique');
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsetop:e1',li);
end

try
  c = fsetop(1:6,1);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsetop:e2',li);
end

try
  c = fsetop('fo',1);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsetop:e3',li);
end

try
  c = fsetop('lalalalalalalalallalalalalaallala',1);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsetop:e3',li);
end

try
  c = fsetop('union',ones(3,3,3),ones(2,1));
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsetop:e4',li);
end

try
  c = fsetop('union',sparse(1:3,1:3,1),ones(2,1));
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsetop:e4',li);
end

try
  c = fsetop('union',1+i,ones(2,1));
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsetop:e4',li);
end

try
  c = fsetop('union',2,1+i);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsetop:e4',li);
end

try
  c = fsetop('union',1:3);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsetop:e5',li);
end

try
  c = fsetop('union',{'foo'});
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsetop:e5',li);
end

try
  c = fsetop('union',1:3,3,2,1);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsetop:e5',li);
end

try
  c = fsetop('union',1:3,'foba');
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsetop:e6',li);
end

try
  c = fsetop('union',{'bba'},'foba');
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsetop:e6',li);
end

try
  c = fsetop('union',ones(2,3),zeros(3,4));
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsetop:e7',li);
end

try
  c = fsetop('unique',ones(2,3),zeros(3,4));
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsetop:e8',li);
end

try
  c = fsetop('check',ones(2,3),zeros(3,4));
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsetop:e8',li);
end

try
  c = fsetop('unique',{'eee'},{'fff'});
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsetop:e8',li);
end

try
  c = fsetop('unique',{1:3 struct('ba',1)});
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsetop:e11',li);
end

try
  c = fsetop('unique',{char(zeros(1,2,2)) 1+2i});
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsetop:e11',li);
end

try
  c = fsetop('union',{@sin},{@cos})
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsetop:e11',li);
end

try
  c = fsetop('setxor',{'bbb'},{char(zeros(1,2,2)) {'bb'}});
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsetop:e11',li);
end

try
  c = fsetop('setdiff',{'llala'},{1:3 sparse(1,1,1)});
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsetop:e11',li);
end

try
  [c,ia] = fsetop('check',{'llala'});
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsetop:e12',li);
end

try
  [c,ii,j,k,l] = fsetop('intersect',{'llala'},{'aaaa'});
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsetop:e9',li);
end

try
  [c,ii,j,l] = fsetop('setdiff',{'llala'},{'aaaa'});
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsetop:e10',li);
end

try
  [c,ii,k,l] = fsetop('setxor',1:3,1:4);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsetop:e9',li);
end

try
  [c,ii,k,l,m] = fsetop('intersect',1:3,1:4);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsetop:e9',li);
end

try
  [c,ii,k,l,m,n] = fsetop('union',1:3,1:4);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsetop:e13',li);
end

try
  [ii,j,l] = fsetop('ismember',{'laklkalka' 'llala'},{'aaaa'});
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsetop:e10',li);
end

try
  c = fsetop('union',1:3,{'' '' 1});
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsetop:e6',li);
end

try
  c = fsetop('setdiff',{1:3 1:4},'');
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('fsetop:e6',li);
end

%--------------------------------------------------------------------------
function ok = l_spin14
%L_SPIN14 (FSETOP) Some small cases.

ok = 1;

a = [1 3 2 3 4 2];
b = [1 4 6 4 8 1];
A = {int16([1 1]) [3 3] int8([2 2 2]) [3 3] [4] int8([2 2 2])};
B = {int16([1 1]) [4] uint32([6 6 6 6]) [4] single([8 8 8 8 8]) int16([1 1])};

% intersect
c = fsetop('intersect',a,b);
ok = ok && all(c == [1 4]);
[c,ia] = fsetop('intersect',a,b);
ok = ok && all(c == [1 4]) && all(ia == [1 5]);
[c,ia,ib] = fsetop('intersect',a,b);
ok = ok && all(c == [1 4]) && all(ia == [1 5]) && all(ib == [1 2]);
C = fsetop('intersect',A,B);
ok = ok && isequal(C,{int16([1 1]) [4]});
[C,ia] = fsetop('intersect',A,B);
ok = ok && isequal(C,{int16([1 1]) [4]}) && all(ia == [1 5]);
[C,ia,ib] = fsetop('intersect',A,B);
ok = ok && isequal(C,{int16([1 1]) [4]}) && ...
     all(ia == [1 5]) && all(ib == [1 2]);

% setdiff
c = fsetop('setdiff',a,b);
ok = ok && all(c == [3 2]);
[c,ia] = fsetop('setdiff',a,b);
ok = ok && all(c == [3 2]) && all(ia == [2 3]);
C = fsetop('setdiff',A,B);
ok = ok && isequal(C,{[3 3] int8([2 2 2])});
[C,ia] = fsetop('setdiff',A,B);
ok = ok && isequal(C,{[3 3] int8([2 2 2])}) && all(ia == [2 3]);

% setxor
c = fsetop('setxor',a,b);
ok = ok && all(c == [3 2 6 8]);
[c,ia] = fsetop('setxor',a,b);
ok = ok && all(c == [3 2 6 8]) && all(ia == [2 3]);
[c,ia,ib] = fsetop('setxor',a,b);
ok = ok && all(c == [3 2 6 8]) && all(ia == [2 3]) && all(ib == [3 5]);
C = fsetop('setxor',A,B);
ok = ok && isequal(C,{[3 3] int8([2 2 2]) uint32([6 6 6 6]) ...
                    single([8 8 8 8 8])});
[C,ia] = fsetop('setxor',A,B);
ok = ok && isequal(C,{[3 3] int8([2 2 2]) uint32([6 6 6 6]) ...
                    single([8 8 8 8 8])}) && all(ia == [2 3]);
[C,ia,ib] = fsetop('setxor',A,B);
ok = ok && isequal(C,{[3 3] int8([2 2 2]) uint32([6 6 6 6]) ...
                    single([8 8 8 8 8])}) && ...
     all(ia == [2 3]) && all(ib == [3 5]);

% union
c = fsetop('union',a,b);
ok = ok && all(c == [1 3 2 4 6 8]);
[c,ia] = fsetop('union',a,b);
ok = ok && all(c == [1 3 2 4 6 8]) && all(ia == [1 2 3 5]);
[c,ia,ib] = fsetop('union',a,b);
ok = ok && all(c == [1 3 2 4 6 8]) && all(ia == [1 2 3 5]) && ...
     all(ib == [3 5]);
[c,ia,ib,ja] = fsetop('union',a,b);
ok = ok && all(c == [1 3 2 4 6 8]) && all(ia == [1 2 3 5]) && ...
     all(ib == [3 5]) && all(c(ja) == a);
[c,ia,ib,ja,jb] = fsetop('union',a,b);
ok = ok && all(c == [1 3 2 4 6 8]) && all(ia == [1 2 3 5]) && ...
     all(ib == [3 5]) && all(c(ja) == a) && all(c(jb) == b);
C = fsetop('union',A,B);
CC = {int16([1 1]) [3 3] int8([2 2 2]) [4] ...
      uint32([6 6 6 6]) single([8 8 8 8 8])};
ok = ok && isequal(C,CC);
[C,ia] = fsetop('union',A,B);
ok = ok && isequal(C,CC) ...
     && all(ia == [1 2 3 5]);
[C,ia,ib] = fsetop('union',A,B);
ok = ok && isequal(C,CC) ...
     && all(ia == [1 2 3 5]) && all(ib == [3 5]); 
[C,ia,ib,ja] = fsetop('union',A,B);
ok = ok && isequal(C,CC) ...
     && all(ia == [1 2 3 5]) && all(ib == [3 5]) && ...
     isequal(C(ja),A); 
[C,ia,ib,ja,jb] = fsetop('union',A,B);
ok = ok && isequal(C,CC) ...
     && all(ia == [1 2 3 5]) && all(ib == [3 5]) && ...
     isequal(C(ja),A) && isequal(C(jb),B); 

% unique
c = fsetop('unique',a);
ok = ok && all(c == [1 3 2 4]);
[c,ia] = fsetop('unique',a);
ok = ok && all(c == [1 3 2 4]) && all(ia == [1 2 3 5]);
[c,ia,ib] = fsetop('unique',a);
ok = ok && all(c == [1 3 2 4]) && all(ia == [1 2 3 5]) && ...
     all(ib == [1 2 3 2 4 3]);
C = fsetop('unique',A);
ok = ok && isequal(C,{int16([1 1]) [3 3] int8([2 2 2]) [4]});
[C,ia] = fsetop('unique',A);
ok = ok && isequal(C,{int16([1 1]) [3 3] int8([2 2 2]) [4]}) ...
     && all(ia == [1 2 3 5]);
[C,ia,ib] = fsetop('unique',A);
ok = ok && isequal(C,{int16([1 1]) [3 3] int8([2 2 2]) [4]}) ...
     && all(ia == [1 2 3 5]) && ...
     all(ib == [1 2 3 2 4 3]);

% ismember
ia = fsetop('ismember',a,b);
ok = ok && all(ia == [1 0 0 0 1 0]);
[ia,ib] = fsetop('ismember',a,b);
ok = ok && all(ia == [1 0 0 0 1 0]) && all(ib == [1 0 0 0 2 0]);
ia = fsetop('ismember',A,B);
ok = ok && all(ia == [1 0 0 0 1 0]);
[ia,ib] = fsetop('ismember',A,B);
ok = ok && all(ia == [1 0 0 0 1 0]) && all(ib == [1 0 0 0 2 0]);

% from help...just see that it doesn't crash...
% intersection of integers
a = int8(ceil(3*rand(2,10))); b = int8(ceil(3*rand(2,10)));
aandb = fsetop('intersect',a,b);
 
% unique strings in cell-array
strs = {'foobar' 'foo' 'foobar' 'bar' 'barfoo'};
strsunq = fsetop('unique',strs);

% find indices ix in b = a(ix)
a = randperm(6); b = ceil(6*rand(1,6));
[foo,ix] = fsetop('ismember',b,a);

% remove indices from set of indices
a = {[1 2 4] [2 3 5 6] [4 2 1] [2 3 5 6]' int32([1 2 4])};
b = {[4 2 1]};
c = fsetop('setdiff',a,b);

% union of structs
a = struct('foo',1,'bar',NaN,'foobar','hello');
b = struct('faa',2,'bar',Inf,'foobar','goodbye');
[cf,ia,ib] = fsetop('union',fieldnames(a),fieldnames(b));
ac = struct2cell(a); bc = struct2cell(b);
c = cell2struct([ac(ia); bc(ib)],cf);

% a single 32-bit checksum from strings
s = {'check','intersect','setdiff','setxor', ...
     'union','unique','ismember'};
c = fsetop('check',fsetop('check',s)');

%--------------------------------------------------------------------------
function ok = l_spin15
%L_SPIN15 (FSETOP) More serious testing of the set operations.

ok = 1;

% test different types
types = {@double @single @char ...
         @int8 @uint8 @int16 @uint16 ...
         @int32 @uint32 @int64 @uint64};

% to debug, remove operations here
ops = {'intersect' 'setdiff' 'setxor' 'union' 'unique' 'ismember'};
siz = size(ops,2);

% scalar case, various types
a = [1 3 2 1 8 5 3 4 6 3 6 7 8 4 3 2];
b = [9 7 8 5 6 7 6 4 9 8 6];

for i = siz
  op = ops{i};

  switch i
   case 1
    c =  [8 5 4 6 7];
    ia = [5 6 8 9 12];
    ib = [3 4 8 5 2];

    ok = ok && all(c == a(ia)) && all(c == b(ib));
   case 2
    c =  [1 3 2];
    ia = [1 2 3];

    ok = ok && all(c == a(ia));
   case 3
    c = [1 3 2 9];
    ia = [1 2 3];
    ib = 1;

    ok = ok && all(c == [a(ia) b(ib)]);
   case 4
    c = [1 3 2 8 5 4 6 7 9];
    ia = [1 2 3 5 6 8 9 12];
    ib = 1;

    ok = ok && all(c == [a(ia) b(ib)]);
   case 5
    c = [1 3 2 8 5 4 6 7];
    ia = [1 2 3 5 6 8 9 12];
    ib = [1 2 3 1 4 5 2 6  7 2 7 8 4 6 2 3];
 
    ok = ok && all(c == a(ia)) && all(a == c(ib));
   case 6
    ia = logical([0 0 0 0 1 1 0 1 1 0 1 1 1 1 0 0]);
    ib = [0 0 0 0 3 4 0 8 5 0 5 2 3 8 0 0];

    ok = ok && all(a(ia) == b(ib(ib ~= 0)));
  end

  for t = types
    if i == 2
      aa = t{1}(a); bb = t{1}(b);
      cc = fsetop(op,aa,bb);
      ok = ok && all(cc == t{1}(c));
      [cc,iaa] = fsetop(op,aa,bb);
      ok = ok && all(cc == t{1}(c)) && all(iaa == ia);
    elseif i == 5
      aa = t{1}(a);
      cc = fsetop(op,aa);
      ok = ok && all(cc == t{1}(c));
      [cc,iaa] = fsetop(op,aa);
      ok = ok && all(cc == t{1}(c)) && all(iaa == ia);
      [cc,iaa,ibb] = fsetop(op,aa);
      ok = ok && all(cc == t{1}(c)) && all(iaa == ia) && ...
           all(ibb == ib);
    elseif i == 6
      aa = t{1}(a); bb = t{1}(b);
      iaa = fsetop(op,aa,bb);
      ok = ok && all(iaa == ia);
      [iaa,ibb] = fsetop(op,aa,bb);
      ok = ok && all(iaa == ia) && all(ibb == ib);
    else
      aa = t{1}(a); bb = t{1}(b);
      cc = fsetop(op,aa,bb);
      ok = ok && all(cc == t{1}(c));
      [cc,iaa] = fsetop(op,aa,bb);
      ok = ok && all(cc == t{1}(c)) && all(iaa == ia);
      [cc,iaa,ibb] = fsetop(op,aa,bb);
      ok = ok && all(cc == t{1}(c)) && all(iaa == ia) && ...
           all(ibb == ib);
      if i == 4
        [cc,iaa,ibb,jaa] = fsetop(op,aa,bb);
        ok = ok && all(aa(ja) == cc);
        [cc,iaa,ibb,jaa,jbb] = fsetop(op,aa,bb);
        ok = ok && all(aa(ja) == cc) && all(bb(jb) == cc);
      end
    end
  end
end

% logicals
for i = siz
  op = ops{i};
  a = logical([0 0 0 1]);
  b = logical([1 1 1]);
  if i == 2
    [cc,iaa] = fsetop(op,a,b);
  elseif i == 5
    [cc,iaa,ibb] = fsetop(op,a);
  elseif i == 6
    [iaa,ibb] = fsetop(op,a,b);
  else
    [cc,iaa,ibb] = fsetop(op,a,b);
  end
  switch i
   case 1
    ok = ok && cc && iaa == 4 && ibb == 1;
   case 2
    ok = ok && ~cc && iaa == 1;
   case 3
    ok = ok && ~cc && iaa == 1 && all(size(ibb) == [1 0]);
   case 4
    ok = ok && all(cc == [0 1]) && all(iaa == [1 4]) && ...
         all(size(ibb) == [1 0]);
   case 5
    ok = ok && all(cc == [0 1]) && all(iaa == [1 4]) && ...
         all(ibb == [1 1 1 2]);
   case 6
    ok = ok && all(iaa == logical([0 0 0 1])) && all(ibb == [0 0 0 1]);
  end
end

% columns
a = [1 3 2 1 8 5 3 4 6 3 6 7 8 4 3 2; ...
     1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2];
b = [9 7 8 5 6 7 6 4 9 8 6; ...
     1 2 1 2 1 0 1 0 1 2 1];

for i = siz
  op = ops{i};

  switch i
   case 1
    c = [8 5 6 7; ...
         1 2 1 2];
    ia = [5 6 9 12];
    ib = [3 4 5 2];

    ok = ok && norm(c-a(:,ia),1) < eps && norm(c-b(:,ib),1) < eps;
   case 2
    c = [1 3 2 1 3 4 2; ...
         1 2 1 2 1 2 2];
    ia = [1 2 3 4 7 8 16];

    ok = ok && norm(c-a(:,ia),1) < eps;
   case 3
    c = [1 3 2 1 3 4 2 9 7 4 8; ...
         1 2 1 2 1 2 2 1 0 0 2];
    ia = [1 2 3 4 7 8 16];
    ib = [1 6 8 10];

    ok = ok && norm(c-[a(:,ia) b(:,ib)],1) < eps;
   case 4
    c = [1 3 2 1 8 5 3 4 6 7 2 9 7 4 8; ...
         1 2 1 2 1 2 1 2 1 2 2 1 0 0 2];
    ia = [1 2 3 4 5 6 7 8 9 12 16];
    ib = [1 6 8 10];
 
    ok = ok && norm(c-[a(:,ia) b(:,ib)],1) < eps;
   case 5
    c = [1 3 2 1 8 5 3 4 6 7 2; ...
         1 2 1 2 1 2 1 2 1 2 2];
    ia = [1 2 3 4 5 6 7 8 9 12 16];
    ib = [1 2 3 4 5 6 7 8 9 2  9  10 5 8 7 11];
 
    ok = ok && norm(c-a(:,ia),1) < eps && norm(a-c(:,ib),1) < eps;
   case 6
    ia = logical([0 0 0 0 1 1 0 0 1 0 1 1 1 0 0 0]);
    ib = [0 0 0 0 3 4 0 0 5 0 5 2 3 0 0 0];

    ok = ok && norm(a(:,ia)-b(:,ib(ib ~= 0)),1) < eps;
  end

  for t = types
    if i == 2
      aa = t{1}(a); bb = t{1}(b);
      cc = fsetop(op,aa,bb);
      ok = ok && norm(double(cc)-double(t{1}(c)),1) < eps;
      [cc,iaa] = fsetop(op,aa,bb);
      ok = ok && norm(double(cc)-double(t{1}(c)),1) < eps && ...
           all(iaa == ia);
    elseif i == 5
      aa = t{1}(a); bb = t{1}(b);
      cc = fsetop(op,aa);
      ok = ok && norm(double(cc)-double(t{1}(c)),1) < eps;
      [cc,iaa] = fsetop(op,aa);
      ok = ok && norm(double(cc)-double(t{1}(c)),1) < eps && ...
           all(iaa == ia);
      [cc,iaa,ibb] = fsetop(op,aa);
      ok = ok && norm(double(cc)-double(t{1}(c)),1) < eps && ...
           all(iaa == ia) && all(ibb == ib);
    elseif i == 6
      aa = t{1}(a); bb = t{1}(b);
      iaa = fsetop(op,aa,bb);
      ok = ok && all(iaa == ia);
      [iaa,ibb] = fsetop(op,aa,bb);
      ok = ok && all(iaa == ia) && all(ibb == ib);
    else
      aa = t{1}(a); bb = t{1}(b);
      cc = fsetop(op,aa,bb);
      ok = ok && norm(double(cc)-double(t{1}(c)),1) < eps;
      [cc,iaa] = fsetop(op,aa,bb);
      ok = ok && norm(double(cc)-double(t{1}(c)),1) < eps && ...
           all(iaa == ia);
      [cc,iaa,ibb] = fsetop(op,aa,bb);
      ok = ok && norm(double(cc)-double(t{1}(c)),1) < eps && ...
           all(iaa == ia) && all(ibb == ib);
      if i == 4
        [cc,iaa,ibb,jaa] = fsetop(op,aa,bb);
        ok = ok && norm(double(cc(:,ja))-double(aa)) < eps;
        [cc,iaa,ibb,jaa,jbb] = fsetop(op,aa,bb);
        ok = ok && norm(double(cc(:,ja))-double(aa)) < eps && ...
             norm(double(cc(:,jb))-double(bb)) < eps
      end
    end
  end
end

% columns with logicals
for i = siz
  op = ops{i};
  switch i
   case 1
    [cc,iaa,ibb] = fsetop(op,...
                          logical([0 0 0 1; ...
                        0 0 1 1]), ...
                          logical([1 1 1; ...
                        0 0 1]));
    ok = ok && all(cc == true(2,1)) && iaa == 4 && ibb == 3;
   case 2
    [cc,iaa] = fsetop(op,...
                      logical([0 1 0 1; ...
                        0 0 1 1]), ...
                      logical([1 0 1; ...
                        0 0 1]));
    ok = ok && all(cc == logical([0; 1])) && iaa == 3;
   case 3
    [cc,iaa,ibb] = fsetop(op,...
                          logical([0 0 0 1; ...
                        0 0 1 1]), ...
                          logical([1 1 1; ...
                        0 0 1]));
    ok = ok && all(all(cc == logical([0 0 1; ...
                        0 1 0]))) && ...
         all(iaa == [1 3]) && ibb == 1;
   case 4
    [cc,iaa,ibb] = fsetop(op,...
                          logical([0 0 0 1; ...
                        0 0 1 1]), ...
                          logical([1 1 1; ...
                        0 0 1]));
    ok = ok && all(all(cc == logical([0 0 1 1; ...
                        0 1 1 0]))) && ...
         all(iaa == [1 3 4]) && ibb == 1;
   case 5
    [cc,iaa,ibb] = fsetop(op,...
                          logical([0 0 0 1; ...
                        0 0 1 1]));
    ok = ok && all(all(cc == logical([0 0 1; ...
                        0 1 1]))) && ...
         all(iaa == [1 3 4]) && all(ibb == [1 1 2 3]);
   case 6
    [iaa,ibb] = fsetop(op,...
                       logical([0 1 0 1; ...
                        0 0 1 1]), ...
                       logical([1 0 1; ...
                        0 0 1]));
    ok = ok && all(iaa == logical([1 1 0 1])) && all(ibb == [2 1 0 3]);
  end
end

%--------------------------------------------------------------------------
function ok = l_spin16
%L_SPIN16 (FSETOP) Cell-arrays.

ok = 1;

a = {'aba' 'gag' 'aca' 's ' 'quattro' ...
     's' 'st' 'ss' ' ' 'gag2' ' ' 'aba' 'foo' ...
     'baa ' ' baa' ' baa ' '_baa' '*=%&?#' 'seize' ...
     'vingt-et-huite' 'bab' ' ' '  ' 'oof' 'baa_'};
b = {'seize' '_seize' 'bab' ' ' '*' 'quattre' 'Quattro' 'gag' ...
     's s ' 'vingt' '*=' 'neneh' 'nen' 'pour-quoi?' 'parce-que!' ...
     'n''a?' 'ss ' '  ' 'Bab'};

c1 = fsetop('intersect',a,b);
c2 = intersect(a,b);
ok = ok && all(fsetop('ismember',c1,c2)) && all(size(c1) == size(c2));

[c1,ia,ib] = fsetop('intersect',a,b);
ok = ok && isequal(c1,a(ia)) && isequal(c1,b(ib));

c1 = fsetop('intersect',b,a);
c2 = intersect(b,a);
ok = ok && all(fsetop('ismember',c1,c2)) && all(size(c1) == size(c2));

[c1,ia,ib] = fsetop('intersect',b,a);
ok = ok && isequal(c1,b(ia)) && isequal(c1,a(ib));

c1 = fsetop('setdiff',a,b);
c2 = setdiff(a,b);
ok = ok && all(fsetop('ismember',c1,c2)) && all(size(c1) == size(c2));

[c1,ia] = fsetop('setdiff',a,b);
ok = ok && isequal(c1,a(ia));

c1 = fsetop('setdiff',b,a);
c2 = setdiff(b,a);
ok = ok && all(fsetop('ismember',c1,c2)) && all(size(c1) == size(c2));

[c1,ia] = fsetop('setdiff',b,a);
ok = ok && isequal(c1,b(ia));

c1 = fsetop('setxor',a,b);
c2 = setxor(a,b);
ok = ok && all(fsetop('ismember',c1,c2)) && all(size(c1) == size(c2));

[c1,ia,ib] = fsetop('setxor',a,b);
ok = ok && isequal(c1,[a(ia) b(ib)]);

c1 = fsetop('setxor',b,a);
c2 = setxor(b,a);
ok = ok && all(fsetop('ismember',c1,c2)) && all(size(c1) == size(c2));

[c1,ia,ib] = fsetop('setxor',b,a);
ok = ok && isequal(c1,[b(ia) a(ib)]);

c1 = fsetop('union',a,b);
c2 = union(a,b);
ok = ok && all(fsetop('ismember',c1,c2)) && all(size(c1) == size(c2));

[c1,ia,ib] = fsetop('union',a,b);
ok = ok && isequal(c1,[a(ia) b(ib)]);

[c1,ia,ib,ja,jb] = fsetop('union',a,b);
ok = ok && isequal(c1,[a(ia) b(ib)]) && isequal(c1(ja),a) && isequal(c1(jb),b);
[c1,ib,ia,jb,ja] = fsetop('union',b,a);
ok = ok && isequal(c1,[b(ib) a(ia)]) && isequal(c1(ja),a) && isequal(c1(jb),b);

c1 = fsetop('union',b,a);
c2 = union(b,a);
ok = ok && all(fsetop('ismember',c1,c2)) && all(size(c1) == size(c2));

[c1,ia,ib] = fsetop('union',b,a);
ok = ok && isequal(c1,[b(ia) a(ib)]);

c1 = fsetop('unique',a);
c2 = unique(a);
ok = ok && all(fsetop('ismember',c1,c2)) && all(size(c1) == size(c2));

[c1,ia,ib] = fsetop('unique',a);
ok = ok && isequal(c1,a(ia)) && isequal(a,c1(ib));

c1 = fsetop('unique',b);
c2 = unique(b);
ok = ok && all(fsetop('ismember',c1,c2)) && all(size(c1) == size(c2));

[c1,ia,ib] = fsetop('unique',b);
ok = ok && isequal(c1,b(ia)) && isequal(b,c1(ib));

c1 = fsetop('ismember',a,b);
c2 = ismember(a,b);
% ismember does not differ between 'ss' and 'ss '
ok = ok && all(c1([1:7 9:end]) == c2([1:7 9:end]));

[ia,ib] = fsetop('ismember',a,b);
ok = ok && isequal(a(ia),b(ib(ib ~= 0)));

c1 = fsetop('ismember',b,a);
c2 = ismember(b,a);
ok = ok && all(c1([1:16 18:end]) == c2([1:16 18:end]));

[ia,ib] = fsetop('ismember',b,a);
ok = ok && isequal(b(ia),a(ib(ib ~= 0)));

% various mixed cases
A = {[1:3]' int8(1) int16(1) uint32(1) ones(2,2) ones(2,1,2) ...
     int8(ones(2,1,2)) ones(2,2) int32(2) 'end'}';
B = {[1:4] int16(1) int8(2) uint16(1) uint32(1) ones(1,2) ...
     ones(2,1,3) 'begin' int16(2) int8(ones(2,1,2))}';

% check of indices only
[C,ia,ib] = fsetop('intersect',A,B);
ok = ok && isequal(C',A(ia)) && isequal(C',B(ib));
[C,ia] = fsetop('setdiff',A,B);
ok = ok && isequal(C',A(ia));
[C,ia,ib] = fsetop('setxor',A,B);
ok = ok && isequal(C',[A(ia); B(ib)]);
[C,ia,ib] = fsetop('union',A,B);
ok = ok && isequal(C',[A(ia); B(ib)]);
[C,ia,ib,ja,jb] = fsetop('union',A,B);
ok = ok && isequal(C',[A(ia); B(ib)]) && ...
     isequal(C(ja)',A) && isequal(C(jb)',B);
[C,ia,ib] = fsetop('unique',A);
ok = ok && isequal(C',A(ia)) && isequal(A,C(ib)');
[ia,ib] = fsetop('ismember',A,B);
ok = ok && isequal(A(ia),B(ib(ia)));
 
%--------------------------------------------------------------------------
function ok = l_spin17
%L_SPIN17 (FSETOP) Empty cases.

ok = 1;

[c,ia,ib] = fsetop('intersect',1:3,4:6);
ok = ok && all([size(c) == [1 0] size(ia) == [1 0] size(ib) == [1 0]]);
[c,ia,ib] = fsetop('intersect',zeros(2,0),[1:2]');
ok = ok && all([size(c) == [2 0] size(ia) == [1 0] size(ib) == [1 0]]);
[c,ia,ib] = fsetop('intersect',reshape(1:6,3,2),zeros(3,0));
ok = ok && all([size(c) == [3 0] size(ia) == [1 0] size(ib) == [1 0]]);
[c,ia,ib] = fsetop('intersect',[],[]);
ok = ok && all([size(c) == [0 0] size(ia) == [1 0] size(ib) == [1 0]]);
% from GNATS
[c,ia,ib] = fsetop('intersect',zeros(0,3),zeros(0,5));
ok = ok && all(size(c) == [0 1]) && ia == 1 && ib == 1;

[c,ia,ib] = fsetop('intersect',{1 2 3},{4 5 6});
ok = ok && all([size(c) == [1 0] size(ia) == [1 0] size(ib) == [1 0]]);
[c,ia,ib] = fsetop('intersect',cell(2,0),{1 2}');
ok = ok && all([size(c) == [1 0] size(ia) == [1 0] size(ib) == [1 0]]);

[c,ia] = fsetop('setdiff',1:3,1:4);
ok = ok && all([size(c) == [1 0] size(ia) == [1 0]]);
[c,ia] = fsetop('setdiff',zeros(2,0),reshape(1:4,2,2));
ok = ok && all([size(c) == [2 0] size(ia) == [1 0]]);
[c,ia] = fsetop('setdiff',zeros(5,0),zeros(5,0));
ok = ok && all([size(c) == [5 0] size(ia) == [1 0]]);
[c,ia] = fsetop('setdiff',[],[]);
ok = ok && all([size(c) == [0 0] size(ia) == [1 0]]);
% from GNATS
[c,ia] = fsetop('setdiff',zeros(0,3),zeros(0,5));
ok = ok && all([size(c) == [0 0] size(ia) == [1 0]]);

[c,ia] = fsetop('setdiff',{1 2 3},{1 2 3 4});
ok = ok && all([size(c) == [1 0] size(ia) == [1 0]]);
[c,ia] = fsetop('setdiff',cell(2,0),reshape({1 2 3 4},2,2));
ok = ok && all([size(c) == [1 0] size(ia) == [1 0]]);

[c,ia,ib] = fsetop('setxor',1:3,1:3);
ok = ok && all([size(c) == [1 0] size(ia) == [1 0] size(ib) == [1 0]]);
[c,ia,ib] = fsetop('setxor',zeros(2,0),zeros(2,0));
ok = ok && all([size(c) == [2 0] size(ia) == [1 0] size(ib) == [1 0]]);
[c,ia,ib] = fsetop('setxor',[],[]);
ok = ok && all([size(c) == [0 0] size(ia) == [1 0] size(ib) == [1 0]]);
% from GNATS
[c,ia,ib] = fsetop('setxor',zeros(0,3),zeros(0,5));
ok = ok && all([size(c) == [0 0] size(ia) == [1 0] size(ib) == [1 0]]);

[c,ia,ib] = fsetop('setxor',{1 2 3},{1 2 3});
ok = ok && all([size(c) == [1 0] size(ia) == [1 0] size(ib) == [1 0]]);
[c,ia,ib] = fsetop('setxor',cell(2,0),cell(2,0));
ok = ok && all([size(c) == [1 0] size(ia) == [1 0] size(ib) == [1 0]]);

[c,ia,ib] = fsetop('union',zeros(3,0),zeros(3,0));
ok = ok && all([size(c) == [3 0] size(ia) == [1 0] size(ib) == [1 0]]);
[c,ia,ib] = fsetop('union',[],[]);
ok = ok && all([size(c) == [0 0] size(ia) == [1 0] size(ib) == [1 0]]);
% from GNATS
[c,ia,ib] = fsetop('union',zeros(0,3),zeros(0,5));
ok = ok && all([size(c) == [0 1] size(ib) == [1 0]]) && ia == 1;
[c,ia,ib,ja,jb] = fsetop('union',zeros(0,3),zeros(0,5));
ok = ok && all([size(c) == [0 1] size(ib) == [1 0]]) && ia == 1;
ok = ok && isequal(c(:,ja),zeros(0,3)) && isequal(c(:,jb),zeros(0,5));

[c,ia,ib] = fsetop('union',cell(3,0),cell(3,0));
ok = ok && all([size(c) == [1 0] size(ia) == [1 0] size(ib) == [1 0]]);
[c,ia,ib] = fsetop('union',{},{});
ok = ok && all([size(c) == [1 0] size(ia) == [1 0] size(ib) == [1 0]]);

[c,ia,ib] = fsetop('unique',zeros(2,0));
ok = ok && all([size(c) == [2 0] size(ia) == [1 0] size(ib) == [1 0]]);
[c,ia,ib] = fsetop('unique',[]);
ok = ok && all([size(c) == [0 0] size(ia) == [1 0] size(ib) == [1 0]]);
% from GNATS
[c,ia,ib] = fsetop('unique',zeros(0,3));
ok = ok && all(size(c) == [0 1]) && ia == 1 && isequal(ib,[1 1 1]);

[c,ia,ib] = fsetop('unique',cell(2,0));
ok = ok && all([size(c) == [1 0] size(ia) == [1 0] size(ib) == [1 0]]);
[c,ia,ib] = fsetop('unique',{});
ok = ok && all([size(c) == [1 0] size(ia) == [1 0] size(ib) == [1 0]]);

[ia,ib] = fsetop('ismember',zeros(2,0),reshape(1:6,2,3));
ok = ok && all([size(ia) == [1 0] size(ib) == [1 0]]);
[ia,ib] = fsetop('ismember',[],[]);
ok = ok && all([size(ia) == [1 0] size(ib) == [1 0]]);
% from GNATS
[ia,ib] = fsetop('ismember',zeros(0,5),zeros(0,3));
ok = ok && isequal(ia,[1 1 1 1 1]) && isequal(ib,[1 1 1 1 1]);

[ia,ib] = fsetop('ismember',cell(2,0),{1 2 3 4 5 6});
ok = ok && all([size(ia) == [1 0] size(ib) == [1 0]]);
[ia,ib] = fsetop('ismember',[],[]);
ok = ok && all([size(ia) == [1 0] size(ib) == [1 0]]);

c = fsetop('check',[]);
ok = ok && all(size(c) == [1 0]);
c = fsetop('check',{});
ok = ok && all(size(c) == [1 0]);

%--------------------------------------------------------------------------
function ok = l_spin18
%L_SPIN18 (CLENSHAW) Illegal input.

ok = 1;

try
  clenshaw(1,2,3,4,5);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshaw:e1',li);
end

try
  clenshaw(1);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshaw:e1',li);
end

try
  [a,b] = clenshaw(1,2,3);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshaw:e1',li);
end

try
  clenshaw('lala',1);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshaw:e2',li);
end

try
  clenshaw(1,1,{1});
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshaw:e2',li);
end

try
  clenshaw(ones(1,2),ones(1,3));
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshaw:e3',li);
end

try
  clenshaw(ones(2),ones(3,2),ones(3));
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshaw:e3',li);
end

try
  clenshaw(1,2,3,{1});
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshaw:e1',li);
end

try
  clenshaw(1,2,3,{1});
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshaw:e1',li);
end

try
  clenshaw(1,2,3,'lalalalalalalalaala');
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshaw:e1',li);
end

try
  clenshaw(1,2,'la');
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshaw:e2',li);
end

try
  clenshaw({1 2},{2});
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshaw:e2',li);
end

try
  clenshaw({ones(1,2) ones(1,3)},2,'point');
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshaw:e2',li);
end

try
  clenshaw({ones(1,2) ones(1,2)},2,ones(1,1,3),'point');
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshaw:e1',li);
end

try
  clenshaw(ones(4),ones(4,1),ones(3,4));
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshaw:e4',li);
end

try
  clenshaw({ones(4) ones(3)},{ones(4,1) ones(1,3)},ones(4,2));
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshaw:e2',li);
end

try
  clenshaw(1,ones(2),ones(3,2));
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshaw:e4',li);
end

try
  clenshaw({},{},single(1));
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshaw:e2',li);
end

try
  clenshaw({},{},ones(2,1));
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshaw:e2',li);
end

try
  clenshaw({1 2},{1 2},ones(3,3,3));
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshaw:e2',li);
end

try
  clenshaw({1 2},{1 2},ones(3,3,3,2),'point');
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshaw:e1',li);
end

try
  [a,b,c] = clenshaw({1 2},{1 2});
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshaw:e1',li);
end

try
  [a,b] = clenshaw(1,{1 2 3 4});
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshaw:e1',li);
end

try
  a = clenshaw(ones(2,3),ones(2,3),ones(2,2,3,2));
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshaw:e5',li);
end

try
  a = clenshaw(ones(2,3),ones(2,3),ones(2,3,4));
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshaw:e6',li);
end

%--------------------------------------------------------------------------
function Y = l_recurrence(A,B,flag)
%L_RECURRENCE Mimics the 'recurrence'-syntax of CLENSHAW.

if ~iscell(A) && ~iscell(B)
  sza = size(A);
  szb = size(B);
  sz = max(sza,szb);
  A = frepmat(A,sz./max(1,sza));
  B = frepmat(B,sz./max(1,szb));
  Y = zeros(sz);

  for j = 1:sz(2)
    r1 = 1; r2 = 1;
    for i = 1:sz(1)
      r3 = r2;
      r2 = r1;
      r1 = A(i,j)*r2+B(i,j)*r3;
      Y(i,j) = r1;
    end
  end
else
  if ~iscell(A), A = frepmat({A},size(B)); end
  if ~iscell(B), B = frepmat({B},size(A)); end
  len = size(A,2);
  YY = cell(1,len);
  for i = 1:len
    YY{i} = l_recurrence(A{i},B{i});
  end

  if nargin < 3 || strcmp(flag,'tensor')
    Y = YY{1};
    for i = 2:len
      Y = Y(:)*YY{i}(:).';
    end
    sz = [cellfun('size',YY,1); cellfun('size',YY,2)];
    Y = permute(reshape(Y,sz(:)'),[1:2:2*len 2:2:2*len]);
  else
    Y = YY{1}; n = size(Y,2);
    for i = 2:len
      Z = zeros(size(Y,1)*size(YY{i},1),n);
      for j = 1:n
        Z(:,j) = reshape(Y(:,j)*YY{i}(:,j).',[],1);
      end
      Y = Z;
    end
    Y = reshape(Y,[cellfun('size',YY,1) n]);
  end
end

%--------------------------------------------------------------------------
function ok = l_spin19
%L_SPIN19 (CLENSHAW) Small cases, 1-D syntax.

ok = 1;

% Legendre polynomials
for c = [1 0.54+0.37i 0.8        0.45+0.33i; ...
         1 1          0.34+0.45i 0.45+0.89i]
  x = c(1)*linspace(-1,1,17);
  for n = 1:6
    % 'recurrence'-syntax
    j = [1:n-1]';
    A = [ones(size(x)); (2*j-1)./j*x];
    B = [0; (1-j)./j];
    y1 = clenshaw(A,B);
    y = l_recurrence(A,B);
    ok = ok && norm(y(:)-y1(:),inf) < 1e-15;

    B2 = frepmat(B,size(x));
    y2 = clenshaw(A,B2);
    ok = ok && norm(y1(:)-y2(:),inf) < 1e-15;

    % 'coefficient'-syntax
    C = c(2)*exp(n-1:-1:0)';
    y3 = clenshaw(A,B,C);
    yy = C.'*y1;
    ok = ok && norm(y3-yy,inf) < 1e-12;

    C2 = frepmat(C,[1 1 size(x,2)]);
    y4 = clenshaw(A,B,C2);
    ok = ok && norm(y4-yy,inf) < 1e-12;

    y5 = clenshaw(A,B2,C);
    ok = ok && norm(y5-yy,inf) < 1e-12;

    y6 = clenshaw(A,B2,C2);
    ok = ok && norm(y6-yy,inf) < 1e-12;
  end
end

% from help, just check so that it doesn't crash...
% Fibonacci numbers
n = 10; clenshaw(1,[0; 0; ones(n-2,1)])';

% Legendre polynomials
x = linspace(-1,1); n = 4; j = [0.5 1:n]';
A = (2*j-1)./j*x;
B = (1-j)./j;
y = clenshaw(A,B);

% two functions defined by Chebyshev coefficients
x = linspace(-1,1,50);
n = 4; j = [0; 1; 2*ones(n-1,1)];
A = j*x; B = [1; 0; -ones(n-1,1)];
C = [[1 0.75 0.5 0.25 0.125]' [1 -1.75 1.5 -1.25 2.125]'];

%--------------------------------------------------------------------------
function ok = l_spin20
%L_SPIN20 (CLENSHAW) Random recurrence (1-D), all possibillities.

ok = 1;
rand('state',1429);

for c = [1 0.5+0.3i 1        0.5+0.3i; ...
         1 1        0.3+0.4i 0.4+0.8i];
  for n = 1:6
    A = frepmat(c(1)*rand,[n 14]);
    B = frepmat(c(2)*rand,[n 14]);

    % 'recurrence'-syntax
    y = l_recurrence(A,B);
    f1 = {A A(:,1) A(1,:) A(1,1)};
    f2 = {B B(:,1) B(1,:) B(1,1)};
    for j1 = 1:4
      for j2 = 1:4
        % some cases are not interesting:
        siz = max(size(f1{j1}),size(f2{j2}));
        if all(siz(1) == [n 14])
          yy = clenshaw(f1{j1},f2{j2});
          ok = ok && norm(y-yy,inf) < 5e-14 && all(siz == size(yy));
          ok = ok && isequal(clenshaw(f1{j1},f2(j2)),yy);
        end
      end
    end
  end
end

% 'coefficient'-syntax
for c = [1 0.5+0.3i 1        1        0.5+0.3i 0.3+0.6i 1        0.4+0.4i; ...
         1 1        0.3+0.4i 1        0.4+0.8i 1        0.4+0.5i 0.2+0.7i; ...
         1 1        1        0.8+0.4i 1        0.2+0.4i 0.3+0.7i 0.3+0.9i]
  for n = 1:6
    A = frepmat(c(1)*rand,[n 14]);
    B = frepmat(c(2)*rand,[n 14]);
    C = frepmat(c(3)*rand,[n 1 14]);
    
    y = l_recurrence(A,B);
    y = C(:,1,1).'*y;
    f1 = {A A(:,1) A(1,:) A(1,1)};
    f2 = {B B(:,1) B(1,:) B(1,1)};
    f3 = {C C(:,1,1)}; % first dimension of C must be n
    for j1 = 1:4
      for j2 = 1:4
        for j3 = 1:2
          % the 'tensor'- and 'point'-syntaxes are equivalent here
          siz = max(max(size(f1{j1}),size(f2{j2})), ...
                    [size(f3{j3},1) size(f3{j3},3)].*[0 1]);
          if siz(1) == n
            yy = clenshaw(f1{j1},f2{j2},f3{j3});
            ok = ok && norm(y-yy,inf) < 5e-15 && siz(2) == size(yy,2);
          end
        end
      end
    end
  end
end

%--------------------------------------------------------------------------
function ok = l_spin21
%L_SPIN21 (CLENSHAW) Empty cases.

ok = 1;

% 1-D syntax: empty M, N, and both
% 'recurrence'/'coefficient'-syntax
z0 = zeros(0,3); z1 = zeros(1,3);
w0 = zeros(0,1); w1 = zeros(1,1);
w2 = zeros(0,1);
cases = {z0 z1 z0 w0 w1 w0 z0 z1 z0; ...
         z1 z0 z0 z1 z0 z0 w1 w0 w0};
for AB = cases
  y  = clenshaw(AB{1},AB{2});
  ok = ok && all(size(y) == [0 3]);
  y  = clenshaw(AB{1}',AB{2}');
  ok = ok && all(size(y) == [3 0]);
end
y = clenshaw([],[]);
ok = ok && all(size(y) == [0 0]);

cases = {z0 z1 z0 w0 w1 w0 z0 z1 z0 z0 z1 z0 w0 w1 w0 z0 z1 z0; ...
         z1 z0 z0 z1 z0 z0 w1 w0 w0 z1 z0 z0 z1 z0 z0 w1 w0 w0; ...
         w2 w2 w2 w2 w2 w2 w2 w2 w2 w2 w2 w2 w2 w2 w2 w2 w2 w2};
for ABC = cases
  y  = clenshaw(ABC{1},ABC{2},ABC{3});
  ok = ok && all(size(y) == [1 3]) && all(y == 0);
end

% for brevity, some cases are not tested here...
cases = {z0 z1 z0 z0 z1 z0 w0 w1 w0 w0 w1 w0 z0 z1 z0 z0 z1 z0; ...
         z1 z0 z0 z1 z0 z0 z1 z0 z0 z1 z0 z0 w1 w0 w0 w1 w0 w0; ...
         z0 z0 z0 z1 z1 z1 z0 z0 z0 z1 z1 z1 z0 z0 z0 z1 z1 z1};
for ABC = cases
  y  = clenshaw(ABC{1}',ABC{2}',permute(ABC{3}',[1 3 2]));
  ok = ok && all(size(y) == [1 0]);
  if ~ok, break; end
end
y = clenshaw([],[],[]);
ok = ok && all(size(y) == [0 0]);

%--------------------------------------------------------------------------
function ok = l_spin22
%L_SPIN22 (CLENSHAW) Multiple recursions.

ok = 1;

% just check the size
A = {ones(3,2) ones(3,1) ones(1,2)};
B = {ones(3,2) ones(3,2) ones(3,1)};
C = {ones(3,4,1) ones(3,4,2) ones(3,4,1)};

for i = 1:size(A,2)
  y = clenshaw(A{i},B{i},C{i});
  ok = ok && all(size(y) == [4 2]);
end

%--------------------------------------------------------------------------
function ok = l_spin23
%L_SPIN23 (CLENSHAWOLD) Illegal input.

ok = 1;

try
  clenshawold(1,2,3,4,5);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshawold:e1',li);
end

try
  clenshawold(1);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshawold:e1',li);
end

try
  [a,b] = clenshawold(1,2,3);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshawold:e1',li);
end

try
  clenshawold('lala',1);
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshawold:e2',li);
end

try
  clenshawold(1,1,{1});
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshawold:e2',li);
end

try
  clenshawold(ones(1,2),ones(1,3));
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshawold:e3',li);
end

try
  clenshawold(ones(2),ones(3,2),ones(3));
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshawold:e3',li);
end

try
  clenshawold(1,2,3,{1});
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshawold:e4',li);
end

try
  clenshawold(1,2,3,{1});
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshawold:e4',li);
end

try
  clenshawold(1,2,3,'lalalalalalalalaala');
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshawold:e5',li);
end

try
  clenshawold(1,2,'la');
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshawold:e5',li);
end

try
  clenshawold({1 2},{2});
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshawold:e6',li);
end

try
  clenshawold({ones(1,2) ones(1,3)},2,'point');
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshawold:e7',li);
end

try
  clenshawold({ones(1,2) ones(1,2)},2,ones(1,1,3),'point');
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshawold:e7',li);
end

try
  clenshawold(ones(4),ones(4,1),ones(3,4));
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshawold:e8',li);
end

try
  clenshawold({ones(4) ones(3)},{ones(4,1) ones(1,3)},ones(4,2));
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshawold:e8',li);
end

try
  clenshawold(1,ones(2),ones(3,2));
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshawold:e8',li);
end

try
  clenshawold({},{},single(1));
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshawold:e9',li);
end

try
  clenshawold({},{},ones(2,1));
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshawold:e10',li);
end

try
  clenshawold({1 2},{1 2},ones(3,3,3));
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshawold:e10',li);
end

try
  clenshawold({1 2},{1 2},ones(3,3,3,2),'point');
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshawold:e10',li);
end

try
  [a,b,c] = clenshawold({1 2},{1 2});
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshawold:e11',li);
end

try
  [a,b] = clenshawold(1,{1 2 3 4});
  ok = 0;
catch
  [le,li] = lasterr;
  ok = ok && strcmp('clenshawold:e11',li);
end

%--------------------------------------------------------------------------
function ok = l_spin24
%L_SPIN24 (CLENSHAWOLD) Small cases, 1-D syntax.

ok = 1;

% Legendre polynomials
for c = [1 0.54+0.37i 0.8        0.45+0.33i; ...
         1 1          0.34+0.45i 0.45+0.89i]
  x = c(1)*linspace(-1,1,17);
  for n = 1:6
    % 'recurrence'-syntax
    j = [1:n-1]';
    A = [ones(size(x)); (2*j-1)./j*x];
    B = [0; (1-j)./j];
    y1 = clenshawold(A,B);
    y = l_recurrence(A,B);
    ok = ok && norm(y(:)-y1(:),inf) < 1e-15;

    B2 = frepmat(B,size(x));
    y2 = clenshawold(A,B2);
    ok = ok && norm(y1(:)-y2(:),inf) < 1e-15;

    % cell-vectors of length 1
    ok = ok && isequal(clenshawold({A},B),y1);
    ok = ok && isequal(clenshawold(A,{B}),y1);
    ok = ok && isequal(clenshawold({A},{B2}),y2);

    % 'coefficient'-syntax
    C = c(2)*exp(n-1:-1:0)';
    y3 = clenshawold(A,B,C);
    yy = y1.'*C;
    ok = ok && norm(y3-yy,inf) < 1e-12;

    C2 = frepmat(C,size(x));
    y4 = clenshawold(A,B,C2);
    ok = ok && norm(y4-yy,inf) < 1e-12;

    y5 = clenshawold(A,B2,C);
    ok = ok && norm(y5-yy,inf) < 1e-12;

    y6 = clenshawold(A,B2,C2);
    ok = ok && norm(y6-yy,inf) < 1e-12;

    % cell-vectors of length 1
    ok = ok && isequal(clenshawold({A},B,C),y3);
    ok = ok && isequal(clenshawold(A,{B},C),y3);
    ok = ok && isequal(clenshawold({A},{B},C2,'point'),y4);
    ok = ok && isequal(clenshawold({A},B2,C),y5);
    ok = ok && isequal(clenshawold({A},B2,C2,'point'),y6);
  end
end

% from help, just check so that it doesn't crash...
% Fibonacci numbers
n = 10; clenshawold(1,[0; 0; ones(n-2,1)])';
 
% Legendre polynomials
x = linspace(-1,1); n = 4; j = [1:n]';
A = [ones(size(x)); (2*j-1)./j*x];
B = [0; (1-j)./j];
y = clenshawold(A,B);
 
% Chebyshev (tensor-) polynomials in 2 dimensions
x = linspace(-1,1,50); y = linspace(-1,1,40); n = 4;
A1 = [ones(size(x)); x; frepmat(2*x,n-1)];
A2 = [ones(size(y)); y; frepmat(2*y,n-1)];
B = [0; 0; -ones(n-1,1)];
% plot the quadratic times the cubic polynomial only
C = zeros(n+1); C(3,4) = 1;
Z = clenshawold({A1 A2},B,C);

%--------------------------------------------------------------------------
function ok = l_spin25
%L_SPIN25 (CLENSHAWOLD) Random recurrence (1-D), all possibillities.

ok = 1;
rand('state',1429);

for c = [1 0.5+0.3i 1        0.5+0.3i; ...
         1 1        0.3+0.4i 0.4+0.8i];
  for n = 1:6
    A = frepmat(c(1)*rand,[n 14]);
    B = frepmat(c(2)*rand,[n 14]);

    % 'recurrence'-syntax
    y = l_recurrence(A,B);
    f1 = {A A(:,1) A(1,:) A(1,1)};
    f2 = {B B(:,1) B(1,:) B(1,1)};
    for j1 = 1:4
      for j2 = 1:4
        % some cases are not interesting:
        siz = max(size(f1{j1}),size(f2{j2}));
        if all(siz(1) == [n 14])
          yy = clenshawold(f1{j1},f2{j2});
          ok = ok && norm(y-yy,inf) < 5e-14 && all(siz == size(yy));
          ok = ok && isequal(clenshawold(f1{j1},f2(j2)),yy);
        end
      end
    end
  end
end

% 'coefficient'-syntax
for c = [1 0.5+0.3i 1        1        0.5+0.3i 0.3+0.6i 1        0.4+0.4i; ...
         1 1        0.3+0.4i 1        0.4+0.8i 1        0.4+0.5i 0.2+0.7i; ...
         1 1        1        0.8+0.4i 1        0.2+0.4i 0.3+0.7i 0.3+0.9i]
  for n = 1:6
    A = frepmat(c(1)*rand,[n 14]);
    B = frepmat(c(2)*rand,[n 14]);
    C = frepmat(c(3)*rand,[n 14]);
    
    y = l_recurrence(A,B);
    y = y.'*C(:,1);
    f1 = {A A(:,1) A(1,:) A(1,1)};
    f2 = {B B(:,1) B(1,:) B(1,1)};
    f3 = {C C(:,1)}; % first dimension of C must be n
    for j1 = 1:4
      for j2 = 1:4
        for j3 = 1:2
          % the 'tensor'- and 'point'-syntaxes are equivalent here
          siz = max(max(size(f1{j1}),size(f2{j2})),size(f3{j3}).*[0 1]);
          if siz(1) == n
            yy = clenshawold(f1{j1},f2{j2},f3{j3});
            ok = ok && norm(y-yy,inf) < 5e-15 && siz(2) == size(yy,1);
            ok = ok && isequal(clenshawold(f1(j1),f2(j2),f3{j3},'point'),yy);
          end
        end
      end
    end
  end
end

%--------------------------------------------------------------------------
function ok = l_spin26
%L_SPIN26 (CLENSHAWOLD) Small cases in 2-D.

ok = 1;

% Chebyshev polynomials in 2-D
for c = [1 0.54+0.37i 0.8        0.45+0.33i; ...
         1 1          0.34+0.45i 0.45+0.89i]
  x = c(1)*linspace(-1,1,5); y = c(2)*linspace(-1,1,7);
  for n = 1:3
    for m = 1:3
      A1 = [ones(size(x)); x; repmat(2*x,[n-1 1])];
      A2 = [ones(size(y)); y; repmat(2*y,[m-1 1])];
      B1 = [0; 0; -ones(n-1,1)];
      B2 = [0; 0; -ones(m-1,1)];
      C = c(1)*exp(0:-1:-n)'*(1./sqrt(1:m+1));

      % 'tensor'-evaluation
      z1 = clenshawold({A1 A2},{B1 B2});
      z = l_recurrence({A1 A2},{B1 B2});
      ok = ok && norm(z(:)-z1(:),inf) < 5e-14 && all(size(z) == size(z1));

      % multiple recursions
      w = cell(1,2);
      [w{:}] = clenshawold({A1 A2},{B1 B2});
      w = permute(reshape(w{1}(:)*w{2}(:).',[n+1 5 m+1 7]),[1 3 2 4]);
      ok = ok && norm(z(:)-w(:),inf) < 5e-14 && all(size(z) == size(w));

      z2 = clenshawold({A1 A2},{B1 B2},C);
      zz = reshape(C(:).'*reshape(z,[],5*7),5,7);
      ok = ok && norm(zz(:)-z2(:),inf) < 5e-14 && all(size(zz) == size(z2));

      B2 = frepmat(B2,size(y));
      z3 = clenshawold({A1 A2},{B1 B2},'tensor');
      ok = ok && norm(z(:)-z3(:),inf) < 5e-14 && all(size(z) == size(z3));

      z4 = clenshawold({A1 A2},{B1 B2},C);
      ok = ok && norm(zz(:)-z4(:),inf) < 5e-14 && all(size(zz) == size(z4));

      % 'point'-evaluation
      A2 = A2(:,1:5);
      B2 = B2(:,1);
      z1 = clenshawold({A1 A2},{B1 B2},'point');
      z = l_recurrence({A1 A2},{B1 B2},'point');
      ok = ok && norm(z(:)-z1(:),inf) < 5e-14 && all(size(z) == size(z1));

      z2 = clenshawold({A1 A2},{B1 B2},C,'point');
      zz = reshape(C(:).'*reshape(z,[],5),5,[]);
      ok = ok && norm(zz(:)-z2(:),inf) < 5e-14 && all(size(zz) == size(z2));

      B2 = frepmat(B2,size(x));
      z3 = clenshawold({A1 A2},{B1 B2},'point');
      ok = ok && norm(z(:)-z3(:),inf) < 5e-14 && all(size(z) == size(z3));

      z4 = clenshawold({A1 A2},{B1 B2},C,'point');
      ok = ok && norm(zz(:)-z4(:),inf) < 5e-14 && all(size(zz) == size(z4));

      C = frepmat(C,[1 1 5]);
      z5 = clenshawold({A1 A2},{B1 B2},C,'point');
      ok = ok && norm(zz(:)-z5(:),inf) < 5e-14 && all(size(zz) == size(z5));
    end
  end
end

%--------------------------------------------------------------------------
function ok = l_spin27
%L_SPIN27 (CLENSHAWOLD) Random recurrences in 3-D. Tensor/point-evaluation.

ok = 1;
rand('state',429);

% 'recurrence'-syntax
for c = [1 0.5+0.3i 1        0.5+0.3i; ...
         1 1        0.3+0.4i 0.4+0.8i];
  for n = 1:4
    A = frepmat({frepmat(rand,[n 3])},[1 3]); A{2} = c(1)*A{2};
    B = frepmat({frepmat(rand,[n 3])},[1 3]); B{3} = c(2)*B{3};

    y1 = l_recurrence(A,B);
    y2 = l_recurrence(A,B,'point');
    for j = 1:2
      f1 = {A{j}(:,1) A{j}(1,:) A{j}(1,1)};
      f2 = {B{j+1}(:,1) B{j+1}(1,:) B{j+1}(1,1)};
      for j1 = f1
        for j2 = f2
          AA = [A(1:j-1) j1 A(j+1:end)];
          BB = [B(1:j) j2 B(j+2:end)];
          yy = clenshawold(AA,BB,'tensor');
          ok = ok && norm(y1(:)-yy(:),inf) < 5e-14 && ...
               all(size(y1) == size(yy));
          w = cell(1,3);
          [w{:}] = clenshawold(AA,BB);
          ww = w{1}(:)*w{2}(:).';
          ww = ww(:)*w{3}(:).';
          ww = permute(reshape(ww,[n 3 n 3 n 3]),[1 3 5 2 4 6]);
          ok = ok && norm(y1(:)-ww(:),inf) < 5e-14 && ...
               all(size(y1) == size(ww));
          yy = clenshawold(AA,BB,'point');
          ok = ok && norm(y2(:)-yy(:),inf) < 5e-14 && ...
               all(size(y2) == size(yy));
          AA{2} = AA{2}(:,1); BB{2} = BB{2}(:,1);
          yy = clenshawold(AA,BB,'point');
          ok = ok && norm(y2(:)-yy(:),inf) < 5e-14 && ...
               all(size(y2) == size(yy));
          AA{3} = AA{3}(:,1); BB{3} = BB{3}(:,1);
          yy = clenshawold(AA,BB,'point');
          ok = ok && norm(y2(:)-yy(:),inf) < 5e-14 && ...
               all(size(y2) == size(yy));
        end
      end
    end
  end
end

% 'coefficient'-syntax
for c = [1 0.5+0.3i 1        1        0.5+0.3i 0.3+0.6i 1        0.4+0.4i; ...
         1 1        0.3+0.4i 1        0.4+0.8i 1        0.4+0.5i 0.2+0.7i; ...
         1 1        1        0.8+0.4i 1        0.2+0.4i 0.3+0.7i 0.3+0.9i]
  for n = 1:4
    A = frepmat({frepmat(rand,[n 3])},[1 3]); A{2} = c(1)*A{2};
    B = frepmat({frepmat(rand,[n 3])},[1 3]); B{3} = c(2)*B{3};
    C = frepmat(c(3)*rand,[n n n]);
    C2 = frepmat(C,[1 1 1 3]);

    y1 = l_recurrence(A,B);
    y1 = reshape(C(:).'*reshape(y1,[],3*3*3),3,3,3);
    y2 = l_recurrence(A,B,'point');
    y2 = reshape(C(:).'*reshape(y2,[],3),3,[]);
    for j = 1:2
      f1 = {A{j}(:,1) A{j}(1,:) A{j}(1,1)};
      f2 = {B{j+1}(:,1) B{j+1}(1,:) B{j+1}(1,1)};
      for j1 = f1
        for j2 = f2
          AA = [A(1:j-1) j1 A(j+1:end)];
          BB = [B(1:j) j2 B(j+2:end)];
          yy = clenshawold(AA,BB,C,'tensor');
          ok = ok && norm(y1(:)-yy(:),inf) < 5e-13 && ...
               all(size(y1) == size(yy));
          yy = clenshawold(AA,BB,C,'point');
          ok = ok && norm(y2(:)-yy(:),inf) < 5e-13 && ...
               all(size(y2) == size(yy));
          yy = clenshawold(AA,BB,C2,'point');
          ok = ok && norm(y2(:)-yy(:),inf) < 5e-13 && ...
               all(size(y2) == size(yy));
          AA{2} = AA{2}(:,1); BB{2} = BB{2}(:,1);
          yy = clenshawold(AA,BB,C,'point');
          ok = ok && norm(y2(:)-yy(:),inf) < 5e-13 && ...
               all(size(y2) == size(yy));
          yy = clenshawold(AA,BB,C2,'point');
          ok = ok && norm(y2(:)-yy(:),inf) < 5e-13 && ...
               all(size(y2) == size(yy));
          AA{3} = AA{3}(:,1); BB{3} = BB{3}(:,1);
          yy = clenshawold(AA,BB,C,'point');
          ok = ok && norm(y2(:)-yy(:),inf) < 5e-13 && ...
               all(size(y2) == size(yy));
          yy = clenshawold(AA,BB,C2,'point');
          ok = ok && norm(y2(:)-yy(:),inf) < 5e-13 && ...
               all(size(y2) == size(yy));
        end
      end
    end
  end
end

%--------------------------------------------------------------------------
function ok = l_spin28
%L_SPIN28 (CLENSHAWOLD) Empty cases.

ok = 1;

% 1-D syntax: empty M, N, and both
% 'recurrence'/'coefficient'-syntax
z0 = zeros(0,3); z1 = zeros(1,3);
w0 = zeros(0,1); w1 = zeros(1,1);
cases = {z0 z1 z0 w0 w1 w0 z0 z1 z0; ...
         z1 z0 z0 z1 z0 z0 w1 w0 w0};
for AB = cases
  y  = clenshawold(AB{1},AB{2});
  ok = ok && all(size(y) == [0 3]);
  y  = clenshawold(AB{1}',AB{2}');
  ok = ok && all(size(y) == [3 0]);
  % N-D syntax at the same time
  y  = clenshawold(AB(1),AB{2});
  ok = ok && all(size(y) == [0 3]);
  y  = clenshawold({AB{1}'},{AB{2}'});
  ok = ok && all(size(y) == [3 0]);
  y  = clenshawold(AB(1),AB{2},'point');
  ok = ok && all(size(y) == [0 3]);
  y  = clenshawold({AB{1}'},{AB{2}'},'point');
  ok = ok && all(size(y) == [3 0]);
end
y = clenshawold([],[]);
ok = ok && all(size(y) == [0 0]);
y = clenshawold({[]},{[]});
ok = ok && all(size(y) == [0 0]);

cases = {z0 z1 z0 w0 w1 w0 z0 z1 z0 z0 z1 z0 w0 w1 w0 z0 z1 z0; ...
         z1 z0 z0 z1 z0 z0 w1 w0 w0 z1 z0 z0 z1 z0 z0 w1 w0 w0; ...
         z0 z0 z0 z0 z0 z0 z0 z0 z0 w0 w0 w0 w0 w0 w0 w0 w0 w0};
for ABC = cases
  y  = clenshawold(ABC{1},ABC{2},ABC{3});
  ok = ok && all(size(y) == [3 1]) && all(y == 0);
  y  = clenshawold(ABC{1},ABC(2),ABC{3},'point');
  ok = ok && all(size(y) == [3 1]) && all(y == 0);
end
% for brevity, some cases are not tested here...
cases = {z0 z1 z0 z0 z1 z0 w0 w1 w0 w0 w1 w0 z0 z1 z0 z0 z1 z0; ...
         z1 z0 z0 z1 z0 z0 z1 z0 z0 z1 z0 z0 w1 w0 w0 w1 w0 w0; ...
         z0 z0 z0 z1 z1 z1 z0 z0 z0 z1 z1 z1 z0 z0 z0 z1 z1 z1};
for ABC = cases
  y  = clenshawold(ABC{1}',ABC{2}',ABC{3}');
  ok = ok && all(size(y) == [0 1]);
  y  = clenshawold({ABC{1}'},{ABC{2}'},ABC{3}','point');
  ok = ok && all(size(y) == [0 1]);
end
y = clenshawold([],[],[]);
ok = ok && all(size(y) == [0 1]);
y = clenshawold({[]},[],zeros(0,1));
ok = ok && all(size(y) == [0 1]);
y = clenshawold({[]},[],zeros(0,0),'point');
ok = ok && all(size(y) == [0 1]);
y = clenshawold({[]},[],zeros(0,1),'point');
ok = ok && all(size(y) == [0 1]);

% N-D syntax: empty cells, a few  cases of empty M, N and both
% 'recurrence'/'coefficient' + 'tensor'/'point'
cases = {{}      {} ones(3); ...
         ones(3) {} {}};
for s = {'tensor' 'point'}
  for AB = cases
    y = clenshawold(AB{1},AB{2},s{1});
    ok = ok && all(size(y) == 1);
    y = clenshawold(AB{1},AB{2},17,s{1});
    ok = ok && all(size(y) == 1) && y == 17;
  end
end

y = clenshawold({ones(0,2) ones(2,3)},1);
ok = ok && all(size(y) == [0 2 2 3]);
y = clenshawold({ones(0,2) ones(2,0)},1);
ok = ok && all(size(y) == [0 2 2 0]);
y = clenshawold({ones(0,4) ones(0,0)},1);
ok = ok && all(size(y) == [0 0 4 0]);
y = clenshawold({ones(0,4) ones(0,1)},1);
ok = ok && all(size(y) == [0 0 4]);

y = clenshawold(1,{ones(2,0) ones(3,0)});
ok = ok && all(size(y) == [2 3 0 0]);
y = clenshawold(1,{ones(2,0) ones(3,1)});
ok = ok && all(size(y) == [2 3 0]);
y = clenshawold(1,{ones(0,1) ones(3,1)});
ok = ok && all(size(y) == [0 3]);
y = clenshawold(1,{ones(0,1) ones(3,1) ones(4,1)},'point');
ok = ok && all(size(y) == [0 3 4]);
y = clenshawold(1,{ones(0,0) ones(3,0) ones(4,0)},'point');
ok = ok && all(size(y) == [0 3 4 0]);
y = clenshawold(1,{ones(0,0) ones(3,1) ones(4,0)},'point');
ok = ok && all(size(y) == [0 3 4 0]);
y = clenshawold(1,{ones(0,1) ones(3,0) ones(4,1)},'point');
ok = ok && all(size(y) == [0 3 4 0]);
y = clenshawold(1,{ones(0,0) ones(3,0) ones(1,0)},'point');
ok = ok && all(size(y) == [0 3 1 0]);
y = clenshawold(zeros(1,0),{ones(0,0) ones(0,0) ones(1,0)},'point');
ok = ok && all(size(y) == [0 0 1 0]);

y = clenshawold({ones(2,3) ones(0,2)},{ones(2,1) 1},ones(2,0));
ok = ok && all(size(y) == [3 2]);
y = clenshawold({ones(2,0) ones(0,2)},{ones(2,0) 1},ones(2,0));
ok = ok && all(size(y) == [0 2]);
y = clenshawold({ones(0,0) ones(0,2)},{ones(1,0) 1},ones(0,0));
ok = ok && all(size(y) == [0 2]);

y = clenshawold({ones(3,4) ones(0,4)},1,zeros(3,0,4),'point');
ok = ok && all(size(y) == [4 1]);
y = clenshawold({ones(3,4) ones(0,1)},1,zeros(3,0,4),'point');
ok = ok && all(size(y) == [4 1]);
y = clenshawold({ones(3,0) ones(0,0)},1,zeros(3,0,0),'point');
ok = ok && all(size(y) == [0 1]);
y = clenshawold({ones(3,0) ones(0,1)},1,zeros(3,0,0),'point');
ok = ok && all(size(y) == [0 1]);
y = clenshawold({ones(3,3) ones(2,3) ones(0,3)},1,zeros(3,2,0,3),'point');
ok = ok && all(size(y) == [3 1]);
y = clenshawold({ones(3,1) ones(2,3) ones(0,3)},1,zeros(3,2,0,3),'point');
ok = ok && all(size(y) == [3 1]);
y = clenshawold({ones(3,1) ones(2,3) ones(0,1)},1,zeros(3,2,0,1),'point');
ok = ok && all(size(y) == [3 1]);

% some empty multiple recursions
y = cell(1,3);
[y{:}] = clenshawold({ones(2,0) ones(3,1) ones(0,1)}, ...
                  {ones(1,0) ones(3,0) ones(0,0)});
ok = ok && all(cellfun('size',y,1) == [2 3 0]);
ok = ok && all(cellfun('size',y,2) == [0 0 0]);
y = cell(1,0);
[y{:}] = clenshawold({},{});
ok = ok && all(size(y) == [1 0]);
[y{:}] = clenshawold(1,{});
ok = ok && all(size(y) == [1 0]);

%--------------------------------------------------------------------------
