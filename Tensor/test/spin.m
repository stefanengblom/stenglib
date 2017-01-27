function spin(ix)
%SPIN Tests for TENSOR.

% S. Engblom 2005-08-29 (Revision, tsum added)
% S. Engblom 2005-04-27 (Revision)
% S. Engblom 2005-04-10

ftests = {@l_spin1 @l_spin2 @l_spin3 ...
          @l_spin4 @l_spin5 @l_spin6 @l_spin7 ...
          @l_spin8 @l_spin9 @l_spin10 @l_spin11 @l_spin12};
stests = {'Test #1 (tsize/tndims)' 'Test #2 (tndims)' 'Test #3 (tsize)' ...
          'Test #4 (tprod)' 'Test #5 (tprod)' 'Test #6 (tprod)' ...
          'Test #7 (tprod)' ...
          'Test #8 (tsum)' 'Test #9 (tsum)' 'Test #10 (tsum)' ...
          'Test #11 (tsum)' 'Test #12 (tsum)'};
if nargin == 0, ix = 1:size(ftests,2); end
runtest('SPIN (tensor)',ftests(ix),stests(ix));

%--------------------------------------------------------------------------
function ok = l_spin1
%L_SPIN1 (TSIZE/TNDIMS) Illegal input.

ok = 1;

try
  c = tndims('foo','baa');
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tndims:e1');
end

try
  [c,d] = tndims(rand(2));
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tndims:e1');
end

try
  [c,d] = tsize(rand(2));
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tsize:e1');
end

try
  c = tsize(rand(2),'foo','baa');
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tsize:e1');
end

try
  c = tsize(rand(2),'foo');
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tsize:e2');
end

try
  c = tsize(rand(2),[1+2i 0]);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tsize:e2');
end

try
  c = tsize(rand(2),[1.1 2]);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tsize:e3');
end

try
  c = tsize(rand(2),[0 1]);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tsize:e3');
end

try
  c = tsize(rand(2),[1 NaN]);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tsize:e3');
end

%--------------------------------------------------------------------------
function ok = l_spin2
%L_SPIN2 (TNDIMS) Test of TNDIMS.

ok = 1;

n = tndims(ones(2));
ok = ok && n == 2;
n = tndims(ones(2,1));
ok = ok && n == 1;
n = tndims(1);
ok = ok && n == 0;
n = tndims([]);
ok = ok && n == 2;

n = tndims(int32(ones(2,2,2)));
ok = ok && n == 3;

n = tndims({1,2,3});
ok = ok && n == 2;
n = tndims({1; 2; 3;});
ok = ok && n == 1;
n = tndims({});
ok = ok && n == 2;

n = tndims(struct('foo','baa'));
ok = ok && n == 0;

%--------------------------------------------------------------------------
function ok = l_spin3
%L_SPIN3 (TSIZE) Test of TSIZE.

ok = 1;

s = tsize(rand(2,3));
ok = ok && all(s == [2 3]);
s = tsize(rand(3,2,2));
ok = ok && all(s == [3 2 2]);
s = tsize(rand(2,1));
ok = ok && size(s,2) == 1 && s == 2;
s = tsize(pi);
ok = ok && all(size(s) == [1 0]);

s = tsize(rand(1:4),1:3);
ok = ok && all(s == 1:3);
s = tsize(rand(4:-1:1),1:4);
ok = ok && all(s == 4:-1:1);

s = tsize(pi,1:10);
ok = ok && size(s,2) == 10 && all(s == 1);

s = tsize(sparse(1:10,1,1));
ok = ok && size(s,2) == 1 && s == 10;
s = tsize(sparse(1:10,1,1),1:4);
ok = ok && size(s,2) == 4 && all(s == [10 1 1 1]);

%--------------------------------------------------------------------------
function ok = l_spin4
%L_SPIN4 (TPROD) Illegal input.

ok = 1;

try
  c = tprod(rand(2),rand(3),[1 2],[0 2]);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tprod:e1');
end

try
  c = tprod(rand(2),rand(3),[1 2],[1.1 2]);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tprod:e1');
end

try
  c = tprod(rand(2),rand(3),[1 NaN],[1 2]);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tprod:e1');
end

try
  c = tprod(rand(2,3,2),rand(2,1,2),[1 3 2],[2 1 1]);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tprod:e2');
end

try
  c = tprod(rand(2,3,2),rand(2,1,2),[1 2 2],[2 3 1]);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tprod:e2');
end

try
  c = tprod(rand(2),rand(3),[1 -2],[1 2]);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tprod:e3');
end

try
  c = tprod(rand(2),rand(3),[1 4],[1 2]);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tprod:e4');
end

try
  c = tprod(rand(2),rand(3),[-2 1],[1 -2]);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tprod:e4');
end

try
  c = tprod(rand(2,2,2),rand(3),[1 2],[1 2]);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tprod:e5');
end

try
  c = tprod(rand(1,1,2,2,2),rand(1,1,2,2,2),[1 -1 2 3 4],[1 -1 2 3])
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tprod:e5');
end

try
  c = tprod(rand(3),rand(2,2,2),[1 2],[1 2]);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tprod:e5');
end

try
  c = tprod(rand(2),rand(2,3,2),[1 2],[1 2 3]);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tprod:e6');
end

try
  c = tprod(rand(2,3,2),rand(2),[1 -1 2],[1 -1 3]);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tprod:e6');
end

try
  c = tprod(rand(2),rand(2),[-1 -2],[-2 -5]);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tprod:e7');
end

try
  c = tprod(rand(2),rand(2),[1 7],[2 1]);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tprod:e7');
end

try
  c = tprod(rand(2),rand(2),[1 Inf],[2 1]);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tprod:e1');
end

try
  c = tprod(rand(2),sparse(rand(2)),[1 2],[2 1]);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tprod:e8');
end

try
  c = tprod(1,2);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tprod:e9');
end

try
  [c,d] = tprod(1,2,3,4);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tprod:e9');
end

try
  c = tprod(rand(2),rand(3),int32(1:2),[1 2]);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tprod:e10');
end

try
  c = tprod(rand(2),rand(3),1:2,sparse([1 2]));
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tprod:e10');
end

try
  c = tprod(rand(2),rand(3),[1+i 2+i],1:2);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tprod:e10');
end

%--------------------------------------------------------------------------
function ok = l_spin5
%L_SPIN5 (TPROD) Various small tests.

ok = 1;
rand('state',1234);

for ac = [1 1+2i]
  for bc = [1 1+3i]
    a = ac*reshape(1:24,[4 2 3]);
    b = bc*reshape(1:36,[3 4 3]);
    c = tprod(a,b,[-1 2 1],[1 -1 3]);

    ok = ok && all(size(c) == [3 2 3]);

    % from help
    A = ac*rand(3,4); B = bc*rand(2,5);
    C1 = tprod(A,B,[1 2],[3 4]);
    c1 = reshape(A(:)*B(:).',[3 4 2 5]);
    ok = ok && norm(C1(:)-c1(:),inf) <= 1e-14;

    % Kronecker product...
    C2 = reshape(tprod(A,B,[2 4],[1 3]),size(A).*size(B));
    c2 = kron(A,B);
    ok = ok && norm(C2(:)-c2(:),inf) <= 1e-14;

    A = ac*rand(3,4); B = bc*rand(4,2);
    C2 = tprod(A,B,[1 -1],[-1 2]);
    c2 = A*B;
    ok = ok && norm(C2(:)-c2(:),inf) <= 1e-14;

    A = ac*rand(3,5); B = bc*rand(3,5);
    C3 = tprod(A,B,[1 2],[1 2]);
    c3 = A.*B;
    ok = ok && norm(C3(:)-c3(:),inf) <= 1e-14;

    C4 = tprod(A,1,[2 1],[]);
    c4 = A.';
    ok = ok && norm(C4(:)-c4(:),inf) <= 1e-14;

    C5 = tprod(A,ones(size(A)),[-1 -2],[-1 -2]);
    c5 = sum(A(:));
    ok = ok && norm(C5(:)-c5(:),inf) <= 1e-14;

    A = ac*rand(10);
    C6 = tprod(A,eye(size(A)),[-1 -2],[-1 -2]);
    c6 = sum(diag(A));
    ok = ok && norm(C6(:)-c6(:),inf) <= 1e-14;
  end
end

%--------------------------------------------------------------------------
function ok = l_spin6
%L_SPIN6 (TPROD) Advanced tests.

ok = 1;
rand('state',1084);

for ac = [1 1+2i]
  for bc = [1 1+3i]
    % the general case
    A = ac*rand(2,3,4,3,2);
    B = bc*rand(3,2,2,3,2);
    C1 = tprod(A,B,[4 -2 2 -1 1],[-1 1 3 -2 4]);

    A = permute(A,[5 1 4 2 3]);
    B = permute(B,[2 5 1 4 3]);
    C2 = tprod(A,B,[1 4 -1 -2 2],[1 4 -1 -2 3]);
    ok = ok && all(size(C1) == size(C2)) && norm(C1(:)-C2(:),inf) <= 1e-14;

    A = permute(A,[5 3 4 1 2]);
    B = permute(B,[3 4 5 1 2]);
    C3 = tprod(A,B,[2 -1 -2 1 4],[-1 -2 3 1 4]);
    ok = ok && all(size(C1) == size(C3)) && norm(C1(:)-C3(:),inf) <= 1e-14;

    a = reshape(A,[size(A,1) size(A,2)*size(A,3) size(A,4)*size(A,5)]);
    b = reshape(B,[size(B,1)*size(B,2) size(B,3) size(B,4)*size(B,5)]);
    c4 = tprod(a,b,[1 -1 3],[-1 2 3]);
    C4 = reshape(c4,[size(c4,1) size(c4,2) size(A,4) size(A,5)]);
    C4 = permute(C4,[3 1 2 4]);
    ok = ok && all(size(C1) == size(C4)) && norm(C1(:)-C4(:),inf) <= 1e-14;

    c5 = zeros(size(c4));
    for i = 1:size(c5,3)
      c5(:,:,i) = a(:,:,i)*b(:,:,i);
    end
    ok = ok && all(size(c4) == size(c5)) && norm(c4(:)-c5(:),inf) <= 1e-14;

    % no negative indices
    A = ac*rand(1,3,2,1,2);
    B = bc*rand(3,4,2,1,1);
    C1 = tprod(A,B,[3 4 7 5 1],[4 6 1 2 5]);

    A = permute(A,[1 3 5 2 4]);
    B = permute(B,[4 2 3 1 5]);
    C2 = tprod(A,B,[3 7 1 4 5],[2 6 1 4 5]);
    ok = ok && all(size(C1) == size(C2)) && norm(C1(:)-C2(:),inf) <= 1e-14;

    a = reshape(A,[size(A,1)*size(A,2) size(A,3)*size(A,4)*size(A,5)]);
    b = reshape(B,[size(B,1)*size(B,2) size(B,3)*size(B,4)*size(B,5)]);
    c3 = tprod(a,b,[1 3],[2 3]);
    C3 = reshape(c3,[size(A,1) size(A,2) ...
	    size(B,1) size(B,2) ...
	    size(A,3) size(A,4) size(A,5)]);
    C3 = permute(C3,[5 3 1 6 7 4 2]);
    ok = ok && all(size(C1) == size(C3)) && norm(C1(:)-C3(:),inf) <= 1e-14;

    c4 = zeros(size(c3));
    for i = 1:size(c4,3)
      c4(:,:,i) = a(:,i)*b(:,i).';
    end
    ok = ok && all(size(c3) == size(c4)) && norm(c3(:)-c4(:),inf) <= 1e-14;

    % no private indices
    A = ac*rand(2,3,1,2,1,2);
    B = bc*rand(2,2,2,3,1,1);
    C1 = tprod(A,B,[-1 2 3 1 -3 -2],[-2 1 -1 2 -3 3]);

    A = permute(A,[1 6 5 4 2 3]);
    B = permute(B,[3 1 5 2 4 6]);
    C2 = tprod(A,B,[-1 -2 -3 1 2 3],[-1 -2 -3 1 2 3]);
    ok = ok && all(size(C1) == size(C2)) && norm(C1(:)-C2(:),inf) <= 1e-14;

    sz = size(A);
    a = reshape(A,[prod(sz(1:3)) prod(sz(4:end))]);
    sz = size(B);
    b = reshape(B,[prod(sz(1:3)) prod(sz(4:end))]);
    c3 = tprod(a,b,[-1 1],[-1 1]);
    C3 = reshape(c3,sz(4:end));
    ok = ok && all(size(C1) == size(C3)) && norm(C1(:)-C3(:),inf) <= 1e-14;

    c4 = sum(a.*b,1).';
    ok = ok && all(size(c3) == size(c4)) && norm(c3(:)-c4(:),inf) <= 1e-14;

    % no common indices
    A = ac*rand(3,2,2,1,1);
    B = bc*rand(4,2,1,2,1);
    C1 = tprod(A,B,[2 -2 -3 3 -1],[1 -3 -1 -2 4]);

    A = permute(A,[1 4 5 2 3]);
    B = permute(B,[3 4 2 1 5]);
    C2 = tprod(A,B,[2 3 -1 -2 -3],[-1 -2 -3 1 4]);
    ok = ok && all(size(C1) == size(C2)) && norm(C1(:)-C2(:),inf) <= 1e-14;

    sza = size(A);
    a = reshape(A,[prod(sza(1:2)) prod(sza(3:end))]);
    szb = size(B);
    b = reshape(B,[prod(szb(1:3)) prod(szb(4:end))]);
    c3 = a*b;
    C3 = reshape(c3,[sza([1 2]) szb(4)]);
    C3 = permute(C3,[3 1 2 4]);
    ok = ok && all(size(C1) == size(C3)) && norm(C1(:)-C3(:),inf) <= 1e-14;

    % only common indices
    A = ac*rand(1,2,3);
    B = bc*rand(3,2,1);
    C1 = tprod(A,B,[3 2 1],[1 2 3]);

    A = permute(A,[3 2 1]);
    C2 = reshape(A(:).*B(:),size(A));
    ok = ok && all(size(C1) == size(C2)) && norm(C1(:)-C2(:),inf) <= 1e-14;

    % only private indices
    A = ac*rand(2,1,3);
    B = bc*rand(4,5,1);
    C1 = tprod(A,B,[4 5 1],[3 2 6]);

    C2 = A(:)*B(:).';
    C2 = reshape(C2,[size(A) size(B)]);
    C2 = permute(C2,[3 5 4 1 2 6]);
    ok = ok && all(size(C1) == size(C2)) && norm(C1(:)-C2(:),inf) <= 1e-14;

    % only negative indices
    A = ac*rand(1,2,3);
    B = bc*rand(2,3,1);
    C1 = tprod(A,B,[-2 -3 -1],[-3 -1 -2]);

    C2 = permute(A,[3 1 2]).*permute(B,[2 3 1]);
    C2 = sum(C2(:));
    ok = ok && all(size(C1) == size(C2)) && norm(C1(:)-C2(:),inf) <= 1e-14;
   end
end

%--------------------------------------------------------------------------
function ok = l_spin7
%L_SPIN7 (TPROD) Various pitfalls.

ok = 1;
rand('state',1245);

% 0-D/0-D
a = 3;
b = 2;
c = tprod(a,b,[],[]);
ok = ok && c == 6;

% 0-D/1-D
b = rand(3,1);
c = tprod(a,b,[],[1]);
ok = ok && norm(c-a*b,inf) <= 1e-15;

% 1-D/1-D
a = rand(3,1);
c = tprod(a,b,[1],[1]);
ok = ok && norm(c-a.*b,inf) <= 1e-15;
c = tprod(a,b,[1],[2]);
ok = ok && norm(c-a*b',inf) <= 1e-15;
c = tprod(a,b,[-1],[-1]);
ok = ok && norm(c-a'*b,inf) <= 1e-15;

% 0-D/2-D
a = 3;
b = rand(3,4);
c = tprod(a,b,[],[1 2]);
ok = ok && norm(c-a*b,inf) <= 1e-15;

% 1-D/2-D
a = rand(3,1);
c = tprod(a,b,[1],[1 2]);
ok = ok && norm(c-repmat(a,[1 4]).*b,inf) <= 1e-15;
c = tprod(a,b,[-1],[-1 1]);
ok = ok && norm(c-b'*a,inf) <= 1e-15;
c = tprod(a,b,[2],[1 3]);
c1 = permute(reshape(a*b(:)',[3 3 4]),[2 1 3]);
ok = ok && all(size(c) == size(c1)) && norm(c(:)-c1(:),inf) <= 1e-15;

% empty cases
a = rand(2,0,3,1);
b = rand(0,2,1);
c = tprod(a,b,[1 2 3 4],[2 1 4]);
ok = ok && all(size(c) == [2 0 3]);
c = tprod(a,b,[-1 3 2 1],[3 -1 1]);
ok = ok && all(size(c) == [1 3 0]);
c = tprod(a,b,[2 5 1 4],[3 2 4]);
ok = ok && all(size(c) == [3 2 0 1 0]);

%--------------------------------------------------------------------------
function ok = l_spin8
%L_SPIN8 (TSUM) Illegal inputs.

ok = 1;

try
  c = tsum(1,2,3);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tsum:e1');
end

try
  [c,d] = tsum(1,2);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tsum:e1');
end

try
  c = tsum(1,{2},1,2);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tsum:e2');
end

try
  c = tsum(1,2,1+2i,2);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tsum:e3');
end

try
  c = tsum(1,2,1.1,2);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tsum:e4');
end

try
  c = tsum(1,2,-1,2);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tsum:e4');
end

try
  c = tsum(rand(3),rand(3),[1 2],[4 5]);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tsum:e5');
end

try
  c = tsum(rand(3),rand(3),[1 2],[3 Inf]);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tsum:e4');
end

try
  c = tsum(rand(3),rand(3),[2 2],1);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tsum:e6');
end

try
  c = tsum(rand(3),rand(3),[2 1],1);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tsum:e7');
end

try
  c = tsum(rand(3),rand(3,3,2),[2 1],[1 2]);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tsum:e7');
end

try
  c = tsum(rand(3),rand(2,3),[2 1],[1 2]);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tsum:e8');
end

try
  c = tsum(rand(3),rand(3,3,2),[1 2],[1 3 2]);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tsum:e8');
end

try
  c = tsum(rand(3),0);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tsum:e9');
end

try
  c = tsum(rand(3),NaN);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tsum:e9');
end

try
  c = tsum(single(rand(3)),1);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tsum:e2');
end

try
  c = tsum(rand(3),1+2i);
  ok = 0;
catch
  [le,id] = lasterr;
  ok = ok && strcmp(id,'tsum:e3');
end

%--------------------------------------------------------------------------
function ok = l_spin9
%L_SPIN9 (TSUM) Some small tests.

ok = 1;
rand('state',2387)

% from help
v = rand(3,1); w = rand(4,1);
C1 = tsum(v,w,[1],[2]);
c1 = log(exp(v)*exp(w.'));
ok = ok && norm(C1(:)-c1(:),inf) < 1e-14 && all(size(C1) == size(c1));

A = rand(4);
C2 = 0.5*tsum(A,A,[1 2],[2 1]);
c2 = 0.5*(A+A.');
ok = ok && norm(C2(:)-c2(:),inf) < 1e-14 && all(size(C2) == size(c2));

C3 = tsum(A,0,[2 1],[]);
c3 = A.';
ok = ok && norm(C3(:)-c3(:),inf) < 1e-14 && all(size(C3) == size(c3));

C4 = tsum(A,1:4,[1 2],[3 1]);
c4 = A+repmat([1:4]',[1 4]);
ok = ok && norm(C4(:)-c4(:),inf) < 1e-14 && all(size(C4) == size(c4));

B = rand(3,4,5);
C5 = tsum(B,[1 3]);
ok = ok && all(size(C5) == [1 4]);

C6 = tsum([1:100000].^-4,-2);
ok = ok && abs(C6-pi^4/90) <= eps;

%--------------------------------------------------------------------------
function ok = l_spin10
%L_SPIN10 (TSUM) Serious tests.

ok = 1;
rand('state',340912);

AA = 1+rand(3,3,3,3,3);
BB = 1+rand(3,3,3,3,3);
for c = [1 1+1i 1    1+0.5i; ...
         1 1    1+1i 1+0.5i] 
  A = c(1)*AA;
  B = c(2)*BB;
  % all common indices
  ia = [4 3 2 5 1];
  ib = [5 2 4 3 1];
  C1 = tsum(A,B,ia,ib);
  c1 = ipermute(A,ia)+ipermute(B,ib);
  ok = ok && norm(C1(:)-c1(:),inf) <= 1e-14 && all(size(C1) == size(c1));

  % three singletons
  ia = [3 5 4 2 1 6 8 7];
  ib = [1 3 2 4 5 6 7 8];
  C1 = tsum(A,B,ia,ib);
  c1 = ipermute(A,ia)+ipermute(B,ib);
  ok = ok && norm(C1(:)-c1(:),inf) <= 1e-14 && all(size(C1) == size(c1));

  % two private indices
  ia = [5 2 7 6 3];
  ib = [4 5 6 3 1];
  C2 = tsum(A,B,ia,ib);
  c2 = log(tprod(exp(A),exp(B),ia,ib)); % dirty, but simple
  ok = ok && norm(C2(:)-c2(:),inf) <= 1e-14 && all(size(C2) == size(c2));

  % all private indices
  ia = [3 6 7 2 8];
  ib = [1 9 4 5 10];
  C2 = tsum(A,B,ia,ib);
  c2 = ipermute(reshape(log(exp(A(:))*exp(B(:).')),repmat(3,1,10)),[ia ib]);
  ok = ok && norm(C2(:)-c2(:),inf) <= 1e-14 && all(size(C2) == size(c2));
end

%--------------------------------------------------------------------------
function ok = l_spin11
%L_SPIN11 (TSUM) Empty cases.

ok = 1;
rand('state',11113);

% 0-D/0-D
a = 3;
b = 2;
c = tsum(a,b,[],[]);
ok = ok && c == 5;
c = tsum(a,b,[1],[1]);
ok = ok && c == 5;
c = tsum(a,b,[1:2],[1]);
ok = ok && c == 5;
c = tsum(a,b,[1:2],[1:3]);
ok = ok && c == 5;

% 0-D/1-D
b = rand(3,1);
c = tsum(a,b,[],[1]);
ok = ok && norm(c-(a+b),inf) <= 1e-15;
c = tsum(a,b,[2],[1]);
ok = ok && norm(c-(a+b),inf) <= 1e-15;
c = tsum(a,b,[2],[1 2 3]);
ok = ok && norm(c-(a+b),inf) <= 1e-15;

% 1-D/1-D
a = rand(3,1);
c = tsum(a,b,[1],[1]);
ok = ok && norm(c-(a+b),inf) <= 1e-15;
c = tsum(a,b,[1 2],[1]);
ok = ok && norm(c-(a+b),inf) <= 1e-15;
c = tsum(a,b,[1 2 3],[1 2]);
ok = ok && norm(c-(a+b),inf) <= 1e-15;

% 0-D/2-D
a = 3;
b = rand(3,4);
c = tsum(a,b,[],[1 2]);
ok = ok && norm(c-(a+b),inf) <= 1e-15;
c = tsum(a,b,[3 4],[1 2]);
ok = ok && norm(c-(a+b),inf) <= 1e-15;

% 1-D/2-D
a = rand(3,1);
c = tsum(a,b,[1],[1 2]);
ok = ok && norm(c-(repmat(a,[1 4])+b),inf) <= 1e-15;
c = tsum(a,b,[1 3 4],[1 2]);
ok = ok && norm(c-(repmat(a,[1 4])+b),inf) <= 1e-15;
c = tsum(a,b,[3 4],[1 2]);
cc = repmat(permute(a,[2 3 1]),[3 4])+repmat(b,[1 1 3]);
ok = ok && norm(c(:)-cc(:),inf) <= 1e-15;

% empty dimensions
a = rand(2,0,3,1);
b = rand(0,2,1);
c = tsum(a,b,[1 2 3 4],[2 1 4]);
ok = ok && all(size(c) == [2 0 3]);
c = tsum(a,b,[2 5 1 4],[3 2 4]);
ok = ok && all(size(c) == [3 2 0 1 0]);

% empy sum of one array
c = tsum(a,[]);
ok = ok && all(size(c) == [2 0 3]);
c = tsum(a,[4]);
ok = ok && all(size(c) == [2 0 3]);
c = tsum(a,[-2 -4]);
ok = ok && all(size(c) == [2 1 3]);
c = tsum(a,[1]);
ok = ok && all(size(c) == [1 0 3]);
c = tsum(a,[-3 2]);
ok = ok && all(size(c) == [2 1]);
c = tsum(b,[-5]);
ok = ok && all(size(c) == [0 2]);
c = tsum(b,[3 1]);
ok = ok && all(size(c) == [1 2]);

%--------------------------------------------------------------------------
function ok = l_spin12
%L_SPIN12 (TSUM) Sum of one array.

ok = 1;
rand('state',139847);

% isequal seems to work here (this could depend on the platform)
for c = [1 0.5+0.5i]
  B = c*rand(3,4,5,6,5,4);
  C1 = tsum(B,1);
  c1 = sum(B,1);
  ok = ok && isequal(C1,c1);
  C1 = tsum(B,[1 3]);
  c1 = sum(sum(B,1),3);
  ok = ok && isequal(C1,c1);
  C1 = tsum(B,[1 3 2]);
  c1 = sum(sum(sum(B,1),3),2);
  ok = ok && isequal(C1,c1);

  C2 = tsum(B,6);
  c2 = sum(B,6);
  ok = ok && isequal(C2,c2);
  C2 = tsum(B,[6 4]);
  c2 = sum(sum(B,6),4);
  ok = ok && isequal(C2,c2);
  C2 = tsum(B,[6 4 5 3]);
  c2 = sum(sum(sum(sum(B,6),4),5),3);
  ok = ok && isequal(C2,c2);

  b = B(end:-1:1,:,end:-1:1,:,:,:);
  C3 = tsum(B,[-1 -3]);
  c3 = sum(sum(b,1),3);
  ok = ok && isequal(C3,c3);
  b = B(end:-1:1,:,end:-1:1,:,end:-1:1,end:-1:1);
  C3 = tsum(B,[-1 -3 -5 -6]);
  c3 = sum(sum(sum(sum(b,1),3),5),6);
  ok = ok && isequal(C3,c3);
end
%--------------------------------------------------------------------------
