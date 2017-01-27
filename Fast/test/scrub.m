function scrub(ix)
%SCRUB Tests for FAST.

% FREPMAT.
% The tests in this file indicate that:
%
%   - For doubles (real or complex), FREPMAT is about 15-20% faster
%   than REPMAT (test #1:4).
%
%   - For other numerical types (single precision floats, integers,
%   characters and logicals), FREPMAT is about 50% faster than REPMAT
%   (test #5).
%
%   - For sparse matrices (real, complex or logical), FREPMAT is about
%   85-95% faster than REPMAT and allocates much less memory (test
%   #6).
%
%   - For other data-types, such as cell- and structure-arrays,
%   function-handles and objects, FREPMAT is about 10% faster than
%   REPMAT (test #7).

% FSPARSE.
% The tests in this file indicate that:
%
% - In general, when the input indices and values are of the same
%   sizes, FSPARSE is about 30-70% faster than SPARSE depending on the
%   number of colliding indices, the usage of integer indices and the
%   'nosort' option (test #1,#2,#4).
%
% - For the 'assembly' syntax, FSPARSE is about
%   50-60% faster than any equivalent construction in Matlab.

% S. Engblom 2013-12-03 (l_scrub12)
% S. Engblom 2010-01-13 (fixed 'nosort' problem on mexmaci)
% S. Engblom 2005-04-27 (Revision, merging frepmat/fsparse)
% S. Engblom 2004-10-25

if nargin == 0, ix = 1:11; end
ftests = {@l_scrub1 @l_scrub2 @l_scrub3 @l_scrub4 ...
          @l_scrub5 @l_scrub6 @l_scrub7 @l_scrub8 ...
          @l_scrub9 @l_scrub10 @l_scrub11 @l_scrub12};
stests = {'Test #1 (frepmat)' 'Test #2 (frepmat)' 'Test #3 (frepmat)' ...
          'Test #4 (frepmat)' 'Test #5 (frepmat)' 'Test #6 (frepmat)' ...
          'Test #7 (frepmat)' 'Test #8 (fsparse)' 'Test #9 (fsparse)' ...
          'Test #10 (fsparse)' 'Test #11 (fsparse)' ...
          'Test #12 (fsparse/OpenMP)'};
runtest('SCRUB (fast)',ftests(ix),stests(ix));

%--------------------------------------------------------------------------
function ok = l_scrub1
%L_SCRUB1 (FREPMAT) Large tests.

ok = 1;

% large in, large out
b = 1:1000;
a1 = frepmat(b,1000);
a2 = repmat(b,1000,1);
ok = ok & isequal(a1,a2);

b = [1:1000]';
a1 = frepmat(b,1000);
a2 = repmat(b,1000,1);
ok = ok & isequal(a1,a2);

b = [1:1000]';
a1 = frepmat(b,[1 1000]);
a2 = repmat(b,1,1000);
ok = ok & isequal(a1,a2);

b = [1:250]';
a1 = frepmat(b,[100 100]);
a2 = repmat(b,[100 100]);
ok = ok & isequal(a1,a2);

b = [1:100];
a1 = frepmat(b,[50 50 100]);
a2 = repmat(b,[50 50 100]);
ok = ok & isequal(a1,a2);

b = reshape(1:900,30,30);
a1 = frepmat(b,[100 100]);
a2 = repmat(b,[100 100]);
ok = ok & isequal(a1,a2);

b = reshape(1:600,20,30);
a1 = frepmat(b,[40 30 10]);
a2 = repmat(b,[40 30 10]);
ok = ok & isequal(a1,a2);

%--------------------------------------------------------------------------
function ok = l_scrub2
%L_SCRUB2 (FREPMAT) Small cases, many repetitions and/or many dimensions.

ok = 1;

b = pi;
a1 = frepmat(b,10000);
a2 = repmat(b,10000,1);
ok = ok & isequal(a1,a2);

b = [pi -pi];
a1 = frepmat(b,10000);
a2 = repmat(b,10000,1);
ok = ok & isequal(a1,a2);

b = [pi -pi]';
a1 = frepmat(b,10000);
a2 = repmat(b,10000,1);
ok = ok & isequal(a1,a2);

b = reshape(1:48,3,2,4,2);
a1 = frepmat(b,[40 20 30]);
a2 = repmat(b,[40 20 30]);
ok = ok & isequal(a1,a2);

a1 = frepmat(b,[10 5 10 10 30]);
a2 = repmat(b,[10 5 10 10 30]);
ok = ok & isequal(a1,a2);

a1 = frepmat(b,[10 5 10 5 10 30]);
a2 = repmat(b,[10 5 10 5 10 30]);
ok = ok & isequal(a1,a2);

%--------------------------------------------------------------------------
function ok = l_scrub3
%L_SCRUB3 (FREPMAT) Large cases, few repetitions.

ok = 1;

b = 1:100000;
a1 = frepmat(b,10);
a2 = repmat(b,10,1);
ok = ok & isequal(a1,a2);

b = [1:100000]';
a1 = frepmat(b,10);
a2 = repmat(b,10,1);
ok = ok & isequal(a1,a2);

b = [1:100000]';
a1 = frepmat(b,[1 10]);
a2 = repmat(b,[1 10]);
ok = ok & isequal(a1,a2);

b = reshape(1:100000,[100 1000]);
a1 = frepmat(b,[5 5]);
a2 = repmat(b,[5 5]);
ok = ok & isequal(a1,a2);

b = reshape(1:100000,[100 10 100]);
a1 = frepmat(b,[5 5 5]);
a2 = repmat(b,[5 5 5]);
ok = ok & isequal(a1,a2);
%--------------------------------------------------------------------------
function ok = l_scrub4
%L_SCRUB4 (FREPMAT) Various complex cases.

ok = 1;

c = 3+i;

b = c*[1:1000]';
a1 = frepmat(b,[1 1000]);
a2 = repmat(b,1,1000);
ok = ok & isequal(a1,a2);

b = c*[1:500]';
a1 = frepmat(b,[70 80]);
a2 = repmat(b,[70 80]);
ok = ok & isequal(a1,a2);

b = c*[1:100];
a1 = frepmat(b,[40 50 80]);
a2 = repmat(b,[40 50 80]);
ok = ok & isequal(a1,a2);
  
b = c*[17 -29]';
a1 = frepmat(b,10000);
a2 = repmat(b,10000,1);
ok = ok & isequal(a1,a2);

b = c*reshape(1:48,3,2,4,2);
a1 = frepmat(b,[40 40 50]);
a2 = repmat(b,[40 40 50]);
ok = ok & isequal(a1,a2);
  
b = c*reshape(1:10000,[100 100]);
a1 = frepmat(b,[5 5]);
a2 = repmat(b,[5 5]);
ok = ok & isequal(a1,a2);

%--------------------------------------------------------------------------
function ok = l_scrub5
%L_SCRUB5 (FREPMAT) Various non-double cases.

ok = 1;

for ttype = {'single' 'uint8' 'int16' 'int32'}
  type = ttype{1};
  
  b = feval(type,[1:300]');
  a1 = frepmat(b,[1 300]);
  a2 = repmat(b,1,300);
  ok = ok & isequal(a1,a2);

  b = feval(type,[1:150]');
  a1 = frepmat(b,[50 80]);
  a2 = repmat(b,[50 80]);
  ok = ok & isequal(a1,a2);

  b = feval(type,[1:70]);
  a1 = frepmat(b,[40 40 80]);
  a2 = repmat(b,[40 40 80]);
  ok = ok & isequal(a1,a2);
  
  b = feval(type,[17 -29]');
  a1 = frepmat(b,10000);
  a2 = repmat(b,10000,1);
  ok = ok & isequal(a1,a2);

  b = feval(type,reshape(1:48,3,2,4,2));
  a1 = frepmat(b,[40 20 30]);
  a2 = repmat(b,[40 20 30]);
  ok = ok & isequal(a1,a2);
  
  b = feval(type,reshape(1:10000,[100 100]));
  a1 = frepmat(b,[5 5]);
  a2 = repmat(b,[5 5]);
  ok = ok & isequal(a1,a2);
end

% char
b = 'The original REPMAT is really slow sometimes... ';
a1 = frepmat(b,[1000 100]);
a2 = repmat(b,[1000 100]);
ok = ok & isequal(a1,a2);

% logical
b = 'frepmat' > ' repmat';
a1 = frepmat(b,[500 1000]);
a2 = repmat(b,[500 1000]);
ok = ok & isequal(a1,a2);

%--------------------------------------------------------------------------
function ok = l_scrub6
%L_SCRUB6 (FREPMAT) Sparse cases.

ok = 1;

% very sparse cases
s = sparse(101:110,810:10:900,1:10,1000,1000);
a1 = frepmat(s,[5 5]); % works fine: a1 = frepmat(s,[1000 1000]);
a2 = repmat(s,[5 5]);  % extremely slow!
ok = ok & isequal(a1,a2);

s = sparse(10,81:90,1:10,100,100);
a1 = frepmat(s,[100 10]);
a2 = repmat(s,[100 10]); % extremely slow!
ok = ok & isequal(a1,a2);

% dense sparse cases
rand('state',1789);
s = rand(100);
s(s < 0.1) = 0;
s = sparse(s);
a1 = frepmat(s,[50 10]);
a2 = repmat(s,[50 10]); % extremely slow!
ok = ok & isequal(a1,a2);

s = (1+2*i)/sqrt(3)*rand(20,1000);
s(abs(s) < 0.3) = 0;
s = sparse(s);
a1 = frepmat(s,[5 100]);
a2 = repmat(s,[5 100]); % extremely slow!
ok = ok & isequal(a1,a2);

% sparse logical
s = rand(200,500);
s = abs(s) < 0.1;
s = sparse(s);
a1 = frepmat(s,[15 50]);
a2 = repmat(s,[15 50]); % extremely slow!
ok = ok & isequal(a1,a2);

%--------------------------------------------------------------------------
function ok = l_scrub7
%L_SCRUB7 (FREPMAT) Cell- and structure-arrays.

ok = 1;

b = cell(1,1000); b{3} = 10;
a1 = frepmat(b,1000);
a2 = repmat(b,1000,1);
ok = ok & isequal(a1,a2);

b = cell(1000,1); b{10} = pi;
a1 = frepmat(b,1000);
a2 = repmat(b,1000,1);
ok = ok & isequal(a1,a2);

b = cell(10,10); b{3,7} = 'nan'; b{1,5} = -pi;
a1 = frepmat(b,[20 20 50]);
a2 = repmat(b,[20 20 50]);
ok = ok & isequal(a1,a2);

b = cell(5,6,5,6); b{2,3,4,2} = ones(10);
a1 = frepmat(b,[40 50 1 20]);
a2 = repmat(b,[40 50 1 20]);
ok = ok & isequal(a1,a2);

b = struct('foo',1:10,'bar',ones(50),'foobar',zeros(10,10,10));
a1 = frepmat(b,[15 17 18 20]);
a2 = repmat(b,[15 17 18 20]);
ok = ok & isequal(a1,a2);

%--------------------------------------------------------------------------
function ok = l_scrub8
%L_SCRUB8 (FSPARSE) Inputs of the same sizes.

ok = 1;

rand('state',3123);
tt = zeros(9,5);
nz = [100 500 1000 5000 10000 50000 100000 500000 1000000];
rep = [250 200 150 100 50 20 10 2 1];
for siz = [100 100000]
  for i = 1:9
    for j = 1:rep(i)
      ii = ceil(rand(nz(i),1)*siz);
      jj = ceil(rand(nz(i),1)*siz);
      ss = rand(nz(i),1);

      t = cputime;
      S = sparse(ii,jj,ss,siz,siz);
      t = cputime-t;
      tt(i,1) = tt(i,1)+t;

      t = cputime;
      S1 = fsparse(ii,jj,ss,[siz siz]);
      t = cputime-t;
      tt(i,2) = tt(i,2)+t;
    
      t = cputime;
      S2 = fsparse(ii,jj,ss,[siz siz],'nosort');
      t = cputime-t;
      tt(i,3) = tt(i,3)+t;

      ii = int32(ii);
      jj = int32(jj);

      t = cputime;
      S3 = fsparse(ii,jj,ss,[siz siz]);
      t = cputime-t;
      tt(i,4) = tt(i,4)+t;

      t = cputime;
      S4 = fsparse(ii,jj,ss,[siz siz],'nosort');
      t = cputime-t;
      tt(i,5) = tt(i,5)+t;
      ok = ok & norm(S1-S,inf) < 1e-10 & nnz(S1) == nzmax(S1) ...
           & norm((S2')'-S,inf) < 1e-10 & nnz(S2) == nzmax(S2);
      ok = ok & norm(S3-S,inf) < 1e-10 & nnz(S3) == nzmax(S3) ...
           & norm((S4')'-S,inf) < 1e-10 & nnz(S4) == nzmax(S4);
    end
  end
end

% statistics
tt = tt./repmat(rep',[1 size(tt,2)]);
ss = 1-repmat(tt(:,1),[1 size(tt,2)-1]).\tt(:,2:end)
% (ss > 0 means speedup vs. Matlab, < 0 a slowdown)

%--------------------------------------------------------------------------
function ok = l_scrub9
%L_SCRUB9 (FSPARSE) Cases with one scalar input.

ok = 1;

rand('state',3125);
tt = zeros(7,2);
nz = [100 500 1000 5000 10000 50000 100000];
rep = [250 200 150 100 50 20 10];
for siz = [1000 100000]
  for i = 1:7
    for j = 1:rep(i)
      ii = ceil(rand(nz(i),1)*siz);
      jj = ceil(rand(nz(i),1)*siz);
      ss = rand(nz(i),1);

      t = cputime;
      S = sparse(ii,jj,1,siz,siz);
      t = cputime-t;
      tt(i,1) = tt(i,1)+t;

      t = cputime;
      S1 = fsparse(ii,jj,1,[siz siz]);
      t = cputime-t;
      tt(i,2) = tt(i,2)+t;
    
      ok = ok & norm(S1-S,inf) < 1e-10 & nnz(S1) == nzmax(S1);
      
      t = cputime;
      S = sparse(ii,1,ss,siz,siz);
      t = cputime-t;
      tt(i,1) = tt(i,1)+t;

      t = cputime;
      S1 = fsparse(ii,1,ss,[siz siz]);
      t = cputime-t;
      tt(i,2) = tt(i,2)+t;
    
      ok = ok & norm(S1-S,inf) < 1e-10 & nnz(S1) == nzmax(S1);
      
      t = cputime;
      S = sparse(1,jj,ss,siz,siz);
      t = cputime-t;
      tt(i,1) = tt(i,1)+t;

      t = cputime;
      S1 = fsparse(1,jj',ss',[siz siz]);
      t = cputime-t;
      tt(i,2) = tt(i,2)+t;
    
      ok = ok & norm(S1-S,inf) < 1e-10 & nnz(S1) == nzmax(S1);
    end
  end
end

% statistics
tt = tt./repmat(rep',[1 size(tt,2)]);
ss = 1-repmat(tt(:,1),[1 size(tt,2)-1]).\tt(:,2:end)
% (ss > 0 means speedup vs. Matlab, < 0 a slowdown)

%--------------------------------------------------------------------------
function S = l_sparse(i,j,s,siz);
%L_SPARSE Simulator of the 'assembly' syntax in FSPARSE.

% no error-checking whatsoever
isiz = size(i);
jsiz = size(j);
ssiz = size(s);

j = repmat(j,[isiz(1)/jsiz(1) 1]);
i = repmat(i,[1 jsiz(2)/isiz(2)]);
s = repmat(s,[isiz(1) jsiz(2)]./ssiz);

S = sparse(i(:),j(:),s(:),siz(1),siz(2));

%--------------------------------------------------------------------------
function ok = l_scrub10
%L_SCRUB10 (FSPARSE) Assembling finite differences.

ok = 1;

% a 3-D stencile simulation
N = 70; N2 = N*N; N3 = N*N2;
ix = reshape(1:N3,[N N N]);
ix([1 N],:,:) = [];
ix(:,[1 N],:) = [];
ix(:,:,[1 N],:) = [];
ix = ix(:);
S1 = fsparse(ix,[ix-N2 ix-N ix-1 ix ix+1 ix+N ix+N2], ...
             [0.5 1 2 -4 2 1 0.5],[N3 N3]);
S2 = l_sparse(ix,[ix-N2 ix-N ix-1 ix ix+1 ix+N ix+N2], ...
              [0.5 1 2 -4 2 1 0.5],[N3 N3]);

ok = ok & isequal(S1,S2) & nnz(S1) == nzmax(S1);

% same idea but using a bigger (and complex) molecule
mol = [-0.25 -0.25 0.5 -0.5 -0.5 1 2 -4 ...
       2 1 -0.5 -0.5 0.5 -0.25 -0.25]*(1+2i);
S1 = fsparse(ix,[ix-N2-1 ix-N2+1 ix-N2 ix-N-1 ix-N+1 ix-N ix-1 ix ...
                 ix+1 ix+N ix+N+1 ix+N-1 ix+N2 ix+N2-1 ix+N2-1], ...
             mol,[N3 N3]);
S2 = l_sparse(ix,[ix-N2-1 ix-N2+1 ix-N2 ix-N-1 ix-N+1 ix-N ix-1 ix ...
                  ix+1 ix+N ix+N+1 ix+N-1 ix+N2 ix+N2-1 ix+N2-1], ...
              mol,[N3 N3]);

ok = ok & isequal(S1,S2) & nnz(S1) == nzmax(S1);

%--------------------------------------------------------------------------
function ok = l_scrub11
%L_SCRUB11 (FSPARSE) Many zeros to be squeezed out.

ok = 1;

rand('state',3123);
tt = zeros(10,2);
nz = [100000:100000:1000000];
siz = 100000;
for i = 1:10
  ii = ceil(rand(nz(i),1)*siz);
  jj = ceil(rand(nz(i),1)*siz);
  ss = rand(nz(i),1);
  ss(ss < 0.3) = 0;

  t = cputime;
  S = sparse(ii,jj,ss,siz,siz);
  t = cputime-t;
  tt(i,1) = tt(i,1)+t;

  t = cputime;
  S1 = fsparse(ii,jj,ss,[siz siz]);
  t = cputime-t;
  tt(i,2) = tt(i,2)+t;
  
  ok = ok & norm(S1-S,inf) < 1e-10 & nnz(S1) == nnz(S);

  ss = complex(rand(nz(i),1),rand(nz(i),1));
  iz = imag(ss) < 0.3;
  ss(iz) = ss(iz)-1i*imag(ss(iz));
  iz = real(ss) > 0.7;
  ss(iz) = ss(iz)-real(ss(iz));

  t = cputime;
  S = sparse(ii,jj,ss,siz,siz);
  t = cputime-t;
  tt(i,1) = tt(i,1)+t;

  t = cputime;
  S1 = fsparse(ii,jj,ss,[siz siz]);
  t = cputime-t;
  tt(i,2) = tt(i,2)+t;
  
  ok = ok & norm(S1-S,inf) < 1e-10 & nnz(S1) == nnz(S);
end

% statistics
tt;
ss = 1-repmat(tt(:,1),[1 size(tt,2)-1]).\tt(:,2:end)
% (ss > 0 means speedup vs. Matlab, < 0 a slowdown)

%--------------------------------------------------------------------------
function ok = l_scrub12
%L_SCRUB12 (FSPARSE/OpenMP)

ok = 1;
for nthreads = 1:16
  % dummy-call, the value of nthreads is remembered afterwards
  fsparse(1,1,1,[],[],nthreads);
  ok = ok && spin(9:12);
  ok = ok && spin(8); % there is an nthreads-call in this one...
end

%--------------------------------------------------------------------------

return;

% small performance test
rand('state',123); % http://www.mathworks.se/help/matlab/math/updating-your-random-number-generator-syntax.html
siz = 1e4;         % size
nnz_row = 50;      % number of nonzeros per row
nnz = nnz_row*siz; % number of nonzeros
nrep = 50;         % number of collisions per row

ii = frepmat((1:siz)',[1 nnz_row]);
jj = ceil(rand(siz,nnz_row)*siz);
ii = frepmat(ii(:),[1 nrep]);
jj = frepmat(jj(:),[1 nrep]); % (some jj's might be the same)
p = randperm(numel(ii));
ii = ii(p);
jj = jj(p);
ss = ones(size(ii));

clear functions; % uses nthreads = omp_get_max_threads()

% to get some feeling for it...
for i = 1:3
  tic, S1 = sparse(ii,jj,ss,siz,siz); toc
end
make('openmp',false,'fsparseonly',1);
for i = 1:3
  tic, S2 = fsparse(ii,jj,ss,[siz siz]); toc
end
make('openmp',true,'fsparseonly',1);
for i = 1:3
  tic, S2 = fsparse(ii,jj,ss,[siz siz]); toc
end

% CAUTION: don't cut-and-paste this code many times (crashes due to
% some Matlab problems when recompiling mex-files)

Mtry = 40;
discard = ceil(0.05*Mtry); % discard these outliers
Mtry = Mtry+2*discard;

% serial
make('openmp',0,'fsparseonly',1,'fsparsetime',1)
t1 = zeros(Mtry,7);
for i = 1:Mtry
  tic;
  [S2,t] = fsparse(ii,jj,ss,[siz siz]);
  total_time = toc
  t = [t total_time-sum(t)];
  t1(i,:) = t;
end
[tt,ix] = sort(sum(t1,2));
ix = sort(ix(1+discard:end-discard)); % (sort keeps the original order)
t1 = t1(ix,:);
tt = mean(t1,1);
figure, pie(tt/sum(tt),{'getix()' 'Part 1' 'Part 2' 'Part 3' 'Part 4' ...
                    'sparse\_insert()' 'Other'});

% parallel
make('openmp',1,'fsparseonly',1,'fsparsetime',1)
t2 = zeros(Mtry,7);
for i = 1:Mtry
  tic;
  [S2,t] = fsparse(ii,jj,ss,[siz siz]);
  total_time = toc
  t = [t total_time-sum(t)];
  t2(i,:) = t;
end
[tt,ix] = sort(sum(t2,2));
ix = sort(ix(1+discard:end-discard));
t2 = t2(ix,:);
tt = mean(t2,1);
figure, pie(tt([1:4 6:7])/sum(tt),{'getix()' 'Part 1' 'Part 2' 'Part 3+4' ...
                    'sparse\_insert()' 'Other'});

% (old version)
% $$$ % speedup
% $$$ figure, bar(diag(mean(t1(:,1:end-1)./t2(:,1:end-1),1),1),6);
% $$$ set(gca,'xticklabel',{})
% $$$ legend('getix()','Part 1','Part 2','Part 3','Part 4', ...
% $$$        'sparse\_insert()');
% $$$ 
% $$$ mean(t1(:,1:6)./t2(:,1:6))
% $$$ std(t1(:,1:6)./t2(:,1:6))

% when Part 3+4 are counted together (for compatibility with #define
% FSPARSE_JOINT_3_4):
figure, bar(diag(mean([t1(:,1:3) sum(t1(:,4:5),2) t1(:,6)]./ ...
                      [t2(:,1:3) sum(t2(:,4:5),2) t2(:,6)],1)),6);
set(gca,'xticklabel',{})
legend('getix()','Part 1','Part 2','Part 3+4', ...
       'sparse\_insert()');

mean([t1(:,1:3) sum(t1(:,4:5),2) t1(:,6)]./ ...
     [t2(:,1:3) sum(t2(:,4:5),2) t2(:,6)],1)
std([t1(:,1:3) sum(t1(:,4:5),2) t1(:,6)]./ ...
    [t2(:,1:3) sum(t2(:,4:5),2) t2(:,6)],1)

% (Mtry = 10,  mexmaci64, 7.14.0.739 (R2012a), -O5, 4 cores Intel)
% #undef FSPARSE_DELAYED_IRANK, #undef FSPARSE_JOINT_3_4
% 1.9426    1.4528    1.6564    1.0454    1.4140
% #undef FSPARSE_DELAYED_IRANK, #define FSPARSE_JOINT_3_4
% 1.9689    1.4564    1.6660    0.9500    1.4286
% #define FSPARSE_DELAYED_IRANK, #undef FSPARSE_JOINT_3_4 *** WINNER ***
% 1.9750    1.4135    1.6619    1.1548    1.4105
% #define FSPARSE_DELAYED_IRANK, #define FSPARSE_JOINT_3_4
% 1.9425    1.4096    1.6530    1.1044    1.3872
%
% Conclusion: DELAYED_IRANK seems good, JOINT_3_4 seems bad.
%
% #define FSPARSE_PAR_INSERT 
% #define FSPARSE_DELAYED_IRANK
% #define FSPARSE_JOINT_3_4
% 1.8037    1.4272    1.5936    1.5702    1.2364
%
% Conclusion: poorer scaling on sparse_insert(), but a gain in Part 3+4.
%
% Discarding outliers implemented.
% 
% #undef FSPARSE_SECTIONS
% 1.8299    1.4421    1.5985    1.5560    1.2016
% #define FSPARSE_SECTIONS
% 1.6158    1.4270    1.5995    1.5660    1.0546
% #define FSPARSE_SECTIONS (new version)
% 1.8222    1.4155    1.6239    1.6125    1.1799
% Conclusion: nothing interesting, but I now know what to do!

% (Mtry = 10,  mexa64, 7.13.0.564 (R2011b), -O5, 6 cores Intel):
% #undef FSPARSE_DELAYED_IRANK, #undef FSPARSE_JOINT_3_4
% 3.9902    3.5824    4.6128    1.2460    2.2075
% #undef FSPARSE_DELAYED_IRANK, #define FSPARSE_JOINT_3_4
% 4.1790    3.3879    4.6482    1.1888    2.2196
% #define FSPARSE_DELAYED_IRANK, #undef FSPARSE_JOINT_3_4
% 4.0751    3.5408    4.5791    1.4968    2.2012
% #define FSPARSE_DELAYED_IRANK, #define FSPARSE_JOINT_3_4
% 4.1852    3.5728    4.5465    1.4724    2.1907
%
% Conclusion: DELAYED_IRANK seems good, JOINT_3_4 seems bad (but
% not VERY bad).
%
% #define FSPARSE_PAR_INSERT 
% #define FSPARSE_DELAYED_IRANK
% #define FSPARSE_JOINT_3_4
% 4.2497    3.5600    4.5610    2.4859    1.5150
%
% Conclusion: again a poorer scaling on sparse_insert(), but a gain in
% Part 3+4.
%
% Discarding outliers implemented.
%
% #undef FSPARSE_SECTIONS
% 4.5471    3.5315    4.4511    2.5338    1.6953
% #define FSPARSE_SECTIONS
% 2.5742    3.4451    4.4657    2.5253    1.1184
% #define FSPARSE_SECTIONS (new versions)
% 4.5390    3.6303    4.5855    2.5248    1.6940

%--------------------------------------------------------------------------

% $$$ Results on mexa64, 7.13.0.564 (R2011b), -O5:
% $$$
% $$$ S1:
% $$$ Elapsed time is 4.092237 seconds.
% $$$ Elapsed time is 3.965469 seconds.
% $$$ Elapsed time is 3.970422 seconds.
% $$$
% $$$ S2: make('openmp',false,'fsparseonly',1);
% $$$ Elapsed time is 1.984723 seconds.
% $$$ Elapsed time is 1.901968 seconds.
% $$$ Elapsed time is 1.906296 seconds.
% $$$
% $$$ S2: make('openmp',true,'fsparseonly',1);
% $$$ Elapsed time is 0.732607 seconds.
% $$$ Elapsed time is 0.725243 seconds.
% $$$ Elapsed time is 0.726692 seconds.
% $$$ (6 cores and threads)

% $$$ Results on mexmaci64, 7.14.0.739 (R2012a), -O5
% $$$ 
% $$$ S1:
% $$$ Elapsed time is 4.605108 seconds.
% $$$ Elapsed time is 4.720203 seconds.
% $$$ Elapsed time is 4.616316 seconds.
% $$$ 
% $$$ S2: make('openmp',false,'fsparseonly',1);
% $$$ Elapsed time is 2.180359 seconds.
% $$$ Elapsed time is 2.189002 seconds.
% $$$ Elapsed time is 2.144102 seconds.
% $$$ 
% $$$ S2: make('openmp',true,'fsparseonly',1);
% $$$ Elapsed time is 1.493890 seconds.
% $$$ Elapsed time is 1.513800 seconds.
% $$$ Elapsed time is 1.534640 seconds.
% $$$ (4 cores and threads)

% $$$ Results on mexmaci64, 7.14.0.739 (R2012a), -O5
% $$$ Positive speedup for all parts!
% $$$
% $$$ S1:
% $$$ Elapsed time is 6.345001 seconds.
% $$$ Elapsed time is 6.356427 seconds.
% $$$ Elapsed time is 6.346332 seconds.
% $$$
% $$$ S2: make('openmp',false,'fsparseonly',1);
% $$$ Elapsed time is 3.646333 seconds.
% $$$ Elapsed time is 3.336019 seconds.
% $$$ Elapsed time is 3.354009 seconds.
% $$$
% $$$ S2: make('openmp',true,'fsparseonly',1);
% $$$ Elapsed time is 2.714618 seconds.
% $$$ Elapsed time is 2.680104 seconds.
% $$$ Elapsed time is 2.692401 seconds.
% $$$ (2 cores and threads)