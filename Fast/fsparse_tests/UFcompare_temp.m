%UFcompare Compare assembly of UF-matrices against RANSPARSE.

% S. Engblom 2015-03-08

% parallel
%make('openmp',1,'fsparseonly',1,'fsparsetime',1)

addpath('UFget');
index = UFget;

Mtry = 40
discard = ceil(0.05*Mtry); % discard these outliers
Mtry = Mtry+2*discard;

% Set 1
siz = 1e4;
nnz_row = 50;
nrep = 50;
nnz_ = nnz_row*siz;

% find UF-matrices of similar statistics
points = abs(index.nnz-nnz_)/siz+abs(index.ncols-siz)/siz;
points(index.ncols ~= index.nrows) = inf;
[points,ix] = sort(points);

% considering the above, we select
ixs = [778 1381];
uf50_50 = zeros(7,2);
uf50_50_total = zeros(1,2);

for j = 1:numel(ixs)
  % create input data
  P = UFget([index.Group{ixs(j)} '/' index.Name{ixs(j)}]);

  [iiuf,jjuf,ssuf] = find(P.A);
  ssuf = ssuf/nrep;
  iiuf = repmat(iiuf,1,nrep);
  jjuf = repmat(jjuf,1,nrep);
  ssuf = repmat(ssuf,1,nrep);
  p = randperm(numel(iiuf));
  iiuf = iiuf(p);
  jjuf = jjuf(p);
  ssuf = ssuf(p);

  % time
  t1 = zeros(Mtry,7);
  for i = 1:Mtry
    tic;
    [S2,t] = fsparse(iiuf,jjuf,ssuf,size(P.A));
    total_time = toc;
    t = [t total_time-sum(t)];
    t1(i,:) = t;
  end

  nnz(S2);

  [tt,ix] = sort(sum(t1,2));
  ix = sort(ix(1+discard:end-discard)); % (sort keeps the original order)
  t1 = t1(ix,:);
  tt = mean(t1,1);
  uf50_50(:,j) = tt./sum(tt)
  uf50_50_total(j) = sum(tt)
end

% Random Data
[ii, jj, ss, siz] = ransparse(siz,nnz_row,nrep);
t1 = zeros(Mtry,7);
for i = 1:Mtry
  tic;
  [S2,t] = fsparse(ii,jj,ss,[siz siz]);
  total_time = toc;
  t = [t total_time-sum(t)];
  t1(i,:) = t;
end
[tt,ix] = sort(sum(t1,2));
ix = sort(ix(1+discard:end-discard)); % (sort keeps the original order)
t1 = t1(ix,:);
tt = mean(t1,1);
pt50_50 = tt./sum(tt)
pt50_50_total = sum(tt)

% values obtained on HW2
uf50_50 = [0.1060,   0.0989;
           0.0185,   0.0196;
           0.1406,   0.1544;
           0.3375,   0.3411;
           0     ,   0.0000;
           0.3783,   0.3681;
           0.0190,   0.0178];
uf50_50_total = [0.2329 0.2410];

pt50_50 = [0.1040;    
           0.0185;    
           0.1416;    
           0.3391;    
           0.0000;  
           0.3777;    
           0.0192];
pt50_50_total = 0.2371;

nnz_50_50 = [index.nnz(ixs) nnz_];
total_time_50_50_per_nnz = [uf50_50_total pt50_50_total]./nnz_50_50

figure(1), h = bar(100*[uf50_50([1:4 6],:) pt50_50([1:4 6],:)]); drawnow
l = cell(1,5);
l{1} = 'Pre-processing'; 
l{2} = 'Part 1'; 
l{3} = 'Part 2'; 
l{4} = 'Part 3+4'; 
l{5} = 'Post-processing';
set(gca,'xticklabel', l) 

h = legend(strcat(index.Name{ixs(1)},', time/nnz = ', ...
                  num2str(total_time_50_50_per_nnz(1),3),'sec'), ...
           strcat(index.Name{ixs(2)},', time/nnz = ', ...
                  num2str(total_time_50_50_per_nnz(2),3),'sec'), ...
           strcat('Data 1 (ransparse)',', time/nnz = ', ...
                  num2str(total_time_50_50_per_nnz(3),3),'sec'), ...
           'location','NorthWest');
set(h,'fontsize',12,'interpreter','none');
%print(h,'-deps2','UF1.eps')

% Set 2
siz = 5e4;
nnz_row = 50;
nrep = 10;
nnz_ = nnz_row*siz;

% find UF-matrices
points = abs(index.nnz-nnz_)/siz+abs(index.ncols-siz)/siz;
points(index.ncols ~= index.nrows) = inf;
[points,ix] = sort(points);

%ixs = [1311 957]; % selected
ixs = [1286 852];
uf50_10 = zeros(7,2);
uf50_10_total = zeros(1,2);

for j = 1:numel(ixs)
  % create input data
  P = UFget([index.Group{ixs(j)} '/' index.Name{ixs(j)}]);

  [iiuf,jjuf,ssuf] = find(P.A);
  ssuf = ssuf/nrep;
  iiuf = repmat(iiuf,1,nrep);
  jjuf = repmat(jjuf,1,nrep);
  ssuf = repmat(ssuf,1,nrep);
  p = randperm(numel(iiuf));
  iiuf = iiuf(p);
  jjuf = jjuf(p);
  ssuf = ssuf(p);

  % time
  t1 = zeros(Mtry,7);
  for i = 1:Mtry
    tic;
    [S2,t] = fsparse(iiuf,jjuf,ssuf,size(P.A));
    total_time = toc;
    t = [t total_time-sum(t)];
    t1(i,:) = t;
  end

  nnz(S2);

  [tt,ix] = sort(sum(t1,2));
  ix = sort(ix(1+discard:end-discard)); % (sort keeps the original order)
  t1 = t1(ix,:);
  tt = mean(t1,1);
  uf50_10(:,j) = tt./sum(tt)
  uf50_10_total(j) = sum(tt)
end

% Random Data
[ii, jj, ss, siz] = ransparse(siz,nnz_row,nrep);
t1 = zeros(Mtry,7);
for i = 1:Mtry
  tic;
  [S2,t] = fsparse(ii,jj,ss,[siz siz]);
  total_time = toc;
  t = [t total_time-sum(t)];
  t1(i,:) = t;
end
[tt,ix] = sort(sum(t1,2));
ix = sort(ix(1+discard:end-discard)); % (sort keeps the original order)
t1 = t1(ix,:);
tt = mean(t1,1);
pt50_10 = tt./sum(tt)
pt50_10_total = sum(tt)

uf50_10 = [0.0914,    0.0792;
           0.0181,    0.0203;
           0.1526,    0.2558;
           0.3271,    0.2933;
           0.0000,    0.0000;
           0.3585,    0.3088;
           0.0523,    0.0427]'

uf50_10_total = [0.2561; 0.3091]'

pt50_10 = [0.0727;    0.0203;    0.2877;    0.2766;         0;  0.3033;    0.0394]'

pt50_10_total = 0.3357          

nnz_50_10 = [2454957; 2571768; siz * nnz_row]'

total_time_50_10_per_nnz = [uf50_10_total'; pt50_10_total']./nnz_50_10'

figure(2), bar([uf50_10' pt50_10']); drawnow

% Set 3
siz = 5e4;
nnz_row = 10;
nrep = 50;
nnz_ = nnz_row*siz;

% find UF-matrices
points = abs(index.nnz-nnz_)/siz+abs(index.ncols-siz)/siz;
points(index.ncols ~= index.nrows) = inf;
[points,ix] = sort(points);

%ixs = [356 1287]; % selected
ixs = [1311 955];
uf10_50 = zeros(7,2);
uf10_50_total = zeros(1,2);

for j = 1:numel(ixs)
  % create input data
  P = UFget([index.Group{ixs(j)} '/' index.Name{ixs(j)}]);
  [iiuf,jjuf,ssuf] = find(P.A);
  ssuf = ssuf/nrep;
  iiuf = repmat(iiuf,1,nrep);
  jjuf = repmat(jjuf,1,nrep);
  ssuf = repmat(ssuf,1,nrep);
  p = randperm(numel(iiuf));
  iiuf = iiuf(p);
  jjuf = jjuf(p);
  ssuf = ssuf(p);

  % time
  t1 = zeros(Mtry,7);
  for i = 1:Mtry
    tic;
    [S2,t] = fsparse(iiuf,jjuf,ssuf,size(P.A));
    total_time = toc;
    t = [t total_time-sum(t)];
    t1(i,:) = t;
  end

  nnz(S2);

  [tt,ix] = sort(sum(t1,2));
  ix = sort(ix(1+discard:end-discard)); % (sort keeps the original order)
  t1 = t1(ix,:);
  tt = mean(t1,1);
  uf10_50(:,j) = tt./sum(tt)
  uf10_50_total(j) = sum(tt)
end

[ii, jj, ss, siz] = ransparse(siz,nnz_row,nrep);
t1 = zeros(Mtry,7);
for i = 1:Mtry
  tic;
  [S2,t] = fsparse(ii,jj,ss,[siz siz]);
  total_time = toc;
  t = [t total_time-sum(t)];
  t1(i,:) = t;
end
[tt,ix] = sort(sum(t1,2));
ix = sort(ix(1+discard:end-discard)); % (sort keeps the original order)
t1 = t1(ix,:);
tt = mean(t1,1);
pt10_50 = tt./sum(tt)
pt10_50_total = sum(tt)

uf10_50 = [0.0781,    0.0789;
           0.0213,    0.0226;
           0.2960,    0.3004;
           0.2882,    0.2845;
           0.0000,    0.0000;
           0.3039,    0.2979;
           0.0125,    0.0158]'
           

uf10_50_total = [0.3173;    0.2955]'


pt10_50 = [0.0781;    0.0213;    0.3151;    0.2744;         0; \
           0.2963;    0.0148]'



pt10_50_total =  0.3056


nnz_10_50 = [512084; 486129; siz * nnz_row]'

total_time_10_50_per_nnz = [uf10_50_total'; pt10_50_total']./nnz_10_50'


figure(3), bar([uf10_50' pt10_50']); drawnow
