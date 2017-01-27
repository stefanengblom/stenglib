%SERIALMATLAB Serial performance tests of FSPARSE.

% run manually before calling this script:
% make('openmp',false,'fsparseonly',1);

Mtry = 40;
discard = ceil(0.05*Mtry); % discard these outliers
Mtry = Mtry+2*discard;

% Set 1
[ii, jj, ss, siz] = ransparse(1e4, 50, 50);
t1 = zeros(Mtry,1);
for i = 1:Mtry
  tic;
  S1 = sparse(ii,jj,ss,siz,siz); 
  total_time = toc
  t1(i) = total_time;
end
[tt,ix] = sort(sum(t1,2));
ix = sort(ix(1+discard:end-discard)); % (sort keeps the original order)
t1 = t1(ix,:);
t_matlab_50_50 = mean(t1,1)

t2 = zeros(Mtry,1);
for i = 1:Mtry
  tic;
  S2 = fsparse(ii,jj,ss,[siz siz]);
  total_time = toc
  t2(i) = total_time;
end
[tt,ix] = sort(sum(t2,2));
ix = sort(ix(1+discard:end-discard)); % (sort keeps the original order)
t2 = t2(ix,:);
t_fsparse_50_50 = mean(t2,1)

% Set 2
[ii, jj, ss, siz] = ransparse(5e4, 50, 10);
t1 = zeros(Mtry,1);
for i = 1:Mtry
  tic;
  S1 = sparse(ii,jj,ss,siz,siz); 
  total_time = toc
  t1(i) = total_time;
end
[tt,ix] = sort(sum(t1,2));
ix = sort(ix(1+discard:end-discard)); % (sort keeps the original order)
t1 = t1(ix,:);
t_matlab_50_10 = mean(t1,1)

t2 = zeros(Mtry,1);
for i = 1:Mtry
  tic;
  S2 = fsparse(ii,jj,ss,[siz siz]);
  total_time = toc
  t2(i) = total_time;
end
[tt,ix] = sort(sum(t2,2));
ix = sort(ix(1+discard:end-discard)); % (sort keeps the original order)
t2 = t2(ix,:);
t_fsparse_50_10 = mean(t2,1)

% Set 3
[ii, jj, ss, siz] = ransparse(5e4, 10, 50);
t1 = zeros(Mtry,1);
for i = 1:Mtry
  tic;
  S1 = sparse(ii,jj,ss,siz,siz); 
  total_time = toc
  t1(i) = total_time;
end
[tt,ix] = sort(sum(t1,2));
ix = sort(ix(1+discard:end-discard)); % (sort keeps the original order)
t1 = t1(ix,:);
t_matlab_10_50 = mean(t1,1)

t2 = zeros(Mtry,1);
for i = 1:Mtry
  tic;
  S2 = fsparse(ii,jj,ss,[siz siz]);
  total_time = toc
  t2(i) = total_time;
end
[tt,ix] = sort(sum(t2,2));
ix = sort(ix(1+discard:end-discard)); % (sort keeps the original order)
t2 = t2(ix,:);
t_fsparse_10_50 = mean(t2,1)

S1 = [t_matlab_50_50/t_fsparse_50_50, ...
      t_matlab_50_10/t_fsparse_50_10, ...
      t_matlab_10_50/t_fsparse_10_50]';

S1 % record these values!

% obtained values on HW1
S1 = [2.3310, 2.0044, 2.0856]';   % (updated: Matlab 8.4)
%S1 = [2.0981, 1.7168, 1.7474]';  % (old: Matlab 7.13)

% obtained values on HW2
S2  = [2.1711, 1.4970, 1.9685]'; 

h = figure, bar([S1, S2]);

l{1}='Data 1'; 
l{2}='Data 2'; 
l{3}='Data 3'; 
set(gca,'xticklabel', l) 

legend('Hardware 1', ...
       'Hardware 2', ...
    'Location','NorthEast');

print(h,'-deps2','serial.eps')
