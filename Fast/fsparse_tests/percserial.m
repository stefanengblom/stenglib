
% serial
make('openmp',0,'fsparseonly',1,'fsparsetime',1)

Mtry = 40;
discard = ceil(0.05*Mtry); % discard these outliers
Mtry = Mtry+2*discard;

% Set 1
[ii, jj, ss, siz] = ransparse(1e4, 50, 50);
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

t50_50=tt./sum(tt)
t50_50_total = sum(tt)

% Set 2
[ii, jj, ss, siz] = ransparse(5e4, 10, 50);
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

t10_50=tt./sum(tt)
t10_50_total = sum(tt)

% Set 3
[ii, jj, ss, siz] = ransparse(5e4, 50, 10);
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

t50_10=tt./sum(tt)
t50_10_total = sum(tt)

C=[t50_50(1:6); t50_10(1:6); t10_50(1:6)]';

C = C*100;

% may run from here:
C = [17.1707,   9.3779,  15.5561;...
      2.0901,   2.6542,   4.3923;...
     25.5023,  15.6903,  26.0170;...
     32.3756,  21.5565,  30.6695;...
      1.5936,   2.3958,   4.0670;...
     21.0250,  47.9113,  19.0844];

t50_50_total = 1.6136
t50_10_total = 2.9525
t10_50_total = 1.7788

t50_50_total
t50_10_total
t10_50_total

h = figure, bar(C);

l{1}='Pre-processing'; 
l{2}='Part 1'; 
l{3}='Part 2'; 
l{4}='Part 3'; 
l{5}='Part 4'; 
l{6}='Post-processing';
%l{7}='other';
set(gca,'xticklabel', l) 

legend(strcat('Data 1=', num2str(t50_50_total,3),'sec'), ...
       strcat('Data 2=', num2str(t50_10_total,3),'sec'), ...
       strcat('Data 3=', num2str(t10_50_total,3),'sec'), ...
       'Location','NorthWest');
set(gcf,'PaperPositionMode','auto');
ylabel('% of total time');
print(h,'-depsc2','perc.eps')
