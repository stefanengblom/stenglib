
% parallel
make('openmp',1,'fsparseonly',1,'fsparsetime',1)

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

pt50_50=tt./sum(tt)
pt50_50_total = sum(tt)

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

pt10_50=tt./sum(tt)
pt10_50_total = sum(tt)

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

pt50_10=tt./sum(tt)
pt50_10_total = sum(tt)

C=[[pt50_50(1:3), sum(pt50_50(4:5)), pt50_50(6)]; ...
   [pt50_10(1:3), sum(pt50_10(4:5)), pt50_10(6)];...
   [pt10_50(1:3), sum(pt10_50(4:5)), pt10_50(6)]]';

C = C*100;

% may run from here:
C = [ 10.8335,   7.9822,   8.5441;...
       1.9508,   2.3237,   2.4254;...
      11.1867,  26.5467,  28.3057;...
      37.8685,  31.0315,  30.4131;...
      37.2114,  29.0363,  29.4449];

pt50_50_total = 0.3402
pt50_10_total = 0.4642
pt10_50_total = 0.4355

pt50_50_total
pt50_10_total
pt10_50_total

h = figure, bar(C);

l{1}='Pre-processing'; 
l{2}='Part 1'; 
l{3}='Part 2'; 
l{4}='Part 3+4'; 
%l{5}='Part 4'; 
l{5}='Post-processing';
%l{7}='other';
set(gca,'xticklabel', l) 

legend(strcat('Data 1=', num2str(pt50_50_total,3),'sec'), ...
    strcat('Data 2=', num2str(pt50_10_total,3),'sec'), ...
    strcat('Data 3=', num2str(pt10_50_total,3),'sec'), ...
    'Location','NorthWest');
set(gcf,'PaperPositionMode','auto');
ylabel('% of total time');
print(h,'-depsc2','perc_parallel.eps')
