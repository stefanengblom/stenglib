percserial

percparallel

% serial - parallel - speed-up
sp50_50_speedups = [t50_50(1:3)./pt50_50(1:3), ...
               sum(t50_50(4:5))./sum(pt50_50(4:5)), ...
               t50_50(6)./pt50_50(6)]*t50_50_total/pt50_50_total

sp50_10_speedups = [t50_10(1:3)./pt50_10(1:3), ...
               sum(t50_10(4:5))./sum(pt50_10(4:5)), ...
               t50_10(6)./pt50_10(6)]*t50_10_total/pt50_10_total

sp10_50_speedups = [t10_50(1:3)./pt10_50(1:3), ...
               sum(t10_50(4:5))./sum(pt10_50(4:5)), ...
               t10_50(6)./pt10_50(6)]*t10_50_total/pt10_50_total

C=[sp50_50_speedups(1:5);...
   sp50_10_speedups(1:5);...
   sp10_50_speedups(1:5)]';

% may run from here:
C = [7.5178,   5.0820,  10.8131,   4.2548,   2.6800;...
     7.4722,   7.2648,   3.7591,   4.9092,  10.4944;...
     7.4371,   7.3973,   3.7545,   4.6655,   2.6475];

pt50_50_total = 3402
pt50_10_total = 0.4642
pt10_50_total = 0.4355

t50_50_total = 1.6136
t50_10_total = 2.9525
t10_50_total = 1.7788

h = figure, bar(C');

l{1}='Pre-processing'; 
l{2}='Part 1'; 
l{3}='Part 2'; 
l{4}='Part 3+4'; 
l{5}='Post-processing';
set(gca,'xticklabel', l) 

legend(strcat('Data 1=', num2str(t50_50_total/pt50_50_total,3),'x'), ...
    strcat('Data 2=', num2str(t50_10_total/pt50_10_total,3),'x'), ...
    strcat('Data 3=', num2str(t10_50_total/pt10_50_total,3),'x'), ...
    'Location','NorthWest');
set(gcf,'PaperPositionMode','auto');
ylabel('Speedup');
print(h,'-depsc2','speedup_serial_parallel.eps')
