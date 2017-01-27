function spym(S,n)
%SPYM Visualize magnitude of elements in matrix.
%   SPYM(S) draws 10 filled contourlines of LOG10(ABS(S)). Since SPYM uses
%   FULL(S) this will not work for large sparse matrices.
%
%   SPYM(S,N) draws N contourlines.
%
%   Example:
%     spym(gallery('chebspec',100));
%
%   See also SPY, CONTOURF.

% S. Engblom 2006-10-31

% default number of contour lines
if nargin < 2, n = 10; end

% prepare window
cax = newplot;
next = lower(get(cax,'NextPlot'));
hold_state = ishold;

% plot
[i,j,s] = find(S);
cutoff = min(abs(s(s~= 0)));
s = log10(max(abs(s),cutoff));
S = full(fsparse(i,j,s,size(S),'nosort'));
S(S == 0) = NaN;
contourf(S,n);

% scaling
[m,n] = size(S);
xlabel(['nz = ' int2str(nnz(S))]);
set(gca,'xlim',[0 n+1],'ylim',[0 m+1],'ydir','reverse', ...
        'grid','none','plotboxaspectratio',[n+1 m+1 1]);

if ~hold_state, set(cax,'NextPlot',next); end
