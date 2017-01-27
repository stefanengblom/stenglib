function errorbar_tick(h,w,xtype)
%ERRORBAR_TICK Adjust the width of errorbars.
%   ERRORBAR_TICK(H) adjusts the width of error bars with handle H.
%   Error bars width is given as a ratio of X axis length (1/80).
%  
%   ERRORBAR_TICK(H,W) adjust the width of error bars with handle H.
%   The input W is given as a ratio of X axis length (1/W). The result
%   is independent of the x-axis units. A ratio between 20 and 80 is
%   usually fine.
%  
%   ERRORBAR_TICK(H,W,'UNITS') adjusts the width of error bars with
%   handle H. The input W is given in the units of the current x-axis.
%
%   See also ERRORBAR.

% Author: Arnaud Laurent
% Creation : Jan 29th 2009
% MATLAB version: R2007a
%
% Notes: This function was created from a post on the french forum:
% http://www.developpez.net/forums/f148/environnements-developpement/matlab/
% Author: Jerome Briot (Dut) 
%   http://www.mathworks.com/matlabcentral/newsreader/author/94805
%   http://www.developpez.net/forums/u125006/dut/
% It was further modified by Arnaud Laurent and Jerome Briot.

% check number of arguments and provide missing values
error(nargchk(1,3,nargin))
if nargin == 1, w = 80; end
if nargin < 3, xtype = 'ratio'; end

% calculate width of error bars
if ~strcmpi(xtype,'units')
  % retrieve x limits from current axis
  dx = diff(get(gca,'XLim'));
  w = dx/w;
end

% retrieve data
hh = get(h,'children');
x = get(hh(2),'xdata');

% change xdata
x(4:9:end) = x(1:9:end)-w/2;
x(7:9:end) = x(1:9:end)-w/2;
x(5:9:end) = x(1:9:end)+w/2;
x(8:9:end) = x(1:9:end)+w/2;
set(hh(2),'xdata',x(:))
