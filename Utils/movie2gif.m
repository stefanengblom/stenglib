function movie2gif(M,cdata,file,varargin)
%MOVIE2GIF Create GIF animation from MATLAB movie.
%   MOVIE2GIF works like the MATLAB function MOVIE2AVI except that a
%   GIF animation is created instead of an AVI movie. This usually
%   saves a lot of space for simpler movies.
%
%   MOVIE2GIF(M,CDATA,FILE,...) creates a GIF animation from the movie
%   M using colors in the cell-vector CDATA and writes the result to
%   FILE. The extension '.gif' will be added to FILE if it doesn't
%   already have an extension.
%
%   Usually, CDATA kan be obtained from M itself. For instance, CDATA
%   = {M(1).cdata} might work in many situations, meaning that only
%   the colors of the first frame are used.
%
%   Property     Value               Description
%   -----------------------------------------------------------------------
%   comment      String or cell-     Comment(s) to add to the image.
%                array of strings
%
%   delaytime    Scalar in [0,655]   Specifies the delay in seconds
%                                    before displaying the next image.
% 
%   loopcount    Integer in          Specifies the number of times to repeat
%                [0..65535] | Inf    the animation. If loopcount = 0,
%                                    the animation will be played
%                                    once, if loopcount = 1, the
%                                    animation will be played twice,
%                                    and so on.
%
%   For more available properties, see the GIF-section in the help-text for 
%   IMWRITE.
%
%   Example:
%     tspan = linspace(0,2*pi,30); tspan(end) = [];
%     x = linspace(-pi,pi);
%     figure, j = 0;
%     M = struct('cdata',{},'colormap',{});
%     for t = tspan
%       j = j+1;
%       y1 = sin(x+t); y2 = cos(x-t);
%       plot(x,y1,'b',x,y2,'r');
%       M(j) = getframe;
%     end
%
%     % use cdata from first two frames only
%     movie2gif(M,{M(1:2).cdata},'test.gif','delaytime',0.05,'loopcount',inf);
%
%   See also MOVIE2AVI, AVIREAD, IMWRITE, PRIVATE/WRITEGIF.

% S. Engblom 2010-09-01 (Revision, explicitly adding extension .gif)
% S. Engblom 2010-01-14 (Revision, argument cdata)
% S. Engblom 2007-05-08

% straightforward once you know how to do it...
sz = size(M(1).cdata);
gif = zeros([sz(1:2) 1 size(M,2)],'uint8');

if ~all(cellfun('isempty',{M.colormap}))
  warning('Colormap of movie ignored.');
end

% find single color map from cdata-cells provided
[foo,map] = rgb2ind(cat(1,cdata{:}),256);
if size(map,1) < 2
  warning(['Obtained a colormap with less than two colors. ' ...
           'Try to provide more colors.']);
end

% compute GIF-frames
for j = 1:size(M,2)
  gif(:,:,:,j) = rgb2ind(M(j).cdata,map);
end

% final write
if all(file ~= '.'), file = [file '.gif']; end
imwrite(gif,map,file,'gif',varargin{:});
