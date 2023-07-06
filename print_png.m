function im = print_png(filename, figh)
%% function print_png(filename, figh)
% adapted from SCREEN2JPEG:Generate a JPEG file of the current figure with
%   dimensions consistent with the figure's screen dimensions.
%    Sean P. McCarthy
%    Copyright (c) 1984-98 by MathWorks, Inc. All Rights Reserved
if nargin < 1
     error('Not enough input arguments!')
end
if nargin < 2
    figh = gcf;
end

oldscreenunits = get(figh,'Units');
oldpaperunits = get(figh,'PaperUnits');
oldpaperpos = get(figh,'PaperPosition');
set(figh,'Units','pixels');
scrpos = get(figh,'Position');
newpos = scrpos/100;
set(figh,'PaperUnits','inches',...
     'PaperPosition',newpos)
print('-dpng', filename, '-r100');

drawnow
set(figh,'Units',oldscreenunits,...
     'PaperUnits',oldpaperunits,...
     'PaperPosition',oldpaperpos)