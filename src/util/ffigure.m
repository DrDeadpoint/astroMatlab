function [f,ax] = ffigure(varargin)
xscale = 1;
yscale = 1;
if length(varargin) == 2
    xscale = varargin{1};
    yscale = varargin{2};
end
pos = [1 1 12*xscale 6*yscale];
newvarargin = {'units','inches','position',pos};

f = figure(newvarargin{:});
ax = gca;
end