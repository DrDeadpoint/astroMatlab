function [xCurves,yCurves] = ZVC(C,mu3,varargin)
% [xCurves,yCurves] = ZVC(C,mu3,varargin)
% find x and y coordinates of ZVC
% C is the Jacobi Constant you want the ZVC for
% mu3 is the 3-body mu value
% varargin:
%   (...,'plot',true)
%   this will plot the ZVC on an open figure with gray in forbidden region
%   (...,'x1x2',[x1 x2])
%   (...,'y1y2',[y1 y2])
%   if it isn't plotting big enough, change the lower and upper bounds (x1 < x2)
%   (...,'resolution',resol)
%   this allows to change the resolution of calculated x and y (default 10^-2)
%   (...,'color',col)
%   assuming you also called 'plot', this chooses the color of the
%   forbidden region
%
% [x,y] will be cell arrays, with multiple cells of coordinates based on
% how many separate ZVC sections there are
%
% Alex Hoffman
% 10/04/2020
col = [0.8 0.8 0.8];
isPlot = false;
resol = 10^-2;
x1 = -2; x2 = 2;
y1 = -2; y2 = 2;
for i = 1:length(varargin)/2
   switch varargin{i*2-1}
       case 'plot'
           isPlot = varargin{i*2};
       case 'x1x2'
           x1x2 = varargin{i*2};
           x1 = x1x2(1); x2 = x1x2(2);
       case 'y1y2'
           y1y2 = varargin{i*2};
           y1 = y1y2(1); y2 = y1y2(2);
       case 'resolution'
           resol = varargin{i*2};
       case 'color'
           col = varargin{i*2};
       otherwise
           error('variable argument not found')
   end
end
% 0 = (x^2+y^2) + 2*(1-mu3)/d + 2*mu/r - C
%% calculate coordinates of ZVC
[L1, L2, L3, L4, ~] = lagrangePoints(mu3); %nondimensional coordinates
zeroVel = [0; 0; 0];
C1 = jacobiConstant(zeroVel,L1,mu3);
C2 = jacobiConstant(zeroVel,L2,mu3);
C3 = jacobiConstant(zeroVel,L3,mu3);
C4 = jacobiConstant(zeroVel,L4,mu3);
if C > C1
    curves = 3;
    Ctype = 1;
elseif C >= C2
    curves = 2;
    Ctype = 2;
elseif C >= C3
    curves = 1; %really 1 curve but needs two sets of curves to connect at the end
    Ctype = 3;
elseif C >= C4
    curves = 2;
    Ctype = 4;
else
    Ctype = 5;
    curves = 0;
end
xCurves = cell(1,curves);
yCurves = cell(1,curves);
% sweep through x-values and find all zeros 
% start at x=-10, work my way right

xVec = linspace(x1,x2,(x2-x1)/resol);
yVec = linspace(y1,y2,(y2-y1)/resol);
xCurve1 = []; yCurve1 = []; %outside or top curve
xCurve2 = []; yCurve2 = []; %P1 or bottom curve
xCurve3 = []; yCurve3 = []; %P2 curve
xCurve4 = []; yCurve4 = []; %for when the very right bows inward a little
lastNumZeros = 0; %only used in C < C1 case
pastP1 = false;
pastP2 = false;
for i = 1:length(xVec)
    x = xVec(i);
    posMat = [ones(1,length(yVec))*x; yVec; zeros(1,length(yVec))];
    zeroVel = zeros(3,length(yVec));
    Cvec = jacobiConstant(zeroVel,posMat,mu3); %this will be equal to C along the ZVC
    Cdiff = find(diff(sign(Cvec - C))); %finds indices of sign changes
    if isempty(Cdiff) %no ZVC along this x-value
       continue  
    end
    % now use iterative solver to find exact locations of zeros
    if length(Cdiff) ~= 2 && length(Cdiff) ~= 4
        continue %might get length(yzeros) = 3 if really unlucky, that would break the code
    end
    yzeros = zeros(1,length(Cdiff)/2);
    func = @(y) (x^2+y^2) + 2*(1-mu3)/sqrt((x+mu3)^2+y^2) + 2*mu3/sqrt((x-1+mu3)^2+y^2) - C;
    dfunc = @(y) 2*y - 2*y*(1-mu3)*((x+mu3)^2+y^2)^(-3/2) - 2*y*mu3*((x-1+mu3)^2+y^2)^(-3/2);
    for j = 1:length(Cdiff)/2
        y0 = yVec(Cdiff(j));
        [yzero,~] = newtonMethod(func,dfunc,y0,10^-6,'diff');
        yzeros(j) = yzero;
    end
    yzeros = [yzeros -flip(yzeros)];
    switch Ctype
        case 1
            if ~pastP1 && lastNumZeros == 4 && length(yzeros) == 2
                pastP1 = true; %checks if past body P1, so can assign to Curve3
            elseif pastP1 && lastNumZeros == 4 && length(yzeros) == 2
                pastP2 = true;
            end
            lastNumZeros = length(yzeros);
            if length(yzeros) == 2
                xCurve1 = [x xCurve1 x]; %outside curve
                yCurve1 = [yzeros(1) yCurve1 yzeros(2)];
            else
                if pastP2
                    xCurve1 = [x xCurve1 x]; %outside curve
                    yCurve1 = [yzeros(1) yCurve1 yzeros(4)];
                    xCurve4 = [x xCurve4 x]; %bow at far right
                    yCurve4 = [yzeros(2) yCurve4 yzeros(3)]; 
                elseif pastP1
                    xCurve1 = [x xCurve1 x]; %outside curve
                    yCurve1 = [yzeros(1) yCurve1 yzeros(4)];
                    xCurve3 = [x xCurve3 x]; %curve around P2
                    yCurve3 = [yzeros(2) yCurve3 yzeros(3)]; 
                else
                    xCurve1 = [x xCurve1 x]; %outside curve
                    yCurve1 = [yzeros(1) yCurve1 yzeros(4)];
                    xCurve2 = [x xCurve2 x]; %curve around P1
                    yCurve2 = [yzeros(2) yCurve2 yzeros(3)];
                end
            end
        case 2
            if lastNumZeros == 4 && length(yzeros) == 2
                pastP2 = true;
            end
            lastNumZeros = length(yzeros);
            if length(yzeros) == 2
                xCurve1 = [x xCurve1 x]; %outside curve
                yCurve1 = [yzeros(1) yCurve1 yzeros(2)];
            else
                if pastP2
                    xCurve1 = [x xCurve1 x]; %outside curve
                    yCurve1 = [yzeros(1) yCurve1 yzeros(4)];
                    xCurve4 = [x xCurve4 x]; %bow at far right
                    yCurve4 = [yzeros(2) yCurve4 yzeros(3)]; 
                else
                    xCurve1 = [x xCurve1 x]; %outside curve
                    yCurve1 = [yzeros(1) yCurve1 yzeros(4)];
                    xCurve2 = [x xCurve2 x]; %inside curve
                    yCurve2 = [yzeros(2) yCurve2 yzeros(3)];
                end
            end
        case 3
            if length(yzeros) == 2
                xCurve1 = [x xCurve1 x]; %outside curve
                yCurve1 = [yzeros(1) yCurve1 yzeros(2)];
            else
                xCurve1 = [x xCurve1 x]; %outside curve
                yCurve1 = [yzeros(1) yCurve1 yzeros(4)];
                xCurve2 = [x xCurve2 x]; %inside curve, will later match up to the outside curve (really one curve)
                yCurve2 = [yzeros(2) yCurve2 yzeros(3)];
            end  
        case 4
            xCurve1 = [x xCurve1 x]; %top curve (L4)
            yCurve1 = [yzeros(3) yCurve1 yzeros(4)];
            xCurve2 = [x xCurve2 x]; %bottom curve (L5)
            yCurve2 = [yzeros(1) yCurve2 yzeros(2)];
    end
end
switch curves
    case 1 %C3 case, need to combine Curve1 and Curve2
        xCurves{1} = [xCurve1 flip(xCurve2)];
        yCurves{1} = [yCurve1 flip(yCurve2)];
    case 2
        xCurves{1} = [xCurve1 xCurve4]; yCurves{1} = [yCurve1 yCurve4];
        xCurves{2} = [xCurve2 xCurve2(1)]; yCurves{2} = [yCurve2 yCurve2(1)];
    case 3
        xCurves{1} = [xCurve1 xCurve4]; yCurves{1} = [yCurve1 yCurve4];
        xCurves{2} = [xCurve2 xCurve2(1)]; yCurves{2} = [yCurve2 yCurve2(1)];
        xCurves{3} = [xCurve3 xCurve3(1)]; yCurves{3} = [yCurve3 yCurve3(1)];
    case 0
        warning('No ZVC')
end
%%
if isPlot
    pgon = polyshape(xCurves,yCurves);
    hold on
    axis equal
    plot(pgon,'FaceColor',col,'HandleVisibility','off')
end
end
function C = jacobiConstant(varargin)
% C = jacobiConstant(vel,pos,mu3)
% C = jacobiConstant(ST)
% vel is [x'; y'; z'] in nondimensional units, with as many columns as you like
% pos is [x; y; z] in nondimensional units
% mu3 is 3-body mu, or m2/(m1+m2)
% C will be a row vector of C-values
%
% Alex Hoffman
% 9/28/2020
% notation from AAE632 notes
if length(varargin) == 1
    ST_in = varargin{1};
    pos = ST_in.STATE(1:3,:);
    vel = ST_in.STATE(4:6,:);
    mu3 = ST_in.SYSTEM_MODEL.muP1P2;
elseif length(varargin) == 3
    vel = varargin{1};
    pos = varargin{2};
    mu3 = varargin{3};
else
    error('bad inputs')
end

x = pos(1,:); y = pos(2,:); z = pos(3,:);
d = sqrt((x+mu3).^2 + y.^2 + z.^2);
r = sqrt((x-1+mu3).^2 + y.^2 + z.^2);
C = -(vel(1,:).^2 + vel(2,:).^2 + vel(3,:).^2) + x.^2 + y.^2 + 2*(1-mu3)./d + 2*mu3./r;
end