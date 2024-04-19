function [L1, L2, L3, L4, L5, varargout] = lagrangePoints(mu3,varargin)
% [L1, L2, L3, L4, L5, varargout] = lagrangePoints(mu3,varargin)
% mu3 = 3-body mu, or m2/(m1+m2)
% Li are xyz coordinates of lagrange points
%
% varargin:
%   (...,'dim', lstar,...) : dimensional quantities
%       where lstar is the characteristic length of the system
%   [...,gamOut] = lagrangePoints(...,'gamma',...)
%       yields the gamma1, gamma2, gamma3 values
%
% Alex Hoffman
% 9/28/2020
dim = false;
gamOut = false;
for i = 1:length(varargin)/2
   switch varargin{i*2-1}
       case 'dim'
           lstar = varargin{i*2};
           dim = true;
       case 'gamma'
           gamOut = true;
       otherwise
           error('variable argument not found')
   end
end
gamOutVec = [0 0 0]; %initialize
if mu3 < 0 || mu3 > 1
error('mu3 must be between 0 and 1')
end

%% L2
funcL2 = @(g2) (1-mu3)/(1+g2)^2 + mu3/g2^2 - 1 + mu3 - g2;
dfuncL2 = @(g2) -1 - 2*mu3/g2^3 + 2*(mu3-1)/(1+g2)^3;
[g2,~] = newtonMethod(funcL2,dfuncL2,0.5,10^-8,'diff'); %nondimensional
L2 = 1 - mu3 + g2; %nondimensional
if dim
    L2 = L2*lstar; %dimensional, meters
    g2 = g2*lstar; %dimensional, meters
end
gamOutVec(2) = g2;
%% L1
funcL1 = @(g1) -((1-mu3)/(1-g1)^2 - mu3/g1^2 + g1 + mu3 - 1);
% dfuncL1 = @(g1) -(1 + 2*mu3/g1^3 - 2*(1-mu3)/(1-g1)^3);
g1 = fzero(funcL1,1.5);
% [g1nd,~] = newtonMethod(funcL1,dfuncL1,0.5,10^-8,'diff'); %nondimensional
L1 = 1-mu3-g1; %nondimensional
if dim
    L1 = L1*lstar; %dimensional, meters
    g1 = g1*lstar; %dimensional, meters
end
gamOutVec(1) = g1;
%% L3
funcL3 = @(g3) (1-mu3)/g3^2 + mu3/(-g3-1)^2 - mu3 - g3;
% dfuncL3 = @(g3) -1 + 2*mu3/(-g3-1)^3 - 2*(1-mu3)/g3^3;
g3 = fzero(funcL3,0.8);
% [g3nd,~] = newtonMethod(funcL1,dfuncL1,0.8,10^-8,'diff'); %nondimensional
L3 = -mu3-g3; %nondimensional
if dim
    L3 = L3*lstar; %dimensional, meters
    g3 = g3*lstar; %dimensional, meters
end
gamOutVec(3) = g3;
%% L4 and L5
a = 1; %semi major axis
B = mu3*a; %barycenter = (m2*a/(m1+m2))
%equilateral triangles
h = tand(60)*a/2;
L4 = [a/2 - B; h; 0];
L5 = [a/2 - B; -h; 0];
if dim
    L4 = L4.*lstar;
    L5 = L5.*lstar;
end

%% outputs
L1 = [L1; 0; 0];
L2 = [L2; 0; 0];
L3 = [L3; 0; 0];
if gamOut
    varargout{1} = gamOutVec;
end
end