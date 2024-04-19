function [ H ] = plot3DBody(ax, body, r, pos, varargin)
%PLOT3DBODY Plot a 3D representation of a body
%   
%   H = plot3DBody(ax, BODY, R, POS) takes the name of the body (e.g. 'Earth',
%   'Moon', 'Sun'), its radius R, and its position in the plot POS. Note
%   that R and POS should have the same units for accurate representation.
%   The name of the body is used to assign texture; R and BODY are not
%   necessarily linked. The function returns a handle to the body, H.
%
%   H = plot3DBody(ax BODY, R, POS, N) uses N+1 points to create the sphere
%   object for the body.
%   
%   Author: Andrew Cox
%   Version: December 1, 2014

% Add the directory of images for textures
p = mfilename('fullpath');
sepIx = strfind(p, filesep);
path = [p(1:sepIx(end)), 'bodyTextures'];

col_gray = [125 125 125]/255;

N = 50;
if(~isempty(varargin))
    try
        N = varargin{1};
    catch
    end
end

% Create a sphere and plot it with default settings
[sx, sy, sz] = sphere(N);
H = surf(ax, sx*r + pos(1), sy*r + pos(2), -sz*r + pos(3));
set(H, 'FaceColor', col_gray, 'EdgeColor', 'none', 'FaceLighting', 'phong');

switch lower(body)
    case 'sun'
        topo = imread([path, filesep, 'sun.jpg']);
    case 'mercury'
        topo = imread([path, filesep, 'mercury.jpg']);
    case 'venus'
        topo = imread([path, filesep, 'venus.jpg']);
    case 'earth'
        topo = imread([path, filesep, 'earth.jpg']);
    case 'earth-noclouds'
        topo = imread([path, filesep, 'earth-noclouds.bmp']);
    case 'moon'
        topo = imread([path, filesep, 'moon.jpg']);
    case 'mars'
        topo = imread([path, filesep, 'mars.jpg']);
    case 'jupiter'
        topo = imread([path, filesep, 'jupiter.jpg']);
    case 'saturn'
        topo = imread([path, filesep, 'saturn.jpg']);
    case 'neptune'
        topo = imread([path, filesep, 'neptune.jpg']);
    case 'uranus'
        topo = imread([path, filesep, 'uranus.jpg']);
    case 'pluto'
        topo = imread([path, filesep, 'pluto.jpg']);
    otherwise
        topo = imread([path, filesep, 'asteroid.jpg']);
end

set(H, 'FaceColor', 'texture', 'Cdata', topo);

end

