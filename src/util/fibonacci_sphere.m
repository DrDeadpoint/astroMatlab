function [points] = fibonacci_sphere(N)
%https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere
points = zeros(3,N);
phi = pi * (3 - sqrt(5));  % golden angle in radians

for i = 0:N-1
    y = 1 - (i/(N - 1)) * 2;  % y goes from 1 to -1
    radius = sqrt(1 - y^2);  % radius at y

    theta = phi * i;  % golden angle increment

    x = cos(theta) * radius;
    z = sin(theta) * radius;

    points(:,i+1) = [x;y;z];
end
end