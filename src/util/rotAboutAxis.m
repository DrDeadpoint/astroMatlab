function rotm = rotAboutAxis(angle,axis)
%rotm = rotAboutAxis(angle,axis)
%   angle: rotation angle in radians
%   axis: vector defining rotation axis
axis = axis/norm(axis); %unit vector
C = cos(angle);
S = sin(angle);
t = 1 - C;
x = axis(1); y = axis(2); z = axis(3);
rotm = [
  t*x^2+C, t*x*y-S*z, t*x*z+S*y;  
  t*x*y+S*z, t*y^2+C, t*y*z-S*x;
  t*x*z-S*y, t*y*z+S*x, t*z^2+C;
];
end