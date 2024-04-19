function delT = KepTOF(nu1,nu2,e,T)
%delT = KepTOF(nu1,nu2,e,T)
%   nu1 is start angle, in radians
%   nu2 is end angle, in radians
%   e is eccentricity
%   T is time period of orbit

n = 2*pi/T;
if nu1 <= pi
    E1 = acos((e+cos(nu1))/(1+e*cos(nu1)));
else
    E1 = 2*pi - acos((e+cos(nu1))/(1+e*cos(nu1)));
end
M1 = E1-e*sin(E1);
if nu2 <= pi
    E2 = acos((e+cos(nu2))/(1+e*cos(nu2)));
else
    E2 = 2*pi - acos((e+cos(nu2))/(1+e*cos(nu2)));
end
M2 = E2-e*sin(E2);
delT = (M2-M1)/n;
end