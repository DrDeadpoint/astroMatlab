function FPA = flightPathAngle(r,v)
%FPA = flightPathAngle(r,v)
h = cross(r,v); %angular momentum
b = cross(h,r); %local horizontal
FPA = acos(dot(v,b)/norm(v)/norm(b));
if dot(r,v) < 0
   FPA = -FPA; %below horizontal 
end
end