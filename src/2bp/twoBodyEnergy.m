function E = twoBodyEnergy(r,v,mu)
r = norm(r); v = norm(v);
SMA = 1/(2/r - v^2/mu);
E = -mu/2/SMA;
end