function [value,isterminal,direction] = phase2event(t,y,mu,rP1,rP2,rMax)
p1 = [-mu;0;0];
p2 = [1-mu;0;0];
bary = [0;0;0];
pos = y(1:3);
prim1dist = norm(pos - p1) - rP1;
prim2dist = norm(pos - p2) - rP2;
barydist = norm(pos - bary) - rMax;
value = [prim1dist, prim2dist, barydist];  
isterminal = [1, 1, 1];
direction  = [0, 0, 0];
end
