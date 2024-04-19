function [U, Up, Upp] = ppot_bcr4bp_p4b1(pos,t,mu,mu12,aP4,descNodeP4,...
    incP4,theta0P4,nP4,tstar,tstarP1P2,want)
x = pos(1); y = pos(2); z = pos(3);
wP4 = nP4 - 1;
t = t*tstar.value/tstarP1P2.value;
thetaP4 = wP4*t + theta0P4;
r1 = [1-mu-mu12/aP4*cos(thetaP4); -mu12/aP4*sin(thetaP4);0];
r2 = [1-mu-(1-mu12)/aP4*cos(thetaP4); (1-mu12)/aP4*sin(thetaP4);0];
r13 = pos - r1;
r13n = norm(r13);
r23 = pos - r2;
r23n = norm(r23);
r4 = [-mu;0;0];
r43 = pos - r4;
r43n = norm(r43);

U = nan; Up = nan; Upp = nan; %initialize

%U
if want(1)
    U = 0.5*(x^2 + y^2) + (1-mu)/r43n + mu*((1-mu12)/r13n + mu12/r23n);
end

%U prime
if want(2)
    Ux = x - (mu+x)*(1-mu)/r43n^3 - mu*((1-mu12)*r13(1)/r13n^3 + mu12*r23(1)/r23n^3);
    Uy = y - y*(1-mu)/r43n^3 - mu*((1-mu12)*r13(2)/r13n^3 + mu12*r23(2)/r23n^3);
    Uz = -z*(1-mu)/r43n^3 - mu*((1-mu12)*r13(3)/r13n^3 + mu12*r23(3)/r23n^3);   
    Up = [Ux; Uy; Uz];
end

%U prime prime
if want(3)
    Uxx = 1 + 3*(mu+x)^2*(1-mu)/r43n^5 - (1-mu)/r43n^3 ...
            -mu*( (1-mu12)*(-3*r13(1)^2/r13n^5 + 1/r13n^3) + ...
                  mu12*(-3*r23(1)^2/r23n^5 + 1/r23n^3) );
    Uxy = 3*y*(mu+x)*(1-mu)/r43n^5 + 3*mu* ...
            ((1-mu12)*(r13(1)*r13(2)/r13n^5 + mu12*r23(1)*r23(2)/r23n^5));
    Uxz = 3*z*(mu+x)*(1-mu)/r43n^5 + 3*mu* ...
            ((1-mu12)*(r13(1)*r13(3)/r13n^5 + mu12*r23(1)*r23(3)/r23n^5));
    Uyx = Uxy;
    Uyy = 1 + 3*y^2*(1-mu)/r43n^5 - (1-mu)/r43n^3 ...
            + mu*( (1-mu12)*(-3*r13(2)^2/r13n^5 + 1/r13n^3) + ...
                  mu12*(-3*r23(2)^2/r23n^5 + 1/r23n^3) );
    Uyz = 3*y*z*(1-mu)/r43n^5 + 3*mu* ...
            ((1-mu12)*(r13(2)*r13(3)/r13n^5 + mu12*r23(2)*r23(3)/r23n^5));
    Uzx = Uxz;
    Uzy = Uyz;
    Uzz = 3*z^2*(1-mu)/r43n^5 - (1-mu)/r43n^3 ...
            + mu*( (1-mu12)*(-3*r13(3)^2/r13n^5 + 1/r13n^3) + ...
                  mu12*(-3*r23(3)^2/r23n^5 + 1/r23n^3) );
    Upp = [Uxx Uxy Uxz; Uyx Uyy Uyz; Uzx Uzy Uzz];
end
end