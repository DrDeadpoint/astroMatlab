function [U, Up, Upp] = ppot_cr3bp(pos,mu,want)
x = pos(1); y = pos(2); z = pos(3);
d = sqrt((x+mu)^2 + y^2 + z^2);
r = sqrt((x-1+mu)^2 + y^2 + z^2);

U = nan; Up = nan; Upp = nan; %initialize

%U
if want(1)
    U = 0.5*(x^2 + y^2) + (1-mu)/d + mu/r;
end

%U prime
if want(2)
    Ux = x - (mu+x)*(1-mu)/d^3 - mu*(x-1+mu)/r^3;
    Uy = y - y*(1-mu)/d^3 - mu*y/r^3;
    Uz = -z*(1-mu)/d^3 - mu*z/r^3;   
    Up = [Ux; Uy; Uz];
end

%U prime prime
if want(3)
    Uxx = 1 + 3*(mu+x)^2*(1-mu)/d^5 - (1-mu)/d^3 + 3*mu*(x-1+mu)^2/r^5 - mu/r^3;
    Uxy = 3*y*(mu+x)*(1-mu)/d^5 + 3*mu*y*(x-1+mu)/r^5;
    Uxz = 3*z*(mu+x)*(1-mu)/d^5 + 3*mu*z*(x-1+mu)/r^5;
    Uyx = Uxy;
    Uyy = 1 + 3*y^2*(1-mu)/d^5 - (1-mu)/d^3 + 3*mu*y^2/r^5 - mu/r^3;
    Uyz = 3*y*z*(1-mu)/d^5 + 3*mu*y*z/r^5;
    Uzx = Uxz;
    Uzy = Uyz;
    Uzz = 3*z^2*(1-mu)/d^5 - (1-mu)/d^3 + 3*mu*z^2/r^5 - mu/r^3;
    Upp = [Uxx Uxy Uxz; Uyx Uyy Uyz; Uzx Uzy Uzz];
end
end