function [U, Up, Upp,varargout] = ppot_bcr4bp_p1p2(pos,t,mu,mP4,aP4,descNodeP4,incP4,theta0P4,nP4,want)
x = pos(1); y = pos(2); z = pos(3);
r13 = sqrt((x+mu)^2 + y^2 + z^2);
r23 = sqrt((x-1+mu)^2 + y^2 + z^2);
wP4 = nP4 - 1;
thetaP4 = wP4*t + theta0P4;
x4 = aP4*(cos(thetaP4-descNodeP4)*cos(descNodeP4) - sin(thetaP4-descNodeP4)*sin(descNodeP4)*cos(incP4));
y4 = aP4*(cos(thetaP4-descNodeP4)*sin(descNodeP4) + sin(thetaP4-descNodeP4)*cos(descNodeP4)*cos(incP4));
z4 = aP4*sin(thetaP4-descNodeP4)*sin(incP4);
r4 = [x4;y4;z4];
r43 = norm(r4 - pos);
varargout{1} = r4;
varargout{2} = r13;
varargout{3} = r23;
varargout{4} = r43;

U = nan; Up = nan; Upp = nan; %initialize

%U
if want(1)
    U = 0.5*(x^2 + y^2) + (1-mu)/r13 + mu/r23...
        + mP4/r43 - mP4/aP4^3 * dot(r4,pos);
end

%U prime
if want(2)
    Ux = x - (mu+x)*(1-mu)/r13^3 - mu*(x-1+mu)/r23^3 - mP4*(x-x4)/r43^3 - mP4*x4/aP4^3;
    Uy = y - y*(1-mu)/r13^3 - mu*y/r23^3 - mP4*(y-y4)/r43^3 - mP4*y4/aP4^3;
    Uz = -z*(1-mu)/r13^3 - mu*z/r23^3 - mP4*(z-z4)/r43^3 - mP4*z4/aP4^3;   
    Up = [Ux; Uy; Uz];
end

%U prime prime
if want(3)
    Uxx = 1 + 3*(mu+x)^2*(1-mu)/r13^5 - (1-mu)/r13^3 ...
            + 3*mu*(x-1+mu)^2/r23^5 - mu/r23^3 ...
            + 3*(x-x4)^2*mP4/r43^5 - mP4/r43^3;
    Uxy = 3*y*(mu+x)*(1-mu)/r13^5 + 3*mu*y*(x-1+mu)/r23^5 ...
            + 3*mP4*(y-y4)*(x-x4)/r43^5;
    Uxz = 3*z*(mu+x)*(1-mu)/r13^5 + 3*mu*z*(x-1+mu)/r23^5 ...
            + 3*mP4*(z-z4)*(x-x4)/r43^5;
    Uyx = Uxy;
    Uyy = 1 + 3*y^2*(1-mu)/r13^5 - (1-mu)/r13^3 ...
            + 3*mu*y^2/r23^5 - mu/r23^3 ...
            + 3*mP4*(y-y4)^2/r43^5 - mP4/r43^3;
    Uyz = 3*y*z*(1-mu)/r13^5 + 3*mu*y*z/r23^5 ...
            + 3*mP4*(y-y4)*(z-z4)/r43^5;
    Uzx = Uxz;
    Uzy = Uyz;
    Uzz = 3*z^2*(1-mu)/r13^5 - (1-mu)/r13^3 ...
            + 3*mu*z^2/r23^5 - mu/r23^3 ...
            + 3*mP4*(z-z4)^2/r43^5 - mP4/r43^3;
    % need factors of nP4, since thetaP4 = nP4*t
    % derivatives ought to show that (d/dthetaP4) = (d/dnP4*t) = 1/nP4*(d/dt)
    % actually maybe not????
%     Uxt = mP4 * aP4 * (-sin(thetaP4-descNodeP4)*cos(descNodeP4) ...
%             - cos(thetaP4-descNodeP4)*sin(descNodeP4)*cos(incP4))...
%             * (1/r43^3 - 3*(x-x4)^2/r43^5);
%     Uyt = mP4 * aP4 * (-sin(thetaP4-descNodeP4)*sin(descNodeP4) ...
%             + cos(thetaP4-descNodeP4)*cos(descNodeP4)*cos(incP4))...
%             * (1/r43^3 - 3*(y-y4)^2/r43^5);
%     Uzt = mP4 * aP4 * (cos(thetaP4-descNodeP4)*sin(incP4))...
%             * (1/r43^3 - 3*(z-z4)^2/r43^5);
%     Upp = [Uxx Uxy Uxz Uxt; Uyx Uyy Uyz Uyt; Uzx Uzy Uzz Uzt];
    Upp = [Uxx Uxy Uxz; Uyx Uyy Uyz; Uzx Uzy Uzz];
end
end