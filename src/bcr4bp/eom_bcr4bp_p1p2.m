function [xdd, ydd, zdd, varargout] = eom_bcr4bp_p1p2(state,mu3,t,t0,m4,a4,descNodeP4,incP4,theta0P4,nP4,Tmax,controlLaw,wantUpp)
% state is 6x1
% mu is 3B mu
% m is mass at time t
% t is time
% t0 is initial time of propagation
% Tmax is max thrust (nd) --- 0 if ballistic
% controlLaw is a c_control_law containing claw coeffs
% wantUpp is logical stating whether Upp is wanted
% Alex Hoffman
% 03/21/22
if size(state,1) ~= 7 || size(state,2) ~= 1
    error('Input state must be 7x1')
end
varargout{1} = [];
pos = state(1:3);
x = pos(1); y = pos(2); z = pos(3);
xd = state(4); yd = state(5); zd = state(6);
vel = [xd;yd;zd];
m = state(7);

%% low thrust section
if Tmax > 0
    switch controlLaw.law
        case 'SOC'
            cf = controlLaw.coeffs;
            dt = t-t0;
            alpha = cf.alpha0 + cf.alphadot*dt + cf.alphaddot*dt^2/2 + cf.alphaAmp*sin(cf.alphaFreq*dt + cf.alphaPhase);
            beta = cf.beta0 + cf.betadot*dt + cf.betaddot*dt^2/2 + cf.betaAmp*sin(cf.betaFreq*dt + cf.betaPhase);
            uhat = [cos(alpha)*cos(beta);...
                    sin(alpha)*cos(beta);...
                    sin(beta)];
        case 'fixed_dir'
            uhat = controlLaw.coeffs;
        otherwise
            error('unknown control law')
    end
    switch controlLaw.frame
        case 'B1centP1P2rot' % no rotations needed
            %good to go
        case 'B1centP1P2rot_VUW' %VUW defined using same r and v as propagation frame
            angMom = cross(pos,vel);
            V = vel/norm(vel);
            W = angMom/norm(angMom);
            U = cross(W,V);            
            VUW_R_EMrot = [V, U, W];
            uhat = VUW_R_EMrot * uhat;
        otherwise
            error(['Control law frame **' controlLaw.frame '** not implemented.'])
    end
    if abs(norm(uhat) - 1) > 1e-5
       error('Uhat magnitude too large') 
    end
    amag = Tmax/m;
    acc_LT = amag*uhat;
else
    uhat = [0;0;0];
    acc_LT = [0;0;0];
end

%% calculate accelerations
if wantUpp
    [~,Up,Upp,r4,r13,r23,r43] = ppot_bcr4bp_p1p2(pos,t,mu3,m4,a4,descNodeP4,incP4,theta0P4,nP4,[0 1 1]); %calculate pseudopotential
    varargout{1} = Upp;
    sP1 = [-mu3;0;0;0;0;0];
    sP2 = [1-mu3;0;0;0;0;0];
    v4 = a4*[-sin(nP4*t-descNodeP4)*nP4*cos(descNodeP4) - cos(nP4*t-descNodeP4)*nP4*sin(descNodeP4)*cos(incP4);
        -sin(nP4*t-descNodeP4)*nP4*sin(descNodeP4) + cos(nP4*t-descNodeP4)*nP4*cos(descNodeP4)*cos(incP4)
        cos(nP4*t-descNodeP4)*nP4*sin(incP4)];
    sP4 = [r4;v4];
    varargout{4} = sP1;
    varargout{5} = sP2;
    varargout{6} = sP4;
    dsd_dr1 = zeros(6,3);
    dsd_dr2 = dsd_dr1;
    dsd_dr4 = dsd_dr1;
    %dsd_dr1
    x1 = -mu3;
    dsd_dr1(4,1) = -(1-mu3)*(-1/r13^3 + 3*(x-x1)/r13^5);
    dsd_dr1(5,1) = -(1-mu3)*3*y*(x-x1)/r13^5;
    dsd_dr1(6,1) = -(1-mu3)*3*z*(x-x1)/r13^5;
    %dsd_dr2
    x2 = 1-mu3;
    dsd_dr2(4,1) = -mu3*(-1/r23^3 + 3*(x-x2)/r23^5);
    dsd_dr2(5,1) = -(1-mu3)*3*y*(x-x2)/r23^5;
    dsd_dr2(6,1) = -(1-mu3)*3*z*(x-x2)/r23^5;
    %dsd_dr4
    x4 = r4(1); y4 = r4(2); z4 = r4(3);
    dsd_dr4(4,1) = -m4*(-1/r43^3 + 3*(x-x4)^2/r43^5 + 1/a4^3);
    dsd_dr4(4,2) = -m4*(x-x4)*3*(y-y4)/r43^5;
    dsd_dr4(4,3) = -m4*(x-x4)*3*(z-z4)/r43^5;
    dsd_dr4(5,1) = dsd_dr4(4,2);
    dsd_dr4(5,2) = -m4*(-1/r43^3 + 3*(y-y4)^2/r43^5 + 1/a4^3);
    dsd_dr4(5,3) = -m4*(y-y4)*3*(z-z4)/r43^5;
    dsd_dr4(6,1) = dsd_dr4(4,3);
    dsd_dr4(6,2) = dsd_dr4(5,3);
    dsd_dr4(6,3) = -m4*(-1/r43^3 + 3*(z-z4)^2/r43^5 + 1/a4^3);
    %output
    varargout{7} = dsd_dr1;
    varargout{8} = dsd_dr2;
    varargout{9} = dsd_dr4;
else
    [~,Up,~] = ppot_bcr4bp_p1p2(pos,t,mu3,m4,a4,descNodeP4,incP4,theta0P4,nP4,[0 1 0]); %calculate pseudopotential
end

xdd = 2*yd + Up(1) + acc_LT(1);
ydd = -2*xd + Up(2) + acc_LT(2);
zdd = Up(3) + acc_LT(3);

varargout{2} = acc_LT;
varargout{3} = uhat;
end