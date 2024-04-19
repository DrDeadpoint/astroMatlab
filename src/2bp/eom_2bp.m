function [xdd, ydd, zdd, varargout] = eom_2bp(state,mu,J2,Re,t,t0,Tmax,controlLaw,wantUpp)
% state is 6x1
% mu is 2B mu
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
r = norm(pos);
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
        case 'P1centinert' % no rotations needed
            %good to go
        case 'P1centinert_VUW' %VUW defined using same r and v as propagation frame
            angMom = cross(pos,vel);
            V = vel/norm(vel);
            W = angMom/norm(angMom);
            U = cross(W,V);            
            VUW_R_EMrot = [V, U, W];
            uhat = VUW_R_EMrot * uhat;
            error('Need to get correct Upp since EOMs have r3,v3')
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
% if wantUpp
%     [~,Up,Upp] = ppot_cr3bp(pos,mu,[0 1 1]); %calculate pseudopotential
%     varargout{1} = Upp;
% else
%     [~,Up,~] = ppot_cr3bp(pos,mu,[0 1 0]); %calculate pseudopotential
% end

if J2 ~= 0
    a_J2 = -1.5*J2*(mu/r^2)*(Re/r)^2.*...
        [(1-5*(z/r)^2)*x/r;...
        (1-5*(z/r)^2)*y/r;...
        (3-5*(z/r)^2)*z/r];
else
    a_J2 = 0;
end
dd = -mu/(norm(pos)^3)*pos + a_J2;

if wantUpp
    Upp = [-mu/r^3 + 3*mu*x^2/r^5,  3*mu*x*y/r^5,           3*mu*x*z/r^5
           3*mu*y*x/r^5,            -mu/r^3 + 3*mu*y^2/r^5, 3*mu*y*z/r^5
           3*mu*z*x/r^5,            3*mu*z*y/r^5,           -mu/r^3 + 3*mu*z^2/r^5
          ];
    varargout{1} = Upp;
end

xdd = dd(1) + acc_LT(1);
ydd = dd(2) + acc_LT(2);
zdd = dd(3) + acc_LT(3);

varargout{2} = acc_LT;
varargout{3} = uhat;
end