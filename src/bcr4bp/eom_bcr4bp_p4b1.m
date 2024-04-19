function [xdd, ydd, zdd, varargout] = eom_bcr4bp_p4b1(state,mu,t,t0,mu12,a4,...
    descNodeP4,incP4,theta0P4,nP4,Tmax,controlLaw,tstar,tstarP1P2,wantUpp)
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
xd = state(4); yd = state(5); zd = state(6);
vel = [xd;yd;zd];
m = state(7);

%% low thrust section
if Tmax > 0
    error('something is wrong here, was getting really large acc')
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
        case 'B2centP4B1rot' % no rotations needed
            %good to go
        case 'B2centP4B1rot_VUW' %VUW defined using same r and v as propagation frame
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
if wantUpp
    [~,Up,Upp] = ppot_bcr4bp_p4b1(pos,t,mu,mu12,a4,descNodeP4,incP4,...
        theta0P4,nP4,tstar,tstarP1P2,[0 1 1]); %calculate pseudopotential
    varargout{1} = Upp;
else
    [~,Up,~] = ppot_bcr4bp_p4b1(pos,t,mu,mu12,a4,descNodeP4,incP4,...
        theta0P4,nP4,tstar,tstarP1P2,[0 1 0]); %calculate pseudopotential
end

xdd = 2*yd + Up(1) + acc_LT(1);
ydd = -2*xd + Up(2) + acc_LT(2);
zdd = Up(3) + acc_LT(3);

varargout{2} = acc_LT;
varargout{3} = uhat;
end