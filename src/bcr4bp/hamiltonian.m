function [ham,varargout] = hamiltonian(traj,hamInWhichFrame)
% [ham,varargout] = hamiltonian(traj,hamInWhichFrame)
% hamInWhichFrame lets you change which frame you get the hamiltonian in
% subtracting off solarEffects should yield a scaled JC
% Alex Hoffman
if nargin == 1
   hamInWhichFrame = traj.system_model.frame; 
end
state = [traj.pos.value; traj.vel.value];
time = traj.time.value;
[r,~] = size(time);
if r > 1
   time = time';
   traj.time.value = time;
end
sysModel = traj.system_model.char;
ms = sysModel.massP4.value;
as = sysModel.distP4.value;
sunTheta0 = sysModel.theta0P4.value;
mu3 = sysModel.mu.value;
n_sb1 = sqrt((ms+1) / as^3);
sunThetaDot = n_sb1 - 1;
switch hamInWhichFrame
    case 'B1centP1P2rot'
        if strcmp(traj.system_model.frame,'B2centP4B1rot')
            ST_EM = CONVframe(traj,'B1centP1P2rot');
            state = ST_EM.state;
            time = ST_EM.time;
            [r,~] = size(time);
            if r > 1
               time = time';
               traj.TIME = time;
            end
        end
        sunAngle = sunTheta0 + time*sunThetaDot;
        sunX = as*cos(sunAngle);
        sunY = as*sin(sunAngle);
        x = state(1,:); y = state(2,:); z = state(3,:);
        xd = state(4,:); yd = state(5,:); zd = state(6,:);
        r13 = sqrt((x+mu3).^2 + y.^2 + z.^2);
        r23 = sqrt((x-1+mu3).^2 + y.^2 + z.^2);
        %get location of sun at each time
        r43 = sqrt((x - sunX).^2 + (y - sunY).^2 + z.^2);
        solarEffects = 2*(ms./r43 - ms/as^3.*(x.*sunX + y.*sunY));
        pseudoPotential = (x.^2 + y.^2) + 2*(1-mu3)./r13 + 2*mu3./r23 + solarEffects; 
        Hdot = 2.*sunThetaDot.*ms.*as.*(-x.*sin(sunAngle) + y.*cos(sunAngle)).*(1./r43.^3 - 1./as.^3); %scheurle MS pg 114
        varargout{1} = Hdot;
    case 'B2centP4B1rot'
        if strcmp(traj.system_model.frame,'B1centP1P2rot')
            ST_sb1 = CONVframe(traj,'B2centP4B1rot');
            state = ST_sb1.state;
            time = ST_sb1.time;
            [r,~] = size(time);
            if r > 1
               time = time';
               traj.TIME = time;
            end
        end
        x = state(1,:); y = state(2,:); z = state(3,:);
        xd = state(4,:); yd = state(5,:); zd = state(6,:);
        EMtheta0 = pi - sunTheta0;
        mu4 = sysModel.muP4B1;
        thetaEM = EMtheta0 + time*sunThetaDot;
        x1 = 1 - mu4 - 1/as*mu3*cos(thetaEM);
        y1 = -1/as*mu3*sin(thetaEM);
        x2 = 1 - mu4 + 1/as*(1-mu3)*cos(thetaEM);
        y2 = 1/as*(1-mu3)*sin(thetaEM);
        r43 = sqrt((x-mu4).^2 + y.^2 + z.^2);
        r13 = sqrt((x-x1).^2 + (y-y1).^2 + z.^2);
        r23 = sqrt((x-x2).^2 + (y-y2).^2 + z.^2);
        pseudoPotential = (x.^2 +y.^2) + 2*(1-mu4)./r43 + 2*mu4*(1-mu3)./r13 + 2*mu4*mu3./r23;
end
ham = -(xd.^2 + yd.^2 + zd.^2) + pseudoPotential;
end