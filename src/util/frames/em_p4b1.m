function traj_out = em_p4b1(traj_in)
% traj_out = em_p4b1(traj_in)
%
% Alex Hoffman
% Kenza's masters thesis, page 41 (55 in pdf)
traj_out = traj_in;
sm1 = traj_in.system_model; %orig sysmodel
switch traj_in.system_model.frame
    case 'B1centP1P2rot'
        traj_out.system_model = c_system_model(sm1.system_dynamics,...
            'B2centP4B1rot',sm1.P1,sm1.P2,sm1.P4,'descNodeP4',...
            sm1.char.descNodeP4,'incP4',sm1.char.incP4,...
            'theta0P4',sm1.char.theta0P4);
        sysModelChar = traj_in.system_model.char;
        l_starEM = sysModelChar.lstar.value;
        t_starEM = sysModelChar.tstar.value;
        l_starP4B1 = sysModelChar.lstarP4B1.value; % aSun * l_s
        t_starP4B1 = sysModelChar.tstarP4B1.value;
        mu4 = sysModelChar.muP4B1.value;
        P4B1_2_P1P2 = false;
    case 'B2centP4B1rot'
        traj_out.system_model = c_system_model(sm1.system_dynamics,...
            'B1centP1P2rot',sm1.P1,sm1.P2,sm1.P4,'descNodeP4',...
            sm1.char.descNodeP4,'incP4',sm1.char.incP4,...
            'theta0P4',sm1.char.theta0P4);
        sysModelChar = traj_in.system_model.char;
        l_starEM = sysModelChar.lstarP1P2.value;
        t_starEM = sysModelChar.tstarP1P2.value;
        l_starP4B1 = sysModelChar.lstar.value; % aSun * l_s
        t_starP4B1 = sysModelChar.tstar.value;
        mu4 = sysModelChar.mu.value;
        P4B1_2_P1P2 = true;
end
state = traj_in.getAllStates();
time = traj_in.getAllTime();
oldTdir = traj_in.low_thrust.thrust_dir_history;
angVelP4 = sysModelChar.angVelP4.value; %rad/s
descNodeP4 = sysModelChar.descNodeP4.value;
sunInc = sysModelChar.incP4.value;
theta0P4 = sysModelChar.theta0P4.value;
newState = state;
newTdir = oldTdir;

EM2dim = [eye(3)*l_starEM, zeros(3);
        zeros(3), eye(3)*l_starEM/t_starEM];
dim2P4B1 = [1/l_starP4B1*eye(3), zeros(3);
        zeros(3), t_starP4B1/l_starP4B1*eye(3)];
EM2P4B1 = EM2dim*dim2P4B1;
if P4B1_2_P1P2
    for i = 1:length(time)
       thetaP4 = theta0P4 + time/t_starEM*t_starP4B1*angVelP4;
       thetaP1P2 = pi - thetaP4;
       thEM = thetaP1P2(i);
       X = state(1:6,i)'; %row vector
       k = X - [1-mu4,0,0,0,0,0]; %B2 to B1
       k = k/EM2P4B1;
       C1 = DCM(thEM - descNodeP4,'z',false)';
       C2 = DCM(sunInc,'y',false);
       C3 = DCM(descNodeP4,'z',false);
       C = C1*C2*C3;
       Cd = angVelP4*[sin(thEM - descNodeP4), cos(thEM - descNodeP4), 0;
                  -cos(thEM - descNodeP4), sin(thEM - descNodeP4), 0;
                  0, 0, 0] * C2 * C3; %derivative of C
       Ctilde = [C, Cd'; zeros(3), C];
       STATEOut = k * inv(Ctilde);
       STATEOut = STATEOut'; %column vector
       newState(1:6,i) = STATEOut;
       if ~isempty(oldTdir)
          newTdir(:,i) = (newTdir(:,i)' * transpose(C))'; 
       end
    end
else
    for i = 1:length(time)
       thetaP4 = theta0P4 + time*angVelP4;
       thetaP1P2 = pi - thetaP4;
       thEM = thetaP1P2(i);
       X = state(1:6,i)'; %row vector
       C1 = DCM(thEM - descNodeP4,'z',false)'; %why transposed? because thEM is different from Kenza?
       C2 = DCM(sunInc,'y',false);
       C3 = DCM(descNodeP4,'z',false);
       C = C1*C2*C3;
       Cd = angVelP4*[sin(thEM - descNodeP4), cos(thEM - descNodeP4), 0;
                  -cos(thEM - descNodeP4), sin(thEM - descNodeP4), 0;
                  0, 0, 0] * C2 * C3; %derivative of C
       Ctilde = [C, Cd'; zeros(3), C];
       k = X * Ctilde;
       k = k*EM2P4B1;
       STATEOut = k + [1-mu4,0,0,0,0,0]; %B1 to B2
       STATEOut = STATEOut'; %column vector
       newState(1:6,i) = STATEOut;
       if ~isempty(oldTdir)
          newTdir(:,i) = (newTdir(:,i)' * C)'; 
       end
    end 
end
traj_out.time = traj_out.time.change_unit('sec',sm1);
traj_out.time = traj_out.time.change_unit('nd_t',traj_out.system_model);
traj_out.pos.value = newState(1:3,:);
traj_out.vel.value = newState(4:6,:);
traj_out.low_thrust.thrust_dir_history = newTdir;
end