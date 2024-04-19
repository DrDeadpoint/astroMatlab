function traj_out = rot_inert(traj_in)
% traj_out = rot_inert(traj_in)
% rotates into inertial FRAME based on TIME (or other way around)
% ST_in is a STATE-TIME with STATEs to be rotated
% thetadot = 1 for nondimensional mean motion, 1/t_star for dimensional motion
% looks at ST_in.FRAME to see how to rotate
%
% Alex Hoffman
% 06/29/21
traj_out = traj_in;
pos = traj_in.pos;
vel = traj_in.vel;
tvec = traj_in.time;
thetadot = 1;
switch traj_in.system_model.frame
    case 'B2centP4B1inert'
        outFrame = 'B2centP4B1rot';
        inert2rot = true;
    case 'B2centP4B1rot'
        outFrame = 'B2centP4B1inert';
        inert2rot = false;
    case 'B1centP1P2inert'
        outFrame = 'B1centP1P2rot';
        inert2rot = true;
    case 'B1centP1P2rot'
        outFrame = 'B1centP1P2inert';
        inert2rot = false;
    case 'P2centinert' 
        outFrame = 'P2centP1P2rot';
        inert2rot = true;
    case 'P2centP1P2rot' 
        outFrame = 'P2centinert';
        inert2rot = false;
    case 'P1centinert'
        outFrame = 'P1centP1P2rot';
        inert2rot = true;
    case 'P1centP1P2rot'
        outFrame = 'P1centinert';
        inert2rot = false;
    otherwise
        error('Frame rotation not found')
end
traj_out.system_model.frame = outFrame;
tvec = tvec.change_unit('nd_t', traj_in.system_model);
pos = pos.change_unit('nd_l', traj_in.system_model);
vel = vel.change_unit('nd_v', traj_in.system_model);
[~,c] = size(pos.value);
if inert2rot
    angVelSys = -thetadot;
else
    angVelSys = thetadot;
end
OMEGA_ORIG_IN_NEW = [0;0;angVelSys];
outpos = zeros(3,c);
outvel = zeros(3,c);
newTdir = traj_in.low_thrust.thrust_dir_history;
for j = 1:c
    thisT = tvec.value(j);
    rotmat = DCM(thisT*angVelSys,'z',false); %false for negative angle to move backwards about z
    thisPos = pos.value(:,j);
% ATD method (Follows most people)
    outpos(:,j) = rotmat*thisPos;
    thisVel = vel.value(:,j);
    newVel = thisVel + cross(OMEGA_ORIG_IN_NEW,thisPos); %original FRAME unit vectors
    outvel(:,j) = rotmat*newVel; %new FRAME unit vectors
% %% What I think might be the correct method
%     newPos = rotmat*thisPos + patch_val;
%     thisState = newPos;
%     if r==6
%         thisVel = inStates(4:6,j);
%         newVel = rotmat'*v_cent_bary + thisVel; %puts it into inertial FRAME fixed at barycenter
%         newVel = rotmat*(newVel + cross(OMEGA_ORIG_IN_NEW,newPos));
%         thisState = [thisState; newVel];
%     end
%     outStates(:,j) = thisState;

    if ~isempty(traj_in.low_thrust.thrust_dir_history)
       newTdir(:,j) = rotmat*newTdir(:,j); 
    end
end
pos.value = outpos;
vel.value = outvel;
traj_out.pos = pos;
traj_out.vel = vel;
traj_out.time = tvec;
traj_out.low_thrust.thrust_dir_history = newTdir;
end