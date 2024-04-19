function traj_out = eci_ecef(traj_in)
% traj_out = eci_ecef(traj_in)
% TIME is in seconds
traj_out = traj_in;
switch traj_in.system_model.frame
    case 'P1centinert'
        ECEF2ECI = false;
        traj_out.system_model.frame = 'PcentPfixed';
    case 'PcentPfixed'
        ECEF2ECI = true;
        traj_out.system_model.frame = 'P1centinert';
end
time = traj_in.time;
pos = traj_in.pos.value;
% oldTdir = ST_in.THRUST_DIRECTION_HISTORY;
time = time.change_unit('sec', traj_in.system_model);
time = time.value;
EarthAngularVelocity = 2*pi/86164; %rad/s, sidereal day
EarthRotationAngle = EarthAngularVelocity.*time; %radians
% newTdir = oldTdir;
for i = 1:length(time) %for doing many calculations at once
    if ECEF2ECI
        rotmat = DCM(EarthRotationAngle(i),'z',true)';
    else %eci2ecef
        rotmat = DCM(EarthRotationAngle(i),'z',true);
    end
    pos(:,i) = rotmat*pos(:,i);
%     newTdir(:,i) = rotmat*newTdir(:,i);
end
traj_out.pos.value = pos;
% ST_out.THRUST_DIRECTION_HISTORY = newTdir;
end