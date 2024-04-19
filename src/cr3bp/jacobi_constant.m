function C = jacobi_constant(varargin)
% C = jacobi_constant(pos,vel,mu3)
% C = jacobi_constant(traj)
% vel is [x'; y'; z'] in nondimensional units, with as many columns as you like
% pos is [x; y; z] in nondimensional units
% mu3 is 3-body mu, or m2/(m1+m2)
% C will be a row vector of C-values
%
% Alex Hoffman
% 09/04/2022
% notation from AAE632 notes
if length(varargin) == 1
    traj_in = varargin{1};
    if ~strcmp(traj_in.system_model.frame,'B1centP1P2rot')
        error('Must put states into B1centP1P2rot frame')
    end
    traj_in.pos.change_unit('nd_l');
    traj_in.vel.change_unit('nd_v');
    pos = traj_in.pos.value;
    vel = traj_in.vel.value;
    mu3 = traj_in.system_model.char.mu;
elseif length(varargin) == 3
    pos = varargin{1};
    pos = pos.change_unit('nd_l');
    pos = pos.value;
    vel = varargin{2};
    vel = vel.change_unit('nd_v');
    vel = vel.value;
    mu3 = varargin{3};
else
    error('bad inputs')
end

x = pos(1,:); y = pos(2,:); z = pos(3,:);
d = sqrt((x+mu3).^2 + y.^2 + z.^2);
r = sqrt((x-1+mu3).^2 + y.^2 + z.^2);
C = -(vel(1,:).^2 + vel(2,:).^2 + vel(3,:).^2) + x.^2 + y.^2 + 2*(1-mu3)./d + 2*mu3./r;
end