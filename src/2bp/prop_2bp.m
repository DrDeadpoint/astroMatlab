function [traj_out] = prop_2bp(traj_in,options,oblate,prop)
if oblate(2) ~= 0
    warning('STM is incorrect when using J2 terms')
end
traj_out = traj_in;
sysModel = traj_in.system_model;
tspan = traj_in.getTimeSpan().value;
r0 = traj_in.getInitPos().value;
v0 = traj_in.getInitVel().value;
spacecraft = traj_in.low_thrust.spacecraft;
mass0 = traj_in.getInitMass();
mass0 = mass0.change_unit('nd_m',spacecraft);
mass0 = mass0.value;
mu = sysModel.char.mu.value;
Re = oblate(1);
J2 = oblate(2);
Tmax = spacecraft.TmaxND.value;
Isp = spacecraft.IspND.value;
g0 = spacecraft.g0ND.value;
cLaw = traj_in.low_thrust.control_law;
throttle = spacecraft.throttle;
Tmax = Tmax*throttle;
switch prop
    case 'manual'
        STM_0 = eye(7);
        STM_0 = reshape(STM_0,49,1);
        [tout,stateVec] = ode45(@(t,s) twobody(t,s,tspan(1),Tmax,Isp,g0,cLaw),tspan,[r0; v0; mass0; STM_0],options);
        stateVec = transpose(stateVec);
        STM_f = stateVec(8:56,end);
        STM_f = reshape(STM_f,7,7);
        stateVec = stateVec(1:7,:);
    otherwise
        error('Not implemented')
end
traj_out.time = c_dim_quant(tout,'nd_t');
traj_out.pos = c_dim_quant(stateVec(1:3,:),'nd_l');
traj_out.vel = c_dim_quant(stateVec(4:6,:),'nd_v');
traj_out.low_thrust.mass = c_dim_quant(stateVec(7,:),'nd_m');
STMkey = {'x', 'y', 'z', 'xd', 'yd', 'zd', 'mass', 't0'};
STM_f = [STM_f, zeros(7,1); zeros(1,8)]; %for t0
traj_out.stm = c_stm(STM_f,STMkey);
%% Nested Functions for sharing variables
% ----------------------------------------------------
function dsstmdt = twobody(t,s_stm,t0,Tmax,Isp,g0,controlLaw)
    %y = [X V]
    %dydt = [V A]
    dsstmdt = zeros(56,1);
    s = s_stm(1:7);
    stm = s_stm(8:end);
    dsdt(1:3) = s(4:6);

    [xdd, ydd, zdd, A21, acc_LT, uhat] = eom_2bp(s,mu,J2,Re,t,t0,Tmax,controlLaw,true);
    dsdt(4:6) = [xdd;ydd;zdd];
    dsdt(7) = -Tmax/Isp/g0;

    A23 = -acc_LT/s(7); %3x1
    A = [zeros(3,3),eye(3), zeros(3,1);
        A21, zeros(3,3), A23;
        zeros(1,7)];
    stm = reshape(stm,[7 7]);
    dstmdt = A*stm;
    dstmdt = reshape(dstmdt,[49 1]);
    dsstmdt(8:56) = dstmdt;
    dsstmdt(1:7) = dsdt;
end 
end