function [variableStates,defectStates,defectDerivs,traj,m_i,Ci] = coll_def_nodes(LTvars,varStates,problem,traj,m_i,d_dt)
p = problem.param;
Ainv = p.Ainv;
B = p.B;
D = p.D;
if problem.isBallistic
    traj.low_thrust.spacecraft.throttle = 0;
else
    %define low thrust
    traj.low_thrust.control_law.coeffs = LTvars(1:3)./norm(LTvars(1:3)); %normalize so it doesn't error
    m_f = LTvars(4);
    traj.low_thrust.mass.value = [m_i, m_f];
    throttle = (m_i - m_f) / (traj.low_thrust.spacecraft.mdotND.value*d_dt);
    traj.low_thrust.spacecraft.throttle = throttle; %define throttle based on mass variables
    m_i = m_f; %for next segment
end
%define states ----
variableNodes = 1:2:problem.poly_order;
variableStates = zeros(6,length(variableNodes));
variableDerivs = zeros(6,length(variableNodes));
for j = 1:length(variableNodes) %build variableStates and variableDerivs
    traj.pos.value = varStates(1:3);
    traj.vel.value = varStates(4:6);
    if j ~= length(variableNodes)
        varStates(1:6) = []; %pop off this node, unless using for next segment
    end
    variableStates(1:3,j) = traj.pos.value;
    variableStates(4:6,j) = traj.vel.value;
    traj.time.value = p.nd_node_times(variableNodes(j));
    variableDerivs(1:3,j) = traj.vel.value;
    variableDerivs(4:6,j) = traj.EOM();
end % yields 9xorder matrix of nodes along a segment
variableDerivs = variableDerivs .* d_dt/2; %normalize derivatives wrt time
Ci = [variableStates, variableDerivs] * Ainv;
defectStates = Ci*B;
defectStates = defectStates(:,2:end-1);
defectDerivs = Ci*D;
end