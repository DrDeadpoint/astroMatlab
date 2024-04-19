function conFunc = coll_con_func(problem)
conFunc = @(x) con(x,problem);
%for each segment, get variable and defect node states
end

function [ineqCons, eqCons] = con(optVars, problem)
p = problem.param;
W = p.W;
nd_node_times = p.nd_node_times;
d_node_times = p.d_node_times;
desInitNode = problem.desInitNode;
desFinalNode = problem.desFinalNode;
variableNodeInds = 1:2:problem.poly_order;
defectNodeInds = 2:2:problem.poly_order;
defectEqns = [];
initState = optVars(1:6);
m0 = p.init_mass;
uhatEqns = [];
thrEqns = [];
for i = 1:length(problem.legs)
    traj = problem.legs{i};
    id_node_times = d_node_times{i};
    for k = 1:p.num_segs_per_traj(i)
        tspan = id_node_times(k,[1 end]); %dt is per segment
        d_dt = tspan(2) - tspan(1); %dimensional dt

        % draw out the optimization variables for this segment
        if problem.isBallistic
            LTvars = [];
        else
            LTvars = optVars(7:10); optVars = optVars([1:6,11:end]);
        end
        varStates = optVars(1:6*length(variableNodeInds));
        optVars(1:6*(length(variableNodeInds)-1)) = []; %leave last node on there for the next segment

        % get defect states, replace traj and m0
        [~,defectStates,defectDerivs,traj,m0] = coll_def_nodes(LTvars,varStates,problem,traj,m0,d_dt);
        if ~problem.isBallistic
            % add some low-thrust constraint equations
            uhatEqns = [uhatEqns; norm(LTvars(1:3)) - 1];
            throttle = traj.low_thrust.spacecraft.throttle;
            thrEqns = [thrEqns; throttle - 1; -throttle]; % c <= 0
        end

        % create defect node constraint equations
        trueDefectDerivs = zeros(6,length(defectNodeInds));
        for j = 1:length(defectNodeInds)
            trueDefectDerivs(1:3,j) = defectStates(4:6,j);
            traj.pos.value = defectStates(1:3,j);
            traj.vel.value = defectStates(4:6,j);
            traj.time.value = nd_node_times(defectNodeInds(j));
            trueDefectDerivs(4:6,j) = traj.EOM();
            trueDefectDerivs(1:6,j) = trueDefectDerivs(1:6,j)  * d_dt/2;
        end
        defectEqns = [defectEqns; (defectDerivs - trueDefectDerivs)*W];
    end
end
% fix beginning and end points
finalStateEqn = desFinalNode - optVars; %all that is left should be a 6x1 vector
initStateEqn = desInitNode - initState;
finalStateEqn = finalStateEqn(~isnan(finalStateEqn));
initStateEqn = initStateEqn(~isnan(initStateEqn));
eqCons = [defectEqns(:); uhatEqns; finalStateEqn; initStateEqn];
% NOTE: don't need to enforce continuity between segments because I
% manually match the final node from segment_i to the first node on
% segment_i+1

% ineq constraints:
%keep thrust between 0 and 1 for each segment
ineqCons = thrEqns;
%planetary radius minimum
end