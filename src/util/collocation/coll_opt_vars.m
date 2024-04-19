function [optVars, problem, lb, ub] = coll_opt_vars(problem)
% d_node_times, num_segs_per_traj, init_mass
p = problem.param;
nd_node_times = p.nd_node_times; %just to see the length of this
num_legs = length(problem.legs);
d_node_times = cell(1,num_legs);
num_segs_per_traj = zeros(1,num_legs);

traj = problem.legs{1};
init_mass = traj.getInitMass().value;
optVars = [traj.getInitPos().value; traj.getInitVel().value];
LTlb = [-1;-1;-1;0]; %x,y,z,mass
LTub = [1;1;1;init_mass] + eps;
lb = -inf(6,1);
ub = inf(6,1);

%loop through all legs
for i = 1:num_legs
    traj = problem.legs{i};

    % produce segment nodes and number of segments in this leg
    mesh = problem.mesh{i};
    if isnumeric(mesh) && (iscolumn(mesh) || isrow(mesh))
        meshSegs = length(mesh)-1;
        segTimes = mesh;
    elseif ischar(mesh)
        [meshSegs,method] = strtok(mesh,'_');
        meshSegs = str2double(meshSegs);
        method = method(2:end);
        switch method
            case 'equalTime'
                tspan = traj.getTimeSpan().value;
                segTimes = linspace(tspan(1),tspan(2),meshSegs+1);
            otherwise
                error('unknown mesh method')
        end
    else
        error('unknown mesh')
    end
    num_segs_per_traj(i) = meshSegs;

    timeunit = traj.time.unit;
    segTimes = segTimes(:); %column vector
    defectNodes = 2:2:problem.poly_order;
    for k = 1:meshSegs
        iLT = zeros(4,1);
        inodeStates = [];
        t1 = segTimes(k);
        t2 = segTimes(k+1);
        dimNodes = t1 + (t2-t1)*(nd_node_times+1)/2;
        id_node_times(k,:) = dimNodes;
        iLT(1:3,1) = traj.low_thrust.control_law.coeffs;
        for j = 1:length(nd_node_times)-1
            nodeSpan = dimNodes(j:j+1);
            traj.time = c_dim_quant(nodeSpan,timeunit);
            traj = traj.prop;
            traj.pos = traj.getFinalPos;
            traj.vel = traj.getFinalVel;
            if any(j==defectNodes) %only append if a variable node (end of a propagation from defect node)
                inodeStates = [inodeStates; traj.pos.value];
                inodeStates = [inodeStates; traj.vel.value];
            end
            traj.low_thrust.mass = traj.getFinalMass(); %start with new mass
            if j == length(nd_node_times)-1
                iLT(4,1) = traj.getFinalMass().value;
            end
        end

        %append
        if ~problem.isBallistic
            optVars = [optVars; iLT(:)];
            lb = [lb; LTlb];
            ub = [ub; LTub];
        end
        
        lb = [lb; ones(size(inodeStates)) * (-inf)];
        ub = [ub; ones(size(inodeStates)) * inf];
        optVars = [optVars; inodeStates];
    end
    d_node_times{i} = id_node_times;
end

p.init_mass = init_mass;
p.d_node_times = d_node_times;
p.num_segs_per_traj = num_segs_per_traj;
problem.param = p;
end