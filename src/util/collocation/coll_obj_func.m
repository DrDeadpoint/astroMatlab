function objFunc = coll_obj_func(problem)
p = problem.param;
m0 = p.init_mass;
if problem.isBallistic
    objFunc = @(x) 1;
else
    objFunc = @(x) mass_obj_func(x,m0,problem);
end 
end
%%
function obj = mass_obj_func(optVars,m0,problem)
    p = problem.param;
    variableNodeInds = 1:2:problem.poly_order;
    obj = 0;
    for j = 1:length(problem.legs)
        for i = 1:p.num_segs_per_traj(j) %number of segments
            LTdetails = optVars(7:10);
            m_f = LTdetails(4);
            obj = obj + (m0 - m_f);
            m0 = m_f;
            optVars(1:6*(length(variableNodeInds)-1)+4) = []; % -1 leave one node before LTdetails
        end
    end
    %     obj = initMass - LTdetails(end);
end