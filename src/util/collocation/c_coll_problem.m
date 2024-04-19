classdef c_coll_problem
    %C_COLL_PROBLEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
       legs
       mesh
       node_method
       poly_order
       desInitNode
       desFinalNode
       isBallistic
       param %Ainv, B, D, W, nd_node_times, d_node_times, num_segs_per_traj, init_mass
       converged_opt_vars
    end
    
    methods
        function problem = c_coll_problem(traj_in,varargin)   
            % Ensure all inputs become a cell array
            if isa(traj_in,'c_traj')
                legsCell = {traj_in};
            else
                legsCell = traj_in;
            end
            checkClass(legsCell,'cell')

            % check for time continuity between trajectories
            for i = 1:length(legsCell)-1
                tfa = legsCell{i}.getFinalTime().value;
                t0b = legsCell{i+1}.getInitTime().value;
                if  tfa ~= t0b
                    error(['Trajectories must be continuous in time. traj #'...
                        num2str(i) ', ' legsCell{i}.name ', has final time: '...
                        num2str(tfa) '; traj #' numstr(i+1) ', ' ...
                        legsCell{i+1}.name ' has init time ' num2str(t0b)])
                end
            end
            
            % ---- Defaults ---- %
            nodeMethod = 'LGL';
            order = 7;
            allMesh = '5_equalTime';
            isProblemBallistic = false;
            desiredInitNode = legsCell{1}.getStateByIndex(1);
            desiredInitNode = desiredInitNode(1:6); %fully constrained init state
            desiredFinalNode = nan(6,1);
            % ------------------ %
            
            % User inputs
            for i = 1:length(varargin)/2
                switch varargin{i*2-1}
                    case 'mesh'
                        allMesh = varargin{i*2};
                    case 'nodeMethod'
                        nodeMethod = varargin{i*2};
                    case 'order'
                        order = varargin{i*2};
                    case 'isBallistic'
                        isProblemBallistic = varargin{i*2};
                    case 'desInitNode'
                        desiredInitNode = varargin{i*2};
                    case 'desFinalNode'
                        desiredFinalNode = varargin{i*2};
                    otherwise
                        error('unknown varargin')
                end
            end
    
            % Fix mesh to match shape of traj
            if ~iscell(allMesh)
                problem_mesh = cell(1,length(legsCell));
                for j = 1:length(legsCell)
                    problem_mesh{j} = allMesh;
                end
            elseif length(allMesh) == length(legsCell)
                problem_mesh = allMesh;
            else
                error('Length of mesh must equal length of traj')
            end
    
            % Ensure order of polynomial is odd and within reason
            if mod(order,2) == 0
                error('order must be odd');
            end
            if order > 9
                warning('Order greater than 9 likely to not yield better results. Try 7')
            end
            if order < 3
                error('Order must be 3 or greater')
            end

            %check low thrust control law
            if ~isProblemBallistic
                for i = 1:length(legsCell)
                    if ~strcmp(legsCell{i}.low_thrust.control_law.law,'fixed_dir')
                        error('Must used fixed direction control law')
                    end
                end
            end
            
            % Build problem
            problem.legs = legsCell;
            problem.mesh = problem_mesh;
            problem.poly_order = order;
            problem.node_method = nodeMethod;
            problem.desInitNode = desiredInitNode;
            problem.desFinalNode = desiredFinalNode;
            problem.isBallistic = isProblemBallistic;

            % define constant matrices A, B, D, W, and nd_node_times        
            oddInds = 1:2:order;
            evenInds = 2:2:order;
            switch nodeMethod
                case 'LGL'
                    [nodes,weights,~]=legendreGaussLobatto(order);
                    nodes = flip(nodes');
                    W = diag(weights(evenInds));
                otherwise
                    error('Have not implemented anything except LGL')
            end
            Apos = zeros(order+1,length(oddInds));
            Avel = Apos;
            for i = 0:order
                Apos(i+1,:) = nodes(oddInds).^i;
            end
            for i = 1:order
                Avel(i+1,:) = i*nodes(oddInds).^(i-1);
            end
            A = [Apos, Avel];
            Ainv = inv(A);
            % B with padding for tau = -1,1
            B = ones(order+1,length(evenInds)+2);
            B(2:2:end,1) = -1;
            for i = 0:order
                B(i+1,2:end-1) = nodes(evenInds).^i;
            end
            % D is derivative of B (minus padding)
            D = zeros(order+1,length(evenInds));
            for i = 1:order
                D(i+1,:) = i*nodes(evenInds).^(i-1);
            end  
        
            % assign constants
            p.Ainv = Ainv;
            p.B = B;
            p.D = D;
            p.W = W;
            p.nd_node_times = nodes;
            problem.param = p;
        end

        function [solved_problem, varargout] = solve(problem,varargin)
            checkClass(problem,'c_coll_problem')
            
            % ---- Default ---- %
            optMethod = 'fmincon'; %SNOPT, IPOPT
            meshRefineMethod = 'none'; %CEP, hybrid, de Boor
            options = optimoptions('fmincon','Display','iter');
            % ----------------- %

            if nargout > 1
                debugPlot = true;
            else
                debugPlot = false;
            end

            % User inputs
            for i = 1:length(varargin)/2
                switch varargin{i*2-1}
                    case 'optMethod'
                        optMethod = varargin{i*2};
                    case 'meshRefineMethod'
                        meshRefineMethod = varargin{i*2};
                    case 'options'
                        options = varargin{i*2};
                    otherwise
                        error('Unknown varargin')
                end
            end

            % build initial set of optimization variables
            [optVars, problem, lb, ub] = coll_opt_vars(problem);

            % plot initial guess if desired
            if debugPlot
                dps = defaultPlotting;
                dps.lineColor = colour('b');
                dps.lineWidth = 1;
                debugProblem = problem;
                debugProblem.converged_opt_vars = optVars; %not converged, for debugging only
                axDebug = debugProblem.plot('dps',dps,'DisplayName','Init Guess','plotPoly',false,'plotLT',false);
                varargout{1} = axDebug;
            end
            % define objective function
            objFunc = coll_obj_func(problem);
            
            % define constraint function
            conFunc = coll_con_func(problem);

            %% loop: solve optimization, mesh refinement
            % solve optimization
            switch optMethod
                case 'fmincon'
                    opt_soln = fmincon(objFunc,optVars,[],[],[],[],lb,ub,conFunc,options);
                otherwise
                    error('optMethod not yet implemented')
            end
            % mesh refinement
            switch meshRefineMethod
                case 'none'
                case 'de Boor'
                    error('Not implemented')
                case 'CEP'
                    error('Not implemented')
                case 'hybrid'
                    error('Not implemented')
            end

            %% assign output
            solved_problem = problem;
            solved_problem.converged_opt_vars = opt_soln;
        end

        function ax = plot(problem,varargin)
            % turn optimization variables into long list of variable nodes
            % and defect nodes.
            optVars = problem.converged_opt_vars;
            dps = [];
            trajName = 'Propagated Nodes';
            plotPoly = false;
            plotLT = ~problem.isBallistic;
            for i = 1:length(varargin)/2
                switch varargin{i*2-1}
                    case 'dps'
                        dps = varargin{i*2};
                    case 'ax'
                        ax = varargin{i*2};
                    case 'DisplayName'
                        trajName = varargin{i*2};
                    case 'plotPoly'
                        plotPoly = varargin{i*2};
                    case 'plotLT'
                        plotLT = varargin{i*2}; %only used for debug call
                    otherwise
                        error('unknown varargin')
                end
            end
            if plotPoly
                figure
                axVel = gca;
                hold on
                axis equal
                xlabel(axVel,'xVel')
                ylabel(axVel,'yVel')
                zlabel(axVel,'zVel')
            end
            if plotLT
                ffigure;
                subplot(3,1,1)
                axTdir = gca;
                subplot(3,1,2)
                axThr = gca;
                subplot(3,1,3)
                axMass = gca;
                alltdirHist = [];
                allTimes = [];
                throttleHist = [];
                massHist = c_dim_quant([],'kg');
            end
            p = problem.param;
            order = problem.poly_order;
            nd_node_times = p.nd_node_times;
            d_node_times = p.d_node_times;
            variableNodeInds = 1:2:order;
            defectNodeInds = 2:2:order;
            m0 = p.init_mass;
            tau = linspace(-1,1,1000);

            for i = 1:length(problem.legs)
                traj = problem.legs{i};
                traj.name = trajName;
                id_node_times = d_node_times{i};
                for k = 1:p.num_segs_per_traj(i)
                    tspan = id_node_times(k,[1 end]); %dt is per segment
                    d_dt = tspan(2) - tspan(1);
            
                    % draw out the optimization variables for this segment
                    if problem.isBallistic
                        LTvars = [];
                    else
                        LTvars = optVars(7:10); optVars = optVars([1:6,11:end]);
                    end
                    varStates = optVars(1:6*length(variableNodeInds));
                    optVars(1:6*(length(variableNodeInds)-1)) = []; %leave last node on there for the next segment
            
                    % get defect states, replace traj and m0 with mf
                    [variableStates,defectStates,~,traj,mf,Ci] = coll_def_nodes(LTvars,varStates,problem,traj,m0,d_dt);
                    allStates = [];
                    allStates(1:3,variableNodeInds) = variableStates(1:3,:);
                    allStates(4:6,variableNodeInds) = variableStates(4:6,:);
                    allStates(1:3,defectNodeInds) = defectStates(1:3,:);
                    allStates(4:6,defectNodeInds) = defectStates(4:6,:);
                    for j = 1:length(nd_node_times)-1 %plot each propagated node, grab histories
                        traj.pos.value = allStates(1:3,j);
                        traj.vel.value = allStates(4:6,j);
                        traj.time.value = [id_node_times(k,j), id_node_times(k,j+1)];
                        traj.low_thrust.mass.value = m0;
                        traj = traj.prop;
                        if i == 1 && j == 1 && k == 1
                            if ~isempty(dps)
                                traj.plotting = dps;
                            end
                            traj.plotting.HandleVisibility = 'on';
                        else
                            traj.plotting.HandleVisibility = 'off';
                        end
                        if exist('ax','var')==0
                            ax = traj.plot;
                        else
                            traj.plot('ax',ax);
                        end 
                        if plotPoly
                            if i == 1 && k == 1
                                plot3(axVel,traj.vel.value(1,:),traj.vel.value(2,:),traj.vel.value(3,:),'k-','DisplayName','Propagated Nodes')
                            else
                                plot3(axVel,traj.vel.value(1,:),traj.vel.value(2,:),traj.vel.value(3,:),'k-','HandleVisibility','off')
                            end
                        end
                        if plotLT
                            alltdirHist = [alltdirHist, traj.low_thrust.thrust_dir_history];
                            allTimes = [allTimes; traj.time.value];
                            throttleHist = [throttleHist; traj.low_thrust.spacecraft.throttle * ones(size(traj.time.value))];
                            massHist.value = [massHist.value; traj.low_thrust.mass.value'];
                        end
                        m0 = traj.low_thrust.mass.value(end);
                    end
                    m0 = mf;
                    
                    % plot polynomials
                    states = zeros(6,length(tau));
                    for j = 1:length(tau) %define states based on polynomial
                        jtau = tau(j);
                        jtauVec = zeros(order + 1,1);
                        for jj = 1:order + 1
                            jtauVec(jj) = jtau ^ (jj-1);
                        end
                        states(:,j) = Ci * jtauVec;
                    end
                    if plotPoly
                        plotStates(ax,states,'Color',colour('y'),'DisplayName','poly');
                        plot3(axVel,states(4,:),states(5,:),states(6,:)-0.01,'-','DisplayName','poly','Color',colour('y'))
                    end
                end
            end
            axis(ax,'equal')
            view(ax,0,90)
            if plotLT
                subplot(axTdir) %thrust direction
                    hold on
                    xlabel(['Time [$' traj.time.unit '$]'])
                    ylabel('Thrust Magnitude')
                    LTmag = vecnorm(alltdirHist);
                    plot(allTimes,LTmag,'k-','DisplayName','Thrust Magnitude')
                    plot(allTimes,alltdirHist(1,:),'Color',colour('r'),'DisplayName','$T_x$')
                    plot(allTimes,alltdirHist(2,:),'Color',colour('g'),'DisplayName','$T_y$')
                    plot(allTimes,alltdirHist(3,:),'Color',colour('b'),'DisplayName','$T_z$')
                    title('Thrust Direction')
                    legend
                subplot(axThr) %thrust magnitude
                    hold on
                    xlabel(['Time [$' traj.time.unit '$]'])
                    ylabel('Throttle')
                    plot(allTimes,throttleHist,'k-')
                    title('Throttle History')
                subplot(axMass) %mass history
                    sysModel = traj.system_model;
                    sysModel.char.mstar = traj.low_thrust.spacecraft.M0;
                    massHist = massHist.change_unit('kg', sysModel);
                    hold on
                    xlabel(['Time [$' traj.time.unit '$]'])
                    ylabel(['Mass [' massHist.unit ']'])
                    plot(allTimes,massHist.value,'k-')
                    title('Mass History')
            end
        end
    end
end