classdef c_coll_problem_v0
    %C_COLL_PROBLEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        traj
        mesh
        node_method
        order
        Ainv
        B
        D
        W
        nodes
        tol_mlp
        tol_err_dist
        tol_poly
        isBallistic
        desInitNode
        desFinalNode
    end
    
    methods
        function problem = c_coll_problem_v0(traj,varargin)
        % create A matrix of normalized times defined by mesh for polynomial of given order
        % define nodes based on nodeMethod
        % also defines weights for quadrature
        
        if iscell(traj)
            temptraj = traj{1};
        else
            temptraj = traj;
        end
        if ~strcmp(temptraj.low_thrust.control_law.law,'fixed_dir')
            error('Must used fixed direction control law')
        end
        
        % ---- Defaults ---- %
        nodeMethod = 'LGL';
        order = 7;
        mesh = '5_equalTime';
        if iscell(traj)
            mesh = cell(1,length(traj));
            for i = 1:length(traj)
                mesh{i} = '5_equalTime';
            end
        end
        Bpad = false;
        isBallistic = false;
        desInitNode = temptraj.getStateByIndex(1);
        desInitNode = desInitNode(1:6);
        desFinalNode = nan(6,1);
        % ------------------ %
        
        % User inputs
        for i = 1:length(varargin)/2
            switch varargin{i*2-1}
                case 'mesh'
                    mesh = varargin{i*2};
                    if iscell(traj) && ~iscell(mesh)
                        allMesh = mesh;
                        mesh = cell(1,length(traj));
                        for j = 1:length(traj)
                            mesh{j} = allMesh;
                        end
                    end
                case 'nodeMethod'
                    nodeMethod = varargin{i*2};
                case 'order'
                    order = varargin{i*2};
                case 'Bpad'
                    Bpad = varargin{i*2};
                case 'isBallistic'
                    isBallistic = varargin{i*2};
                case 'desInitNode'
                    desInitNode = varargin{i*2};
                case 'desFinalNode'
                    desFinalNode = varargin{i*2};
                otherwise
                    error('unknown varargin')
            end
        end
        if mod(order,2) == 0
            error('order must be odd');
        end
        if order > 9
            warning('Order greater than 9 likely to not yield better results. Try 7')
        end
        if order < 3
            error('Order must be 3 or greater')
        end
        oddInds = 1:2:order;
        evenInds = 2:2:order;
        
        %% first step: define A,B,D in each segment of the mesh
        if isnumeric(nodeMethod)
            % user input indices
            error('Have not implemented user input nodes')
        else
            switch nodeMethod
                case 'LGL'
                    [nodes,weights,~]=legendreGaussLobatto(order);
                    nodes = flip(nodes');
                    W = diag(weights(evenInds));
                otherwise
                    error('Have not implemented anything except LGL')
            end
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
        if Bpad
            % B with padding for tau = -1,1
            B = ones(order+1,length(evenInds)+2);
            B(2:2:end,1) = -1;
            for i = 0:order
                B(i+1,2:end-1) = nodes(evenInds).^i;
            end
        else
            % B without padding
            B = ones(order+1,length(evenInds));
            for i = 0:order
                B(i+1,:) = nodes(evenInds).^i;
            end
        end
        D = zeros(order+1,length(evenInds));
        for i = 1:order
            D(i+1,:) = i*nodes(evenInds).^(i-1);
        end        
        
        %% build problem
        problem.traj = traj;
        problem.mesh = mesh;
        problem.order = order;
        problem.node_method = nodeMethod;
        problem.nodes = nodes;
        problem.Ainv = Ainv;
        problem.B = B;
        problem.D = D;
        problem.W = W;
        problem.isBallistic = isBallistic;
        problem.desInitNode = desInitNode;
        problem.desFinalNode = desFinalNode;
        end
        
        function [collSoln, varargout] = solve(problem,varargin)
            % recursive method: calls problem.solve until all tolerances met
            checkClass(problem,'c_coll_problem')
            
            % ---- Default ---- %
            optMethod = 'fmincon'; %SNOPT, IPOPT
            meshRefineMethod = 'de Boor'; %CEP, hybrid
            options = optimoptions('fmincon','Display','iter');
            % ----------------- %

            %% User inputs
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
            if nargout == 2
                debugPlot = true;
            else
                debugPlot = false;
            end
            
            %% solve this problem
            % get states and derivatives at each node in each segment
            % Define segment and node times
            traj_ = problem.traj;
            mesh_ = problem.mesh;
            nodes_ = problem.nodes;
            nodeStates = [];
            nodeTimes = [];
            LTdetails = [];
            if iscell(traj_)
                trajLength = length(traj_);
                for i = 1:trajLength
                    [inodeStates, iNodeTimes, iLTdetails,meshSegs] = getNodes(traj_{i},mesh_{i},nodes_,problem.order);
                    if i > 1 %remove first node since it uses previous traj
                        inodeStates = inodeStates(7:end); 
                        iLTdetails = iLTdetails(2:end);
                    end
                    %each traj: 6*order + 6*(order-1)*(numSeg-1)
                    nodeStates = [nodeStates; inodeStates]; 
                    nodeTimes = [nodeTimes; iNodeTimes];
                    LTdetails = [LTdetails; iLTdetails]; %each segment uses 4 values
                end
            elseif isa(traj_,'c_traj')
                trajLength = 1;
                [nodeStates, nodeTimes, LTdetails,meshSegs] = getNodes(traj_,mesh_,nodes_,problem.order);
            else
                error('traj must either be a c_traj or cell array of c_traj')
            end
            %yields initial guess for nodeStates
            initMass = LTdetails(1);
            LTdetails = LTdetails(2:end);
            LTlb = [-1;-1;-1;0]-eps; %uhat, mass_f
            LTub = [1;1;1;initMass]+eps;
            LTlb = repmat(LTlb,[length(LTdetails)/4,1]);
            LTub = repmat(LTub,[length(LTdetails)/4,1]);
            Sub = inf(size(nodeStates));
            Slb = -inf(size(nodeStates));
            if ~problem.isBallistic
                optVars = [LTdetails; nodeStates]; %variables for optimizer
                ub = [LTub; Sub];
                lb = [LTlb; Slb];
            else
                optVars = nodeStates;
                ub = Sub;
                lb = Slb;
            end
            if debugPlot
                dps = defaultPlotting;
                dps.lineColor = colour('b');
                dps.lineWidth = 1;
                tempSoln = struct('soln',optVars,'nodeTimes',nodeTimes,...
                    'meshSegs',meshSegs,'initMass',initMass); %meshSegs = number of segments per trajectory due to mesh
                ax = problem.plot(tempSoln,'pkm',dps,'DisplayName','Init Guess','plotPoly',false);
                varargout{1} = ax;
            end
            % calculate C_i, build X and F
            collCons = @(x) collConFunc(x,problem,meshSegs,nodeTimes,initMass);
            if problem.isBallistic
                objFun = @(x) ballisticObjFunc(x);
            else
                objFun = @(x) collObjFunc(x,initMass,meshSegs*trajLength);
%                 objFun = @(x) ballisticObjFunc(x);
            end
            switch optMethod
                case 'fmincon'
                    soln = fmincon(objFun,optVars,[],[],[],[],lb,ub,collCons,options);
                otherwise
                    error('optMethod not yet implemented')
            end

            %% check for mesh refinement
            

            %% output soln
            collSoln = struct('soln',soln,'nodeTimes',nodeTimes,'meshSegs',meshSegs,'initMass',initMass);
        end

        function ax = plot(problem,collSoln,varargin)
            dps = defaultPlotting;
            trajName = 'Propagated Nodes';
            plotPoly = false;
            plotLT = false;
            for i = 1:length(varargin)/2
                switch varargin{i*2-1}
                    case 'pkm'
                        dps = varargin{i*2};
                    case 'ax'
                        ax = varargin{i*2};
                    case 'DisplayName'
                        trajName = varargin{i*2};
                    case 'plotPoly'
                        plotPoly = varargin{i*2};
                    case 'plotLT'
                        plotLT = varargin{i*2};
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
            dpsOff = dps; dpsOff.HandleVisibility = 'off';
            dpsBad = dps; dpsBad.LineColor = colour('r'); dpsBad.HandleVisibility = 'off';
            optVars = collSoln.soln;
            nodeTimes = collSoln.nodeTimes;
            meshSegs = collSoln.meshSegs;
            initMass = collSoln.initMass;
            order_ = problem.order;
            traj_ = problem.traj;
            nodes_ = problem.nodes;
            if iscell(traj_)
                trajLength = length(traj_);
                traj_ = traj_{1};
            else
                trajLength = 1;
            end

            if plotLT
                ffigure;
                subplot(1,3,1)
                axTdir = gca;
                subplot(1,3,2)
                axThr = gca;
                subplot(1,3,3)
                axMass = gca;
            end
            traj_.name = trajName;
            traj_.low_thrust.control_law = c_control_law(traj_.low_thrust.control_law.law,...
                            traj_.low_thrust.control_law.frame,[1;0;0]); %temp definition of direction
            totalSegs = trajLength * meshSegs;
            if problem.isBallistic
                optVars = [zeros(4*totalSegs,1); optVars];
            end
            LTdetails = optVars(1:4*totalSegs);
            nodeStates = optVars(4*totalSegs+1:end);
            Ainv_ = problem.Ainv;
            B_ = problem.B;
            variableNodes = 1:2:order_;
            defectNodes = 2:2:order_;
            alltdirHist = [];
            allTimes = [];
            throttleHist = [];
            massHist = traj_.low_thrust.mass; massHist.value = [];
            for i = 1:totalSegs %number of segments
                tspan = nodeTimes(i,[1 end]); %dt is per segment
                dt = tspan(2) - tspan(1);
                isBad = false;
                if problem.isBallistic
                    traj_.low_thrust.spacecraft.throttle = 0;
                    trajMass = 1;
                else
                    %define low thrust ----
                    iLTdetails = LTdetails(1:4);
                    LTdetails(1:4) = []; %pop off these values
                    traj_.low_thrust.control_law.coeffs = iLTdetails(1:3)./norm(iLTdetails(1:3));
                    if i == 1
                        m_i = initMass;
                    end
                    m_f = iLTdetails(4);
                    traj_.low_thrust.mass.value = [m_i, m_f];
                    throttle = (m_i - m_f) / (traj_.low_thrust.spacecraft.mdotND*dt);
                    traj_.low_thrust.spacecraft.throttle = throttle; %define throttle based on mass variables
                    if throttle > 1
                        isBad = true;
                    end
                    trajMass = m_i;
                    m_i = m_f; %for next segment
                end
                %define states ----
                variableStates = zeros(6,length(variableNodes));
                variableDerivs = zeros(6,length(variableNodes));
                for j = 1:length(variableNodes) %build variableStates and variableDerivs
                    traj_.pos.value = nodeStates(1:3);
                    traj_.vel.value = nodeStates(4:6);
                    if j ~= length(variableNodes)
                        nodeStates(1:6) = []; %pop off this node, unless using for next segment
                    end
                    variableStates(1:3,j) = traj_.pos.value;
                    variableStates(4:6,j) = traj_.vel.value;
                    traj_.time.value = nodes_(variableNodes(j));
                    variableDerivs(1:3,j) = traj_.vel.value;
                    variableDerivs(4:6,j) = traj_.EOM();
                end % yields 9xorder matrix of nodes along a segment
                variableDerivs = variableDerivs .* dt/2; %normalize derivatives wrt time
                Ci = [variableStates, variableDerivs] * Ainv_;
                defectStates = Ci*B_; %why isn't this the right value???
                allStates(1:3,variableNodes) = variableStates(1:3,:);
                allStates(4:6,variableNodes) = variableStates(4:6,:);
                allStates(1:3,defectNodes) = defectStates(1:3,:);
                allStates(4:6,defectNodes) = defectStates(4:6,:);
                for j = 1:length(nodes_)-1 %plot each propagated node, grab histories
                    traj_.pos.value = allStates(1:3,j);
                    traj_.vel.value = allStates(4:6,j);
                    traj_.time.value = [nodeTimes(i,j), nodeTimes(i,j+1)];
                    traj_.low_thrust.mass.value = trajMass;
                    traj_ = traj_.prop;
                    if i == 1 && j == 1
                        traj_.plotting = dps;
                    elseif isBad
                        traj_.plotting = dpsBad;
                    else
                        traj_.plotting = dpsOff;
                    end
                    if exist('ax','var')==0
                        ax = traj_.plot;
                    else
                        traj_.plot('ax',ax);
                    end 
                    if plotPoly
                        if j == 1
                            plot3(axVel,traj_.vel.value(1,:),traj_.vel.value(2,:),traj_.vel.value(3,:),'k-','DisplayName','Propagated Nodes')
                        else
                            plot3(axVel,traj_.vel.value(1,:),traj_.vel.value(2,:),traj_.vel.value(3,:),'k-','HandleVisibility','off')
                        end
                    end
                    if plotLT
                        alltdirHist = [alltdirHist, traj_.low_thrust.thrust_dir_history];
                        allTimes = [allTimes; traj_.time.value];
                        throttleHist = [throttleHist; throttle * ones(size(traj_.time.value))];
                        massHist.value = [massHist.value; traj_.low_thrust.mass.value'];
                    end
                    trajMass = traj_.low_thrust.mass.value(end);
                end
                
                % plot polynomials
                tau = linspace(-1,1,1000);
                states = zeros(6,length(tau));
                for j = 1:length(tau) %define states based on polynomial
                    jtau = tau(j);
                    jtauVec = zeros(order_ + 1,1);
                    for k = 1:order_+1
                        jtauVec(k) = jtau ^ (k-1);
                    end
                    states(:,j) = Ci * jtauVec;
                end
                if plotPoly
                    plotStates(ax,states,'Color',colour('y'),'DisplayName','poly');
                    plot3(axVel,states(4,:),states(5,:),states(6,:)-0.01,'-','DisplayName','poly','Color',colour('y'))
                end
            end
            if plotLT
                subplot(axTdir) %thrust direction
                    hold on
                    xlabel(['Time [$' traj_.time.unit '$]'])
                    ylabel(['Thrust Magnitude'])
                    LTmag = vecnorm(alltdirHist);
                    plot(allTimes,LTmag,'k-','DisplayName','Thrust Magnitude')
                    plot(allTimes,alltdirHist(1,:),'Color',colour('r'),'DisplayName','$T_x$')
                    plot(allTimes,alltdirHist(2,:),'Color',colour('g'),'DisplayName','$T_y$')
                    plot(allTimes,alltdirHist(3,:),'Color',colour('b'),'DisplayName','$T_z$')
                    title('Thrust Direction')
                    legend
                subplot(axThr) %thrust magnitude
                    hold on
                    xlabel(['Time [$' traj_.time.unit '$]'])
                    ylabel('Throttle')
                    plot(allTimes,throttleHist,'k-')
                    title('Throttle History')
                subplot(axMass) %mass history
                    sysModel = traj_.system_model;
                    sysModel.char.mstar = traj_.low_thrust.spacecraft.M0;
                    massHist = massHist.change_unit('kg', sysModel);
                    hold on
                    xlabel(['Time [$' traj_.time.unit '$]'])
                    ylabel(['Mass [' massHist.unit ']'])
                    plot(allTimes,massHist.value,'k-')
                    title(['Mass History'])
            end
        end
    end
end

%% objective function
function obj = collObjFunc(optVars,initMass,totalSegs)
    LTdetails = optVars(1:4*totalSegs); %column vector for all segments
    obj = 0;
    for i = 1:totalSegs %number of segments
        m_f = LTdetails(4);
        obj = obj + (initMass - m_f);
        initMass = m_f;
        LTdetails(1:4) = [];
    end
%     obj = initMass - LTdetails(end);
end

function obj = ballisticObjFunc(optVars)
    obj = 1;
end

%% Constraint function
function [ineqCons, eqCons] = collConFunc(optVars,problem,meshSegs,nodeTimes,initMass)
    order = problem.order;
    isBallistic = problem.isBallistic;
    desInitNode = problem.desInitNode;
    desFinalNode = problem.desFinalNode;
    Ainv = problem.Ainv;
    B = problem.B;
    D = problem.D;
    W = problem.W;
    nodes = problem.nodes;
    traj = problem.traj;
    if iscell(traj)
        trajLength = length(traj);
        traj = traj{1};
    else
        trajLength = 1;
    end
    traj.low_thrust.control_law = c_control_law(traj.low_thrust.control_law.law,...
                            traj.low_thrust.control_law.frame,[1;0;0]); %temp direction
    totalSegs = trajLength * meshSegs;
    if isBallistic
        optVars = [zeros(4*totalSegs,1); optVars];
    end
    LTdetails = optVars(1:4*totalSegs);
    nodeStates = optVars(4*totalSegs+1:end);
    defectEqns = [];
    uhatEqns = [];
    initNode = nodeStates(1:6);
    thrEqns = [];
    for i = 1:totalSegs %number of segments
        tspan = nodeTimes(i,[1 end]); %dt is per segment
        dt = tspan(2) - tspan(1);
        if isBallistic
            traj.low_thrust.spacecraft.throttle = 0;
        else
            %define low thrust
            iLTdetails = LTdetails(1:4);
            LTdetails(1:4) = []; %pop off these values
            traj.low_thrust.control_law.coeffs = iLTdetails(1:3)./norm(iLTdetails(1:3)); %normalize so it doesn't error
            uhatEqns = [uhatEqns; norm(iLTdetails(1:3)) - 1];
            if i == 1
                m_i = initMass;
            end
            m_f = iLTdetails(4);
            traj.low_thrust.mass.value = [m_i, m_f];
            throttle = (m_i - m_f) / (traj.low_thrust.spacecraft.mdotND.value*dt);
            traj.low_thrust.spacecraft.throttle = throttle; %define throttle based on mass variables
            thrEqns = [thrEqns; throttle - 1; -throttle]; % c <= 0
            m_i = m_f; %for next segment
        end
        %define states ----
        variableNodes = 1:2:order;
        defectNodes = 2:2:order;
        variableStates = zeros(6,length(variableNodes));
        variableDerivs = zeros(6,length(variableNodes));
        for j = 1:length(variableNodes) %build variableStates and variableDerivs
            traj.pos.value = nodeStates(1:3);
            traj.vel.value = nodeStates(4:6);
            if j ~= length(variableNodes)
                nodeStates(1:6) = []; %pop off this node, unless using for next segment
            end
            variableStates(1:3,j) = traj.pos.value;
            variableStates(4:6,j) = traj.vel.value;
            traj.time.value = nodes(variableNodes(j));
            variableDerivs(1:3,j) = traj.vel.value;
            variableDerivs(4:6,j) = traj.EOM();
        end % yields 9xorder matrix of nodes along a segment
        variableDerivs = variableDerivs .* dt/2; %normalize derivatives wrt time
        Ci = [variableStates, variableDerivs] * Ainv;
        defectStates = Ci*B;
        defectDerivs = Ci*D;
        trueDefectDerivs = zeros(6,length(defectNodes));
        for j = 1:length(defectNodes)
            trueDefectDerivs(1:3,j) = defectStates(4:6,j);
            traj.pos.value = defectStates(1:3,j);
            traj.vel.value = defectStates(4:6,j);
            traj.time.value = nodes(defectNodes(j));
            trueDefectDerivs(4:6,j) = traj.EOM();
            trueDefectDerivs(1:6,j) = trueDefectDerivs(1:6,j)  * dt/2;
        end
        defectEqns = [defectEqns; (defectDerivs - trueDefectDerivs)*W];
    end
    % fix beginning and end points
    finalStateEqn = desFinalNode - nodeStates; %all that is left should be a 6x1 vector
    initStateEqn = desInitNode - initNode;
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

%% obtain nodes from mesh
function [nodeStates,nodeTimes,LTdetails,meshSegs] = getNodes(traj,mesh,nodes,order)
    %mesh is either a set of times along traj, or a string defining how to
    %build the mesh
    %nodes are normalized time [-1,1]
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
    nodeTimes = zeros(meshSegs,length(nodes));
    % propagate traj for each node within segment timespan to get exact
    % states (only propagates trajectory once)
    nodeTraj = traj;
    LTdetails = traj.getInitMass().value;
    nodeStates = [nodeTraj.getInitPos().value; nodeTraj.getInitVel().value];
    timeunit = traj.time.unit;
    segTimes = segTimes(:); %column vector
    defectNodes = 2:2:order;
    for i = 1:meshSegs
        t1 = segTimes(i);
        t2 = segTimes(i+1);
        dimNodes = t1 + (t2-t1)*(nodes+1)/2;
        nodeTimes(i,:) = dimNodes;
        inodeStates = [];
        iLT(1:3) = nodeTraj.low_thrust.control_law.coeffs;
        for j = 1:length(nodes)-1
            nodeSpan = dimNodes(j:j+1);
            nodeTraj.time = c_dim_quant(nodeSpan,timeunit);
            nodeTraj = nodeTraj.prop;
            nodeTraj.pos = nodeTraj.getFinalPos;
            nodeTraj.vel = nodeTraj.getFinalVel;
            if any(j==defectNodes) %only append if a variable node
                inodeStates = [inodeStates; nodeTraj.pos.value];
                inodeStates = [inodeStates; nodeTraj.vel.value];
            end
            nodeTraj.low_thrust.mass = nodeTraj.getFinalMass(); %start with new mass
            if j == length(nodes)-1
                iLT(4) = nodeTraj.getFinalMass().value;
            end
        end
        LTdetails = [LTdetails; iLT(:)];
        nodeStates = [nodeStates; inodeStates];
    end
end