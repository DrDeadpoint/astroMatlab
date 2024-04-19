classdef c_traj
    % obj = c_traj(name, time, pos, vel, sysModel, varargin)
    %
    % time is 1xn c_dim_quant of times matching the states (or a start and end time)
    % pos is 3xn c_dim_quant of positions matching the times
    % vel is 3xn c_dim_quant of velocities matching the times
    % sysModel is a c_system_model
    % varargin used to add extra details:
    %   (...,'STM',STM,...) adds an STM for tf in time. This can be calculated by simply propagating the trajectory
    %   (...,'mass',mass,...) adds a vector of masses for the spacecraft at times in time_vec (LT)
    %   (...,'spacecraft',spacecraft,...) adds a structure with thrust level details for the spacecraft (LT)
    %   (...,'control_law',claw,...) adds a low-thrust control law from the class CONTROL_LAW
    %   (...,'thrust_history',thrust_direction_history,...) adds a vector of thrust directions associated with the state vector
    %   (...,'plotting',plotVars,...) adds a structure containing details for when the trajectory is plotted
    %       default plotVars:   .lineColor = colour('k')
    %                           .lineWidth = 2
    %                           .lineStyle = '-'
    %                           .HandleVisibility = 'on'
    %                           
    % Alex Hoffman
    
    properties
        name char = 'traj'
        system_model c_system_model
        time c_dim_quant
        pos c_dim_quant
        vel c_dim_quant
        stm c_stm
        low_thrust struct
        plotting struct
        etc
    end
    
    methods
        function obj = c_traj(name, time, pos, vel, sysModel, varargin)
            % obj = c_traj(name, time, pos, vel, sysModel, varargin)
            %
            % time is 1xn c_dim_quant of times matching the states (or a start and end time)
            % pos is 3xn c_dim_quant of positions matching the times
            % vel is 3xn c_dim_quant of velocities matching the times
            % sysModel is a c_system_model
            % varargin used to add extra details:
            %   (...,'STM',STM,...) adds an STM for tf in time. This can be calculated by simply propagating the trajectory
            %   (...,'mass',mass,...) adds a vector of masses for the spacecraft at times in time_vec (LT)
            %   (...,'spacecraft',spacecraft,...) adds a structure with thrust level details for the spacecraft (LT)
            %   (...,'control_law',claw,...) adds a low-thrust control law from the class CONTROL_LAW
            %   (...,'thrust_history',thrust_direction_history,...) adds a vector of thrust directions associated with the state vector
            %   (...,'plotting',plotVars,...) adds a structure containing details for when the trajectory is plotted
            %       default plotVars:   .lineColor = colour('k')
            %                           .lineWidth = 2
            %                           .lineStyle = '-'
            %                           .HandleVisibility = 'on'
            %                           
            % Alex Hoffman
            obj.name = name;
            if size(pos.value,1) == 3
                obj.pos = pos.change_unit('nd_l',sysModel);
            else
                error('Position must have dim 3xn')
            end
            if size(vel.value,1) == 3
                obj.vel = vel.change_unit('nd_v',sysModel);
            else
                error('Velocity must have dim 3xn')
            end
            obj.time = time.change_unit('nd_t',sysModel);
            obj.system_model = sysModel;
            %default low-thrust details (defaults to no thrust)
            defMass = c_dim_quant(1, 'kg');
            defTmax = c_dim_quant(1, 'N');
            defIsp = c_dim_quant(1, 'sec');
            defM0 = c_dim_quant(1, 'kg');
            defSpacecraft = c_spacecraft(defTmax,defIsp,defM0,sysModel,0);
            defControlLaw = c_control_law('SOC',sysModel.frame,struct());
            defThrustHistory = [];
            low_thrust = struct('mass',defMass,'spacecraft',defSpacecraft,...
                'control_law',defControlLaw,'thrust_dir_history',defThrustHistory);
            obj.low_thrust = low_thrust;
            %default plotting details
            obj.plotting = defaultPlotting;
            for i = 1:length(varargin)/2
               switch varargin{i*2-1}
                   case 'STM'
                       obj.stm = varargin{i*2};
                   case 'mass'
                       obj.low_thrust.mass = varargin{i*2};
                   case 'spacecraft'
                       obj.low_thrust.spacecraft = varargin{i*2};
                   case 'control_law'
                       obj.low_thrust.control_law = varargin{i*2};
                   case 'thrust_history'
                       obj.low_thrust.thrust_dir_history = varargin{i*2};
                   case 'plotting'
                       obj.plotting = varargin{i*2};
                   case 'ephemSeg'
                       obj.etc.seg = varargin{i*2};
                   case 'ephemSystem'
                       obj.etc.system = varargin{i*2};
                   case 'ephemOptions'
                       obj.etc.ephemOptions = varargin{i*2};
                   otherwise
                       error('Unknown varargin')
               end
            end
        end
        
        %% Equations of Motion
        function [acc] = EOM(obj,varargin)
            if nargin == 1
                index = 1;
            elseif nargin == 2
                index = varargin{1};
            end
            pos_ = obj.getPosByIndex(index);
            vel_ = obj.getVelByIndex(index);
            mass_ = obj.getMassByIndex(index);
            t = obj.getTimeByIndex(index);
            t0 = obj.getTimeByIndex(1);
            state = [pos_.value; vel_.value; mass_.value];
            model = obj.system_model;
            Tmax = obj.low_thrust.spacecraft.TmaxND.value * ...
                obj.low_thrust.spacecraft.throttle; %scale by throttle
            claw = obj.low_thrust.control_law;
            switch model.system_dynamics
                case '2BP'
                    [xdd, ydd, zdd] = eom_2bp(state,...
                        model.char.mu.value,model.char.J2,...
                        model.P1.radius.value,t.value,...
                        t0.value,Tmax,claw,false);
                    acc = c_dim_quant([xdd;ydd;zdd],'nd_a');
                case 'CR3BP'
                    [xdd, ydd, zdd] = eom_cr3bp(state,...
                        model.char.mu.value,t.value,...
                        t0.value,Tmax,claw,false);
                    acc = c_dim_quant([xdd;ydd;zdd],'nd_a');
                case 'BCR4BP'
                    char = model.char;
                    nP4 = char.angVelP4.value+1;
                    [xdd, ydd, zdd] = eom_bcr4bp_p1p2(state,...
                        char.mu.value,t.value,t0.value,char.massP4.value,...
                        char.distP4.value,char.descNodeP4.value,...
                        char.incP4.value,char.theta0P4.value,...
                        nP4,Tmax,claw,false);
                    acc = c_dim_quant([xdd;ydd;zdd],'nd_a');
                case 'EPHEMERIS'
                    ephem_options = obj.etc.ephemOptions;
                    % seg = ephem_options{1};
                    seg = obj.etc.seg;
                    segNo = seg.segNumber;
                    if length(ephem_options) > 5
                        harm = ephem_options{6};
                    else
                        harm = {false;[]}; %don't use harmonics
                    end
                    
                    sysModel = obj.system_model;
                    char = sysModel.char;
                    EMFrame = char.EMFrame;
                    SEFrame = char.SEFrame;
                    EMFrame2 = char.EMFrame2;
                    Body = char.Body;
                    JDFix = char.JDFix;
                    
                    obj.pos = obj.pos.change_unit('nd_l',sysModel);
                    obj.vel = obj.vel.change_unit('nd_v',sysModel);
                    Y0 = obj.getStateByIndex(index)';

                    obj.time = obj.time.change_unit('nd_t',sysModel);
                    t = obj.getTimeByIndex(index);
                    t0 = obj.getTimeByIndex(1);

                    switch seg.prop.centralBody
                        case 'EARTH'
                            propFrame = EMFrame;
                        case 'MOON'
                            propFrame = EMFrame2;
                        case 'SUN'
                            propFrame = SEFrame;
                            error('SE frame not implemented');
                        otherwise
                            error(['Unknown propagation frame on segment ' segNo])
                    end
                    
                    if ~strcmp(seg.IC.stateVars{1},'Rx (km)')
                       error(['Segment ' num2str(ii) ...
                           ' does not have state in position and velocity space'])
                    end
                    
                    % propagation
                    m0seg1 = obj.low_thrust.spacecraft.M0; %kg
                    m0 = obj.low_thrust.mass;
                    use_lowthrust = false;
                    use_harmonics = harm{1};
                    %initialize
                    uCoef = [];
                    SC = [];
                    manFrame = [];
                    manSys = [];
                    harmonics_coefs = harm{2};
                    %declare values
                    if ~isnan(seg.maneuver.frame) % a maneuver takes place
                        use_lowthrust = true;
                        claw = obj.low_thrust.control_law;
                        spacecraft = obj.low_thrust.spacecraft;
                        uCoef = [claw.coeffs.alpha0; claw.coeffs.alphadot;...
                            claw.coeffs.beta0; claw.coeffs.betadot];
                        thrust = spacecraft.Tmax.change_unit('N');
                        thrust = thrust.value; 
                        thrust = thrust * spacecraft.throttle;
                        isp = spacecraft.Isp.value;
                        SC.m0D = m0seg1.value;
                        SC.thrustMaxD = thrust;
                        SC.ispD = isp;
                        SC = setSpacecraft(obj.etc.system,SC);
                        SC = SC{1};
                        % sometimes maneuver is in J2000 frame, sometimes in 2B_rotating
                        switch claw.frame{1}
                            case 'J2000'
                                switch claw.frame{2}
                                    case 'IJK'
                                        manFrame = 'J2000';
                                    case 'VUW'
                                        if strcmp(seg.maneuver.centralBody,propFrame.centralBody)
                                            manFrame = 'J2000 VUW';
                                        else
                                            error('need new prop code')
                                        end
                                end
                            case '2B_ROTATING'
                                manFrame = '2B rotating VUW';
                                manMainBody = seg.maneuver.mainBody;
                                manCentralBody = claw.frame{3};
                                switch manMainBody
                                    case 'SUN'
                                        manSys = SEFrame;
                                    case 'EARTH'
                                        switch manCentralBody
                                            case 'EARTH'
                                                manSys = EMFrame;
                                            case 'MOON'
                                                manSys = EMFrame2;
                                        end
                                end
                        end
                    end
                    dYdt = ephemeris(t, Y0, JDFix, propFrame, Body,...
                        'lowthrust',use_lowthrust,'uCoef',uCoef,'t0',t0,'spacecraft',SC,...
                        'manFrame',manFrame,'manSystem',manSys,'harmonics',use_harmonics,...
                        'harmonics_coefs',harmonics_coefs);
                    acc = c_dim_quant(dYdt(4:6),'nd_a');
                otherwise
                    error('model not implemented')
            end
        end

        %% Propagate
        function [traj_out,myvarargout] = prop(traj_in,varargin)
            %   [traj_out,varargout] = PROP(traj_in,varargin)
            %   traj_in contains an initial time, final time (both in .time), and initial state
            %   varargin:
            %    for all:
            %       (...,'options',options_in,...)
            %       (...,'prop',prop,...) where prop is either 'atd' or 'manual'
            %       (...,'events',events,...) where events could be 'xzcrossing' or 'perilune' (NOT IMPLEMENTED FOR ALL PROPAGATORS)
            %    for n-body only:
            %       (...,'mVec',[m1 m2 ... mn],...)
            %
            %   tout will be n-vector of times for calculated stateVec
            %   stateVec will be 6xn matrix of states for each time in tout
            %   varargout: 
            %       {1} solution struct for event propagation
            %       {2} for bcfbp SB1 only: varargout{2} is the positions of earth and moon
            %
            % Alex Hoffman
            % Date: 11/30/2020
            % Last Revision: 10/22/2021
            prop = 'manual';
            options = odeset('RelTol',1e-12,'AbsTol',1e-12);        
            events = 'none';
            mVec = [];
            wantSTM = false;
            useKepler = false;
            modProp = false;
            for i = 1:length(varargin)/2
               switch varargin{i*2-1}
                   case 'options'
                       options = varargin{i*2};
                   case 'events'
                       events = varargin{i*2};
                   case 'prop'
                       prop = lower(varargin{i*2});
                   case 'mVec'
                       mVec = varargin{i*2};
                   case 'wantSTM'
                       wantSTM = varargin{i*2};
                   case 'useKepler'
                       useKepler = varargin{i*2};
                   case 'modProp'
                       modProp = varargin{i*2};
                   otherwise
                       error('Unknown varargin')
               end
            end
            % call a propagator
            myvarargout{1} = NaN; %{1}
            myvarargout{2} = NaN;
            sysModel = traj_in.system_model;
            switch sysModel.system_dynamics
                case 'CR3BP'
                    switch sysModel.frame
                        case 'B1centP1P2rot'
                            [traj_out,myvarargout] = prop_cr3bp(traj_in,options,prop);
                        otherwise
                            error('Need to implement frame change here')
                    end
                case 'BCR4BP'
                    switch sysModel.frame
                        case 'B1centP1P2rot'
                            [traj_out] = prop_bcr4bp_p1p2(traj_in,options,prop);
                        case 'B2centP4B1rot'
                            [traj_out] = prop_bcr4bp_p4b1(traj_in,options,prop);
%                             myvarargout{2} = {P1vec, P2vec};
                        otherwise
                            error('bad frame')
                    end
                case 'EPHEMERIS'
                    traj_out = prop_ephem(traj_in,wantSTM);
                case 'NBODY'
                    %[traj_out] = propNBODY(traj_in,options,mVec,prop);
                case '2BP'
                    if useKepler
                        traj_out = keplerProblem(traj_in);
                    elseif modProp
                        sysModel = traj_in.system_model;
                        tspan = traj_in.getTimeSpan().value;
                        r0vec = traj_in.getInitPos().value;
                        v0vec = traj_in.getInitVel().value;
                        mu = sysModel.char.mu.value;
                        
                        dt = tspan(2) - tspan(1);
                        r0 = norm(r0vec);
                        v0 = norm(v0vec);
                        a = 1/(2/r0-v0^2/mu); %sma, nd
                        per = 2*pi*sqrt(a^3/mu); %period, nd
                        dt = sign(dt) * mod(abs(dt),per);
                        traj_in.time = c_dim_quant([tspan(1), tspan(1) + dt],'nd_t');
                        oblate = [sysModel.P1.radius.value, sysModel.char.J2];
                        [traj_out] = prop_2bp(traj_in,options,oblate,prop);
                    else
                        oblate = [sysModel.P1.radius.value, sysModel.char.J2];
                        [traj_out] = prop_2bp(traj_in,options,oblate,prop);
                    end
                otherwise
                    error('That propagator has not been implemented yet.')
            end
        end

        %% Plotting
        function [ax, lineHandle] = plot(traj_in,varargin)
            % [ax, lineHandle] = plot(ST_in,varargin)
            % varargin:
            %   'col',[0.8 0.8 0.8] - changes background color
            %   plotbodies - true if no input ax
            %   'lagrange', true
            %   'ax', input axes
            %   'axLim', false - turns off auto scaling for Sun-B1 plots
            %   'addTitle', false
            %   'LTdir', true - adds low thrust arrows if they exist in the ST
            %
            % Alex Hoffman
            lineHandle = {};
            sysModel = traj_in.system_model;
            lagr = false;
            color = [0.8 0.8 0.8];
            pos_unit = strtok(traj_in.pos.unit,'_');
            switch sysModel.frame
                case 'B1centP1P2rot'
                    xlab = ['$x$ [' pos_unit ']'];
                    ylab = ['$y$ [' pos_unit ']'];
                    zlab = ['$z$ [' pos_unit ']'];
                case 'B2centP4B1rot'
                    xlab = ['$\tilde{x}$ [' pos_unit ']'];
                    ylab = ['$\tilde{y}$ [' pos_unit ']'];
                    zlab = ['$\tilde{z}$ [' pos_unit ']'];
                otherwise
                    warning('labels for this frame not yet implemented')
            end
            ax = [];
            axLim = false;
            addTitle = true;
            plotbodies = false;
            override = false;
            useDim = false; %plot with dimensional units
            for i = 1:length(varargin)/2
               switch varargin{i*2-1}
                   case 'col'
                       color = varargin{i*2}; 
                   case 'plotbodies'
                       plotbodies = varargin{i*2};
                   case 'lagrange'
                       lagr = varargin{i*2};
                   case 'ax'
                       ax = varargin{i*2};
                   case 'axLim'
                       axLim = varargin{i*2};
                   case 'addTitle'
                       addTitle = varargin{i*2};
                   case 'useDim'
                       useDim = varargin{i*2};
                   otherwise
                       error(['Variable argument ''' varargin{i*2-1} ''' not found.'])
               end
            end
            if isempty(ax)
                [~, ax] = ffigure;
                if ~override
                    plotbodies = true;
                end
            end
            axis(ax,'equal')
            xlabel(ax,xlab,'interpreter','latex')
            ylabel(ax,ylab,'interpreter','latex')
            zlabel(ax,zlab,'interpreter','latex');
            set(ax,'Color',color)
            hold(ax,'on')
            switch traj_in.system_model.frame %plot everything but the actual states
                case 'B1centP1P2rot'
                    CRplot(ax,sysModel,lagr,traj_in.pos.unit,plotbodies)
                    tit = [camelCase(sysModel.P1.name) '-' camelCase(sysModel.P2.name) ' rotating frame'];
                case 'B2centP4B1rot'
                    [oneMinusMu, P2radius] = SunB1plot(ax,traj_in,lagr,plotbodies);
                    tit = [camelCase(sysModel.P4.name) '-$B_1$ rotating frame'];
                    if axLim
                        zl = zlim(ax);
                        xlim(ax,oneMinusMu + 1.5*P2radius*[-1,1])
                        ylim(ax,1.5*P2radius*[-1,1])
                        zlim(ax,zl);
                    end
                case 'P1centinert'
                    sysModel.P1.plot(ax,sysModel,traj_in.pos.unit);
                    tit = 'ECI';
                case 'PcentPfixed'
                    sysModel.P1.plot(ax,sysModel,traj_in.pos.unit);
                    tit = 'ECEF';
                case 'P2centinert'
                    sysModel.P2.plot(ax,sysModel,traj_in.pos.unit);
                    tit = 'MCI';
                otherwise
                    error('This frame hasn''t been implemented yet. You''re on your own.')
            end
            %grab plotting chars
            dps = traj_in.plotting;
            markerStyle = dps.markerStyle;
            markerSize = dps.markerSize;
            dispName = strtok(traj_in.name,'/');
            posStates = traj_in.pos; 
            if useDim
                posStates = posStates.change_unit('km',sysModel);
            end
            posStates = posStates.value;
            % change marker details if necessary
            if size(posStates,2) == 1
                if strcmp(markerStyle,'.')
                    markerSize = 15;
                end
            else
                markerStyle = 'none';
                markerSize = 1;
            end
            %plot the states
            if size(posStates,2) == 1
                lh = plotStates(ax,posStates,...
                    'Color', dps.lineColor,...
                    'LineStyle', 'none',...
                    'LineWidth', dps.lineWidth,...
                    'Marker',markerStyle,...
                    'MarkerSize',markerSize,...
                    'DisplayName',dispName,...
                    'HandleVisibility',dps.HandleVisibility);
            else
                lh = plotStates(ax,posStates,...
                    'Color', dps.lineColor,...
                    'LineStyle', dps.lineStyle,...
                    'LineWidth', dps.lineWidth,...
                    'Marker',markerStyle,...
                    'MarkerSize',markerSize,...
                    'DisplayName',dispName,...
                    'HandleVisibility',dps.HandleVisibility);
            end
            lineHandle{1} = lh;
            tdir = traj_in.low_thrust.thrust_dir_history;
            throttle = traj_in.low_thrust.spacecraft.throttle;
            LTdir = dps.LTdir;
            if LTdir && throttle > 0
                LTplotmult = dps.LTplotmult;
                LTplotstep = dps.LTplotstep;
                statePlusThrust = posStates(1:3,:) + tdir.*LTplotmult.*throttle;
                lh = plot3(ax,[posStates(1,1:LTplotstep:end); statePlusThrust(1,1:LTplotstep:end)],...
                    [posStates(2,1:LTplotstep:end); statePlusThrust(2,1:LTplotstep:end)],...
                    [posStates(3,1:LTplotstep:end); statePlusThrust(3,1:LTplotstep:end)],...
                    'Color',dps.LTcolor,...
                    'LineStyle',dps.lineStyle,...
                    'LineWidth',dps.lineWidth/2,...
                    'HandleVisibility','off');
                lineHandle{end+1} = lh;
            end
            if addTitle
                title(ax,tit,'Interpreter','latex')
            end
        end

        %%
        function elem = conv2elem(obj)
            if ~strcmp(obj.system_model.system_dynamics,'2BP')
                error('Not a 2-body problem')
            end
            mu = ST.mu;
            rECI = ST.state(1:3,1);
            r = norm(rECI);
            vECI = ST.state(4:6,1);
            v = norm(vECI);
            x = [1; 0; 0];
            z = [0; 0; 1];
            hECI = cross(rECI,vECI);
            h = norm(hECI);
            nodeECI = cross(z,hECI)/h;
            n = norm(nodeECI);
            eECI = (1/mu)*((v^2-mu/r)*rECI-dot(rECI,vECI)*vECI);
            ecc = norm(eECI);
            p = (h^2)/mu;
            a = p/(1-ecc^2);
            T = 2*pi/sqrt(mu)*a^1.5;
            inc = acos(dot(z,hECI)/h);
            RAAN = acos(dot(x,nodeECI)/n);
            if nodeECI(2) < 0
                RAAN = 2*pi - RAAN;
            end
            if isnan(RAAN) %planar motion
                RAAN = 0; %arbitrary
                nodeECI = x;
            end
            argP = acos(dot(nodeECI,eECI)/(n*ecc));
            if dot(eECI,z) < 0
                argP = 2*pi - argP;
            end
            nu = acos(dot(eECI,rECI)/(ecc*r));
            if dot(rECI,vECI) < 0
                nu = 2*pi - nu;
            end
            u = argP + nu;
            l = RAAN + argP + nu;
            if isnan(argP)
                argP = 0;
            end
            if isnan(nu)
               nu = 0; 
            end
            out.a = a;
            out.ecc = ecc;
            % elem.energy = energy;
            out.T = T;
            out.p = p;
            out.inc = inc;
            out.RAAN = RAAN;
            out.argP = argP;
            out.nu = nu;
            out.mu = mu;
            out.u = u;
            out.l = l;
            elem = out;
        end

        %% Change frames
        function outTraj = changeFrame(obj,newFrame)
            currentFrame = obj.system_model.frame;
            if strcmp(currentFrame,newFrame)
                warning('Current frame and new frame are the same')
            end
            %anytime you add a frame, add it to the end of the list and add
            %a new row and column to signify connections (1 for self)
            frameNames = {'P1centinert', 'P2centinert', 'P4centinert',...
                'PcentPfixed', 'B1centP1P2rot', 'B1centP1P2inert',...
                'B2centP4B1rot', 'B2centP4B1inert', 'P1centP1P2rot',...
                'P2centP1P2rot', 'PQW'}; %ordered list of frames for graphFuncs
            graphFuncs = {1, [], [], @(x) eci_ecef(x), [], [], [], [], @(x) rot_inert(x), [], []
                          [], 1, [], @(x) eci_ecef(x), [], [], [], [], [], @(x) rot_inert(x), []
                          [], [], 1, [], [], [], [], [], [], [], []
                          @(x) eci_ecef(x), @(x) eci_ecef(x), [], 1, [], [], [], [], [], [], []
                          [], [], [], [], 1, @(x) rot_inert(x), @(x) em_p4b1(x), [], @(x) change_center(x,'P1'), @(x) change_center(x, 'P2'), []
                          [], [], [], [], @(x) rot_inert(x), 1, [], [], [], [], []
                          [], [], [], [], @(x) em_p4b1(x), [], 1, @(x) rot_inert(x), [], [], []
                          [], [], [], [], [], [], @(x) rot_inert(x), 1, [], [], []
                          @(x) rot_inert(x), [], [], [], @(x) change_center(x,'B1'), [], [], [], 1, [], []
                          [], @(x) rot_inert(x), [], [],  @(x) change_center(x,'B1'), [], [], [], [], 1, []
                          [], [], [], [], [], [], [], [], [], [], 1};
            graphEdges = logical(~cellfun('isempty',graphFuncs)); %define if route exists
            G = graph(graphEdges,frameNames); %create graph
            fpath = shortestpath(G,currentFrame,newFrame); %find path between frames
            if isempty(fpath)
                error(['No path found between ' currentFrame ' and ' newFrame])
            elseif length(fpath) == 1 %same frame
                outTraj = obj;
            else
                for i = 1:length(fpath)-1
                    cf = fpath{i};
                    if ~strcmp(cf, obj.system_model.frame)
                        error(['Expected current frame, ' cf ', does not match true current frame, ' obj.system_model.frame])
                    end
                    nf = fpath{i+1};
                    ci = strcmp(frameNames,cf);
                    ni = strcmp(frameNames,nf);
                    func = graphFuncs{ci,ni};
                    obj = func(obj);
                    if ~strcmp(nf, obj.system_model.frame)
                        error(['Expected new frame, ' nf ', does not match true new frame, ' obj.system_model.frame])
                    end
                end
                outTraj = obj;
            end
        end
        
        %% Functions to get various details
        function outTraj = extractSubTraj(obj,ind1,ind2)
            outTraj = obj;
            outTraj.time.value = outTraj.time.value(ind1:ind2);
            outTraj.pos.value = outTraj.pos.value(:,ind1:ind2);
            outTraj.vel.value = outTraj.vel.value(:,ind1:ind2);
            outTraj.low_thrust.mass.value = outTraj.low_thrust.mass.value(ind1:ind2);
        end

        function initTime = getInitTime(obj)
            initTime = obj.time;
            initTime.value = initTime.value(1);
        end

        function finTime = getFinalTime(obj)
            finTime = obj.time;
            finTime.value = finTime.value(end);
        end

        function tspan = getTimeSpan(obj)
            tspan = obj.time;
            if length(tspan.value) < 2
                error('Length of time must be greater than 1')
            elseif tspan.value(1) == tspan.value(end)
                error('time span must have nonzero span')
            end
            tspan.value = [tspan.value(1), tspan.value(end)];
        end

        function duration = getDuration(obj)
            duration = obj.time;
            if length(duration.value) < 2
                error('Length of time must be greater than 1')
            elseif duration.value(1) == duration.value(end)
                error('time span must have nonzero span')
            end
            duration.value = duration.value(end) - duration.value(1);
        end

        function time = getTimeByIndex(obj,index)
            time = obj.time;
            time.value = time.value(index);
        end

        function initPos = getInitPos(obj)
            initPos = obj.pos;
            initPos.value = initPos.value(:,1);
        end

        function finPos = getFinalPos(obj)
            finPos = obj.pos;
            finPos.value = finPos.value(:,end);
        end

        function initVel = getInitVel(obj)
            initVel = obj.vel;
            initVel.value = initVel.value(:,1);
        end

        function finVel = getFinalVel(obj)
            finVel = obj.vel;
            finVel.value = finVel.value(:,end);
        end

        function state = getStateByIndex(obj,index)
            if index == -1
                index = size(obj.pos.value,2);
            end
            state = [obj.pos.value(:,index); obj.vel.value(:,index);...
                obj.low_thrust.mass.value(index)];
        end

        function state = getPosByIndex(obj,index)
            if index == -1
                index = size(obj.pos.value,2);
            end
            state = c_dim_quant(obj.pos.value(:,index),obj.pos.unit);
        end

        function state = getVelByIndex(obj,index)
            if index == -1
                index = size(obj.vel.value,2);
            end
            state = c_dim_quant(obj.vel.value(:,index),obj.vel.unit);
        end

        function state = getMassByIndex(obj,index)
            if index == -1
                index = size(obj.low_thrust.mass.value,2);
            end
            state = c_dim_quant(obj.low_thrust.mass.value(:,index),obj.low_thrust.mass.unit);
        end

        function initMass = getInitMass(obj)
            initMass = obj.low_thrust.mass;
            initMass.value = initMass.value(1);
        end
        
        function finMass = getFinalMass(obj)
            finMass = obj.low_thrust.mass;
            finMass.value = finMass.value(end);
        end

        function state = getAllStates(obj)
            state = [obj.pos.value(:,:); obj.vel.value(:,:);...
                obj.low_thrust.mass.value(:,:)];
        end

        function time = getAllTime(obj)
            time = obj.time.value(:);
        end

        function len = getLength(obj)
            len = length(obj.time.value);
        end

        function [states, terr] = getStatesAfterDuration(obj,duration)
            timeVec = obj.time;
            duration = duration.change_unit(timeVec.unit,obj.system_model);
            totalDuration = getDuration(obj);
            if abs(duration.value) > abs(totalDuration.value)
                error('Not within trajectory duration')
            end
            [terr,minInd] = min(timeVec.value - timeVec.value(1) - duration.value);
            states = getStateByIndex(obj,minInd);
        end

    end
end
%%
function CRplot(AxObj,sysModel,lagr,state_unit,plotbodies)
if strcmpi(sysModel.P1.name,'sun') && strcmpi(sysModel.P2.name,'earth')
   warning('Did you mean to use Sun-B1 frame?') 
end
if plotbodies
    sysModel.P1.plot(AxObj,sysModel,state_unit);
    sysModel.P2.plot(AxObj,sysModel,state_unit);
end
if lagr
    mu3 = sysModel.char.mu;
    l_s = sysModel.char.lstar;
    if strcmp(state_unit,'km')
        dim_mult = l_s;
    elseif strcmp(state_unit,'nd')
        dim_mult = 1;
    else
        error('this unit not supported yet')
    end
   [L1,L2,L3,L4,L5] = lagrangePoints(mu3);
   lagPoints = [L1, L2, L3, L4, L5];
   plot3(AxObj,lagPoints(1,:).*dim_mult,...
       lagPoints(2,:).*dim_mult,lagPoints(3,:).*dim_mult,...
       'k*','HandleVisibility','off','markers',12)
end
end
%%
function [oneMinusMu, P2radius] = SunB1plot(AxObj,traj_in,lagr,plotbodies)
sysModel = traj_in.system_model;
tvec = traj_in.time.value;
muP1P2 = sysModel.char.muP1P2;
mu = sysModel.char.mu;
lstarP1P2 = sysModel.char.lstarP1P2;
lstar = sysModel.char.lstar;
descNodeP4 = sysModel.char.descNodeP4;
incP4 = sysModel.char.incP4;
linewidth = 1;

% % old method, just made circles even when out of plane
P1radius = muP1P2*lstarP1P2/lstar;
P2radius = (1-muP1P2)*lstarP1P2/lstar;
oneMinusMu = 1 - mu;
rads = linspace(0,2*pi,100);
x = cos(rads); y = sin(rads); z = zeros(1,length(rads));
xP1 = P1radius*x; yP1 = P1radius*y;
xP2 = P2radius*x; yP2 = P2radius*y;
posP1 = [xP1;yP1;z];
posP2 = [xP2;yP2;z];
%rotate P1 and P2 into OoP orbits if necessary
rot1 = DCM(-incP4,'y',false);
rot2 = DCM(2*pi-descNodeP4,'z',false);
for i = 1:length(z)
   posP1(:,i) = rot2*rot1*posP1(:,i);
   posP2(:,i) = rot2*rot1*posP2(:,i);
end
posP1(1,:) = posP1(1,:) + oneMinusMu;
posP2(1,:) = posP2(1,:) + oneMinusMu;

objs = findall(AxObj);
alreadyPlotted = false;
for i = 1:length(objs) %check to see if the Sun is already there
   if isa(objs(i),'matlab.graphics.chart.primitive.Line')
       if strcmp(objs(i).DisplayName,camelCase(sysModel.P1.name))
          alreadyPlotted = true; 
          break
       end
   end
end
if ~alreadyPlotted
    plotStates(AxObj,posP1,'-',...
        'color',colour('blue'),...
        'LineWidth',linewidth/2,...
        'DisplayName',camelCase(sysModel.P1.name),...
        'HandleVisibility','off')
    plotStates(AxObj,posP2,'-',...
        'color',colour('gray'),...
        'LineWidth',linewidth/2,...
        'DisplayName',camelCase(sysModel.P2.name))
end

% % new method, converts points from EM to SB1 frame
% P1EM = repmat([-muP1P2;0;0;0;0;0],1,length(tvec));
% P2EM = repmat([1-muP1P2;0;0;0;0;0],1,length(tvec));
% ST_P1_EM = ST(tvec,P1EM,'nd','nd',sysModel,'B1centP1P2rot');
% ST_P2_EM = ST(tvec,P2EM,'nd','nd',sysModel,'B1centP1P2rot');
% ST_P1_P4B1 = CONVframe(ST_P1_EM,'B2centP4B1rot');
% ST_P2_P4B1 = CONVframe(ST_P2_EM,'B2centP4B1rot');
% P1state = ST_P1_P4B1.STATE;
% P2state = ST_P2_P4B1.STATE;
% plotStates(AxObj,P1state,'-','color',colour('blue'),'LineWidth',linewidth,'DisplayName',camelCase(sysModel.P1))
% plotStates(AxObj,P2state,'-','color',colour('gray'),'LineWidth',linewidth,'DisplayName',camelCase(sysModel.P2))

if plotbodies
    alreadyPlotted = false;
    for i = 1:length(objs) %check to see if the Sun is already there
       if isa(objs(i),'matlab.graphics.chart.primitive.Surface')
           if strcmp(objs(i).DisplayName,camelCase(sysModel.P4.name))
              alreadyPlotted = true; 
              break
           end
       end
    end
    if ~alreadyPlotted
        P4_data = sysModel.P4;
        P4_r_nd = P4_data.radius.value/lstar;
        H = plot3DBody(AxObj,lower(sysModel.P4.name), P4_r_nd, [-mu 0 0]);
        set(H,'HandleVisibility','off','DisplayName',camelCase(sysModel.P4.name)) 
    end
    if lagr
       [L1,L2,L3,L4,L5] = lagrangePoints(mu);
       lagPoints = [L1, L2, L3, L4, L5];
       plot3(AxObj,lagPoints(1,:),lagPoints(2,:),lagPoints(3,:),...
           'k*','HandleVisibility','off','markers',12)
    end
end
end
