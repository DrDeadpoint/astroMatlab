classdef c_problem_setup
    % c_problem_setup(param,solverOptions,fcnOptions)
    %   'param' is a c_problem_parameter of segments and constraints
    %       segments are built using the c_segment class
    %       constraints contains a cell array of anonymous constraint functions {@(alllTrajs) con_XYZ(allTrajs,...)}
    %   'solverOptions' can be empty or contain user defined solver options (for Newton Raphson, for eample)
    %   'fcnOptions'  can be empty or contain options for solving the objective function
    %
    % Alex Hoffman
    % 11/08/2021
        
    properties
        objFunc
        param
        options
    end
    
    methods
        function obj = c_problem_setup(param,varargin)
            solverOptions = struct;
            fcnOptions = struct;
            for i = 1:length(varargin)/2
                switch varargin{i*2-1}
                    case 'solverOptions'
                        solverOptions = varargin{i*2};
                    case 'fcnOptions'
                        fcnOptions = varargin{i*2};
                end
            end
            solverOptionDefaults = struct('maxIter',30,...
                                          'fcnTol',1e-12,...
                                          'verbose',true,...
                                          'maxError',1e+7,...
                                          'notConvError',1e+2,...
                                          'attenuationFactor',1,...
                                          'normType',2,...
                                          'invMethod','qr',...
                                          'weights',[]);
            f = fieldnames(solverOptionDefaults);
            for i = 1:length(f)
                if isfield(solverOptions,f{i})
                    continue
                else
                    solverOptions.(f{i}) = solverOptionDefaults.(f{i});
                end
            end
            % set function default options where the input doesn't have one
            integratorOptions = odeset('reltol', 1e-12, 'abstol', 1e-12);
            fcnOptionDefaults = struct('integratorOptions',integratorOptions,...
                                       'DFmethod','analytic',...
                                       'step_centralDifferencing',1e-7,...
                                       'propMethod','manual' ... % 'atd' or 'manual' propagator method
                                       );
            f = fieldnames(fcnOptionDefaults);
            for i = 1:length(f)
                if isfield(fcnOptions,f{i})
                    continue
                else
                    fcnOptions.(f{i}) = fcnOptionDefaults.(f{i});
                end
            end
            options = struct('solverOptions',solverOptions,'fcnOptions',fcnOptions);
            %% build x0 based on free vars
            segments = param.segments; 
            nodes = param.nodes;
            segsAndNodes = {segments,nodes};
            sanstrings = {'seg';'nod'};
            x0 = [];
            x0key = {};
            for ii = 1:2
                sans = segsAndNodes{ii};
                sanstr = sanstrings{ii};
                for i = 1:length(sans)
                    san = sans{i};
                    fvs = san.free_vars;
                    if ~iscolumn(fvs) && ~isempty(fvs)
                        error('free variables must be a column cell vector')
                    end
                    for j = 1:length(fvs)
                       xLen = length(x0); %length of free var vector before adding more
                       fnames = fieldnames(san.low_thrust.control_law.coeffs);
                       switch fvs{j,1}
                           case 'x'
                               x0 = [x0; san.pos.value(1,1)];
                               x0key{end+1,1} = [sanstr num2str(i) ' x'];
                           case 'y'
                               x0 = [x0; san.pos.value(2,1)];
                               x0key{end+1,1} = [sanstr num2str(i) ' y'];
                           case 'z'
                               x0 = [x0; san.pos.value(3,1)];
                               x0key{end+1,1} = [sanstr num2str(i) ' z'];
                           case 'xd'
                               x0 = [x0; san.vel.value(1,1)];
                               x0key{end+1,1} = [sanstr num2str(i) ' xd'];
                           case 'yd'
                               x0 = [x0; san.vel.value(2,1)];
                               x0key{end+1,1} = [sanstr num2str(i) ' yd'];
                           case 'zd'
                               x0 = [x0; san.vel.value(3,1)];
                               x0key{end+1,1} = [sanstr num2str(i) ' zd'];
                           case 't0'
                               x0 = [x0; san.time.value(1)];
                               x0key{end+1,1} = [sanstr num2str(i) ' t0'];
                           case 'dt'
                               x0 = [x0; san.time.value(end) - san.time.value(1)];
                               x0key{end+1,1} = [sanstr num2str(i) ' dt'];
                           case 'mass'
                               x0 = [x0; san.low_thrust.mass.value(1)];
                               x0key{end+1,1} = [sanstr num2str(i) ' mass'];
                           case 'thrust'
                               x0 = [x0; san.low_thrust.spacecraft.Tmax.value];
                               x0key{end+1,1} = [sanstr num2str(i) ' thrust'];
                           case fnames
                               k = 1;
                               claw = san.low_thrust.control_law.law;
                               while ~strcmp(fvs{j,1},fnames{k})
                                    k = k+1;
                               end
                               x0 = [x0; san.low_thrust.control_law.coeffs.(fnames{k})];
                               x0key{end+1,1} = [sanstr num2str(i) ' ' claw fnames{k}];
                           otherwise
                               error(['Unknown free variable: "' fvs{j,1} '" on ' sanstr ' ' num2str(i)])
                       end
                       fvs{j,2} = (xLen+1) : (xLen+1);
                    end
                    san.free_vars = fvs;
                    switch ii
                        case 1
                            param.segments{i} = san;
                        case 2
                            param.nodes{i} = san;
                    end
                end
            end
            %check for inequality constraints, add slack variables
            for ii = 1:length(param.constraints)
                con = param.constraints{ii};
                if strcmp(con.name(1:10),'Inequality')
                    
                end
            end
            %% build structure 
            param.x0 = x0;
            param.x0key = x0key;
            objFunc = @(x) genObjFunc(x,param,options.fcnOptions);
            obj.objFunc = objFunc;
            obj.param = param;
            obj.options = options;
        end
        
    end
end

