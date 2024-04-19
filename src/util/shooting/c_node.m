classdef c_node < c_traj
    % obj = c_(name, time, pos, vel, sysModel, varargin)
    % varargs:
    % {traj_args, node_args}
    % node_args: free_vars, node_index
    % free_vars: cell of strings: could contain any of {'init_cond','t0','dt','mass','thrust'}
            
    properties
        free_vars
        node_index
    end
    
    methods
        function obj = c_node(name, time, pos, vel, sysModel, varargin)
            % obj = c_(name, time, pos, vel, sysModel, varargin)
            % varargs:
            % {traj_args, node_args}
            % node_args: free_vars, node_index
            % free_vars: cell of strings: could contain any of 
            %   {'x','y','z','xd','yd','zd','t0','mass'}
            if isempty(varargin)
                varargin = {{},{}};
            elseif length(varargin) ~=2
                error('varargin must be ...,{traj_args},{seg_args}')
            end
            obj = obj@c_traj(name, time, pos, vel, sysModel, varargin{1}{:});
            varargin = varargin{2};
            for i = 1:length(varargin)/2
               switch varargin{i*2-1}
                   case 'free_vars'
                       obj.free_vars = varargin{i*2};
                   case 'node_index'
                       obj.node_index = varargin{i*2};
                   otherwise
                       error('Unknown varargin')
               end
            end
        end

        function [node_out,varargout] = prop(traj_in,varargin)
            %don't need to propagate nodes
            node_out = traj_in;
            varargout = {};
        end
        
        function [val,isFree] = checkIfFree(obj,valString,x,fixedValue) 
            fvs = obj.free_vars;
            % used in genObjFunc and when recreating the results from genObjFunc
            if ~isempty(fvs)
                idx = find(ismember(fvs(:,1),valString));
                if isempty(idx)
                   val = fixedValue;
                   isFree = false;
                else
                   xidx = fvs{idx,2};
                   val = x(xidx);
                   isFree = true;
                end
            else
                val = fixedValue;
                isFree = false;
            end
        end
        
        function [traj_out] = convToTraj(obj)
            traj_out = c_traj(obj.name,obj.time,obj.pos,obj.vel,obj.system_model,...
                'STM',obj.stm,'plotting',obj.plotting);
            traj_out.low_thrust = obj.low_thrust;
            traj_out.etc = obj.etc;
        end
    end
end
