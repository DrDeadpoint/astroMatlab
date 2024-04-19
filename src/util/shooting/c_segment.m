classdef c_segment < c_traj
    % obj = c_segment(name, time, pos, vel, sysModel, varargin)
    % varargs:
    % {traj_args, seg_args}
    % seg_args: free_vars, events, seg_index, node_indices
    % free_vars: cell of strings: could contain any of {'init_cond','t0','dt','mass','thrust'}
            
    properties
        free_vars
        event_list
        seg_index
        node_indices
    end
    
    methods
        function obj = c_segment(name, time, pos, vel, sysModel, varargin)
            % obj = c_segment(name, time, pos, vel, sysModel, varargin)
            % varargs:
            % {traj_args, seg_args}
            % seg_args: free_vars, events, seg_index, node_indices
            % free_vars: cell of strings: could contain any of 
            %   {'dt','thrust',claw coefficients}
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
                   case 'seg_index'
                       obj.seg_index = varargin{i*2};
                   case 'node_indices'
                       obj.node_indices = varargin{i*2};
                   case 'events'
                       obj.event_list = varargin{i*2};
                   otherwise
                       error('Unknown varargin')
               end
            end
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
        
    end
end

