classdef c_problem_parameter
    %PROBLEM_PARAMETER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        segments
        nodes
        constraints
        x0
        x0key
    end
    
    methods
        function obj = c_problem_parameter(segments,nodes,constraints,varargin)
            %PROBLEM_PARAMETER(segments,nodes,constraints)
            if iscell(segments)
                obj.segments = segments;
            else
                error('segments must be a cell array of c_segment''s')
            end
            if iscell(nodes)
                obj.nodes = nodes;
            else
                error('segments must be a cell array of c_node''s')
            end
            if iscell(constraints)
                obj.constraints = constraints;
            else
                error('Constraints must be a cell array of function_handles')
            end
        end
       
    end
end

