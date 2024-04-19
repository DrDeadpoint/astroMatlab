classdef c_stm
    %state_transition_matrix Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        matrix
        key
    end
    
    methods
        function obj = c_stm(matrix,key)
            %state_transition_matrix Construct an instance of this class
            %   Detailed explanation goes here
            obj.matrix = matrix;
            obj.key = key;
        end
        
        function STMitem = get(obj,partialOf,withRespectTo)
            %STMitem = get(obj,partialOf,withRespectTo)
            indRow = find(ismember(obj.key,partialOf));
            indCol = find(ismember(obj.key,withRespectTo));
            if isempty(indRow) || isempty(indCol)
                error('Desired STM item does not exist')
            end
            STMitem = obj.matrix(indRow,indCol);
        end
    end
end

