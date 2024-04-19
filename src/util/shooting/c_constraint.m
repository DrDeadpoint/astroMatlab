classdef c_constraint
   
    properties
        name
        func
        partials
        which_segments
        which_nodes
        con_equality
    end
    
    methods
        function con = c_constraint(name,F,partialList,whichSegs,whichNodes,varargin)
            %con = c_constraint(name,F,partialList,whichSegs,whichNodes)
            %F is an anonymous function, generally of the form @(allTrajs) F(allTrajs,...)
            %partialList is a column cell of strings denoting which partials to take
            con.name = name;
            con.func = F;
            con.partials = partialList; %will produce a 2 column cell array of partials
            con.which_segments = whichSegs;
            con.which_nodes = whichNodes;
            if isempty(varargin)
                con.con_equality = [];
            else
                con.con_equality = varargin{1};
            end
        end
        
        function F = eval_f(obj,segs,nodes)
            %Evaluate the constraint
            %segs and nodes will always be column cell arrays (can be empty)
            desSegs = segs(obj.which_segments);
            desNodes = nodes(obj.which_nodes);
            F = obj.func(desSegs,desNodes);
        end

        function DFrow = eval_df(obj,segs,nodes,x0key)
            %determine which partials are needed
            desSegs = segs(obj.which_segments);
            desNodes = nodes(obj.which_nodes);
            desPartials = {};
            desPartialsIndices = [];
            for i = 1:length(x0key)
                thisFV = x0key{i};
                segOrNode = thisFV(1:3);
                thisFV = thisFV(4:end);
                [ind,rem] = strtok(thisFV,' ');
                FVstr = rem(2:end);
                switch segOrNode
                    case 'seg'
                        desSans = desSegs;
                        whichSan = 'which_segments';
                    case 'nod'
                        desSans = desNodes;
                        whichSan = 'which_nodes';
                        segOrNode = 'node'; %looks nicer
                end
                for j = 1:length(desSans)
                    if obj.(whichSan)(j) == str2double(ind)
                        desPartials = [desPartials;...
                            {segOrNode, j, FVstr}];
                        desPartialsIndices = [desPartialsIndices,i];
                    end
                end
            end

            %Evalutate the partials given a list of desired ones
            partials_eval = obj.partials(desSegs,desNodes,desPartials);
            DFrow = zeros(1,length(x0key));
            %search through allPartials for ones listed in x0key
            for i = 1:length(partials_eval)
                DFrow(desPartialsIndices(i)) = partials_eval(i); %x0key is a cell of strings corresponding to fieldnames
            end
        end
    end
end

