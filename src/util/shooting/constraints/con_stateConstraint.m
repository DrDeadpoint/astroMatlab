function con = con_stateConstraint(node,stateIndex,desState)
whichNode = node.node_index;
F = @(segs,nodes) conF(segs,nodes,stateIndex,desState);
partials = @(segs,nodes,desPartials) calc_partials(segs,nodes,desPartials,stateIndex);
conname = ['State Constraint: node ' num2str(whichNode) ...
    ', state ' num2str(stateIndex) ', desState ' num2str(desState)];
con = c_constraint(conname,F,partials,[],whichNode);
end
%% constraint evaluation
function F = conF(segs,nodes,stateIndex,desState)
    s = getState(segs,nodes,stateIndex);
    F = s-desState; 
end

%% partials
function partials = calc_partials(segs,nodes,desPartials,stateIndex)
    partials = nan(1,length(desPartials));
    for i = 1:length(desPartials)
        desiredPartialCell = desPartials(i,:);
        if ~strcmp(desiredPartialCell(1),'node')
            error('Something went wrong')
        end
        desiredPartial = desiredPartialCell{3};
        switch desiredPartial %string of which partial to take
            case {'x','y','z','xd','yd','zd'}
                if stateIndex == find(ismember(...
                        {'x','y','z','xd','yd','zd'},desiredPartial))
                    thisPartial = 1;
                else
                    thisPartial = 0;
                end
            case {'t0','mass','dt','SOCalpha0',...
                    'SOCalphadot','SOCbeta0','SOCbetadot'}
                thisPartial = 0;
            otherwise
                error(['Partial "' desiredPartial '" does not exist for State constraint'])
        end
        partials(i) = thisPartial;
    end    
end
%% utility
function s = getState(segs,nodes,stateIndex)
    %segs should be empty, should only have 1 node
    if ~isempty(segs) || length(nodes) ~= 1
        error('incorrect segments or nodes input to State Constraint')
    end
    s = nodes{1}.getStateByIndex(1);
    s = s(stateIndex);
end