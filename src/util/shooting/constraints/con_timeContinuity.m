function con = con_timeContinuity(segment)
%con = con_timeContinuity(segment)
% input the segment which requires time continuity with its terminal node
whichSeg = segment.seg_index;
whichNodes = segment.node_indices;
F = @(segs,nodes) conF(segs,nodes);
partials = @(segs,nodes,desPartials) calc_partials(segs,nodes,desPartials);
conname = ['Time Continuity Constraint: seg ' num2str(whichSeg)];
con = c_constraint(conname,F,partials,whichSeg,whichNodes);
end
%% constraint evaluation
function F = conF(segs,nodes)
    [T1,T2,dt] = getTimes(segs,nodes);
    currDT = T2.sub(T1);
    F = currDT.sub(dt); %T2 - T1 - dt 
    F = F.value;
end

%% partials
function partials = calc_partials(segs,nodes,desPartials)
    partials = nan(1,length(desPartials));
    for i = 1:length(desPartials)
        desiredPartialCell = desPartials(i,:);
        desiredPartial = desiredPartialCell{3};
        thisPartial = nan;
        switch desiredPartialCell{1}
            case 'seg'
                switch desiredPartial %string of which partial to take
                    case 'dt'
                        thisPartial = -1;
                    case {'thrust','SOCalpha0',...
                            'SOCalphadot','SOCbeta0','SOCbetadot'}
                        thisPartial = 0;
                    otherwise
                        error(['Partial "' desiredPartial...
                            '" does not exist for Time Continuity constraint'])
                end
            case 'node'
                thisNode = nodes(desiredPartialCell{2});
                if segs{1}.node_indices(2) == thisNode{1}.node_index
                    switch desiredPartial
                        case {'x','y','z','xd','yd','zd','mass'}
                            thisPartial = 0;
                        case 't0'
                            thisPartial = 0; %1??
                        otherwise
                            error(['Partial "' desiredPartial ...
                                '" does not exist for Time continuity constraint'])
                    end
                else
                    switch desiredPartial
                        case 't0'
                            thisPartial = 0;
                        case {'x','y','z','xd','yd','zd','mass'}
                            thisPartial = 0;
                        otherwise
                            error(['Partial "' desiredPartial ...
                                '" does not exist for Time continuity constraint'])
                    end
                end
        end
        partials(i) = thisPartial;
    end    
end
%% utility
function [T1,T2,dt] = getTimes(segs,nodes)
    %segs should be empty, should only have 1 node
    if length(segs) ~=1 || length(nodes) ~= 2
        error('incorrect segments or nodes input to Time Continuity')
    end
    finalNode = nodes{2};
    initNode = nodes{1};
    if finalNode.node_index ~= segs{1}.node_indices(2)
        error('Time Continuity constraint not using final node')
    end
    if ~strcmp(initNode.time.unit,finalNode.time.unit)
        error('Mismatched units')
    end
    T1 = initNode.time;
    T2 = finalNode.time;
    dt = segs{1}.getDuration();
    T1 = T1.change_unit(dt.unit,initNode.system_model);
    T2 = T2.change_unit(dt.unit,finalNode.system_model);
end