function con = con_stateContinuity(segment,stateIndex,varargin)
isOffset = false;
desOffset = 0;
if ~isempty(varargin)
    isOffset = true;
    desOffset = varargin{1};
    checkClass(desOffset,'c_dim_quant')
    if ~any(strcmp(desOffset.unit,{'nd_l','nd_v'}))
        error('offset must be defined in nondimensional')
    end
    desOffset = desOffset.value;
end
whichSeg = segment.seg_index;
whichNodes = segment.node_indices;
F = @(segs,nodes) conF(segs,nodes,stateIndex,desOffset);
partials = @(segs,nodes,desPartials) calc_partials(segs,nodes,desPartials,stateIndex);
stateInds = {'x', 'y', 'z', 'xd', 'yd', 'zd', 'mass'};
conname = ['Continuity Constraint: seg ' num2str(whichSeg) ', ' stateInds{stateIndex}];
if isOffset
    conname = [conname ', offset: ' num2str(desOffset)];
end
con = c_constraint(conname,F,partials,whichSeg,whichNodes);
end
%% constraint evaluation
function F = conF(segs,nodes,stateIndex,desOffset)
    [sSeg,sNode] = getStates(segs,nodes);
    if stateIndex == 6
        a=1;
    end
    F = (sSeg(stateIndex) - sNode(stateIndex)) - desOffset; 
end

%% partials
function partials = calc_partials(segs,nodes,desPartials,stateIndex)
    partials = nan(1,size(desPartials,1));
    [sSeg,~] = getStates(segs,nodes);
    STM = segs{1}.stm;
    stateString = STM.key{stateIndex};
    for i = 1:size(desPartials,1)
        desiredPartialCell = desPartials(i,:);
        desiredPartial = desiredPartialCell{3};
        thisPartial = nan;
        switch desiredPartialCell{1}
            case 'seg'
                switch desiredPartial %string of which partial to take
                    case 'dt'
                        if stateIndex <= 3
                            thisPartial = sSeg(stateIndex + 3);
                        elseif stateIndex <= 7
                            thisPartial = sSeg(stateIndex + 4); %account for mass
                        else
                            error('State index must be less than or equal to 7')
                        end
                    case {'thrust','SOCalpha0',...
                            'SOCalphadot','SOCbeta0','SOCbetadot'}
                        thisPartial = nan; %not correct
%                         warning('SOC partials not yet coded')
                    otherwise
                        error(['Partial "' desiredPartial...
                            '" does not exist for Continuity constraint'])
                end
            case 'node'
                thisNode = nodes(desiredPartialCell{2});
                if segs{1}.node_indices(2) == thisNode{1}.node_index %terminal node
                    switch desiredPartial
                        case STM.key
                            if strcmp(stateString,desiredPartial)
                                thisPartial = -1;
                            elseif strcmp(desiredPartial,'t0')
                                if stateIndex <= 3
                                    thisPartial = sSeg(stateIndex + 3);
                                elseif stateIndex <= 7
                                    thisPartial = sSeg(stateIndex + 4); %account for mass
                                else
                                    error('State index must be less than or equal to 7')
                                end
                            else
                                thisPartial = 0;
                                if strcmp(desiredPartial,'mass')
                                    thisPartial = 0;
                                end
                            end
                        otherwise
                            error(['Partial "' desiredPartial ...
                                '" does not exist for continuity constraint'])
                    end
                else %initial node
                    switch desiredPartial
                        case STM.key
                            thisPartial = STM.get(stateString,desiredPartial);
                            if strcmp(desiredPartial,'t0')
                                if stateIndex <= 3
                                    thisPartial = -sSeg(stateIndex + 3);
                                elseif stateIndex <= 7
                                    thisPartial = -sSeg(stateIndex + 4); %account for mass
                                else
                                    error('State index must be less than or equal to 7')
                                end
                                thisPartial = nan;
                            end
                        otherwise
                            error(['Partial "' desiredPartial ...
                                '" does not exist for continuity constraint'])
                    end
                end
        end
        partials(i) = thisPartial;
    end    
end
%% utility
function [sSeg,sNode] = getStates(segs,nodes)
    %segs should be empty, should only have 1 node
    if length(segs) ~=1 || length(nodes) ~= 2
        error('incorrect segments or nodes input to Continuity')
    end    
    if segs{1}.seg_index == 3
        a = 1;
    end
    seg = segs{1};
    sSeg = seg.getStateByIndex(-1);
    acc = seg.EOM(size(seg.pos.value,2)); %for last state
    Tmax = seg.low_thrust.spacecraft.TmaxND.value * ...
        seg.low_thrust.spacecraft.throttle;
    Isp = seg.low_thrust.spacecraft.IspND.value;
    g0 = seg.low_thrust.spacecraft.g0ND.value;
    mdot = -Tmax/Isp/g0;

    sSeg = [sSeg(1:7); acc.value; mdot];
    finalNode = nodes{2};
    if finalNode.node_index ~= segs{1}.node_indices(2)
        error('Continuity constraint not using final node')
    end
    sNode = finalNode.getStateByIndex(1);
end