function con = con_dV(segment,dVmag)
whichSeg = segment.seg_index;
whichNodes = segment.node_indices;
F = @(segs,nodes) conF(segs,nodes,dVmag);
partials = @(segs,nodes,desPartials) calc_partials(segs,nodes,desPartials);
conname = ['DV Constraint: seg ' num2str(whichSeg) ', dV = ' num2str(dVmag)];
con = c_constraint(conname,F,partials,whichSeg,whichNodes);
end
%% constraint evaluation
function F = conF(segs,nodes,dVmag)
    [sSeg,sNode] = getStates(segs,nodes);
    F = norm(sSeg - sNode) - dVmag; 
end

%% partials
function partials = calc_partials(segs,nodes,desPartials)
    partials = nan(1,size(desPartials,1));
end
%% utility
function [sSeg,sNode] = getStates(segs,nodes)
    %segs should be empty, should only have 1 node
    if length(segs) ~=1 || length(nodes) ~= 2
        error('incorrect segments or nodes input to Continuity')
    end
    sSeg = segs{1}.getFinalVel();
    finalNode = nodes{2};
    if finalNode.node_index ~= segs{1}.node_indices(2)
        error('Continuity constraint not using final node')
    end
    sNode = finalNode.getInitVel();
end