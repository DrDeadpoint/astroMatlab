function con = con_primaryDistance(node,whichPrimary,desDistance)
whichNode = node.node_index;
checkClass(desDistance,'c_dim_quant')
if length(desDistance.value) ~= 1
    error('only single value can be input')
end
F = @(segs,nodes) conF(segs,nodes,whichPrimary,desDistance);
partials = @(segs,nodes,desPartials) calc_partials(segs,nodes,desPartials,whichPrimary);
conname = ['Primary Distance Constraint: node ' num2str(whichNode) ...
    ', primary ' num2str(whichPrimary) ', distance ' num2str(desDistance.value)];
con = c_constraint(conname,F,partials,[],whichNode);
end
%% constraint evaluation
function F = conF(segs,nodes,whichPrimary,desDistance)
    [r,r_prim,~] = getState(segs,nodes,whichPrimary);
    desDistance = desDistance.change_unit(r.unit,nodes{1}.system_model);
    F = norm(r.value-r_prim.value)^2 - desDistance.value^2; 
end

%% partials
function partials = calc_partials(segs,nodes,desPartials,whichPrimary)
    partials = nan(1,length(desPartials));
    [r,r_prim,v_prim] = getState(segs,nodes,whichPrimary);
    r = r.value;
    r_prim = r_prim.value;
    for i = 1:length(desPartials)
        desiredPartialCell = desPartials(i,:);
        if ~strcmp(desiredPartialCell(1),'node')
            error('Something went wrong')
        end
        desiredPartial = desiredPartialCell{3};
        switch desiredPartial %string of which partial to take
            case 'x'
                thisPartial = 2*(r(1) - r_prim(1));
            case 'y'
                thisPartial = 2*(r(2) - r_prim(2));
            case 'z'
                thisPartial = 2*(r(3) - r_prim(3));
            case 't0'
                thisPartial = -2*dot(r-r_prim,v_prim);
            case {'xd','yd','zd','mass','dt','SOCalpha0',...
                    'SOCalphadot','SOCbeta0','SOCbetadot'}
                thisPartial = 0;
            otherwise
                error(['Partial "' desiredPartial '" does not exist for Primary Distance constraint'])
        end
        partials(i) = thisPartial;
    end    
end
%% utility
function [r,r_prim,v_prim] = getState(segs,nodes,whichPrimary)
    %segs should be empty, should only have 1 node
    if ~isempty(segs) || length(nodes) ~= 1
        error('incorrect segments or nodes input to Primary Distance')
    end
    node = nodes{1};
    r = node.getInitPos();
    switch node.system_model.frame
        case 'B1centP1P2rot'
            % dot product between relative position and relative velociy vectors
            % F = (nodePos - primaryPos) * (nodeVel - primaryVel)
            mu = node.system_model.char.mu;
            switch whichPrimary
                case 1
                    r_prim = [-mu;0;0];
                    v_prim = [0;0;0];
                case 2
                    r_prim = [1-mu;0;0];
                    v_prim = [0;0;0];
                otherwise
                    error('Not a valid primary')
            end
            r_prim = c_dim_quant(r_prim,'nd_l');
            r_prim = r_prim.change_unit(r.unit,node.system_model);
        case {'P1centinert','J2000'}
            if whichPrimary == 1
                r_prim = c_dim_quant([0;0;0],r.unit);
                v_prim = [0;0;0];
            else
                error('can''t implement inertial frame for other primary yet')
            end
        otherwise
            error('frame not implemented')
    end
end