function con = con_primaryApse(node,whichPrimary)
whichNode = node.node_index;
F = @(segs,nodes) conF(segs,nodes,whichPrimary);
partials = @(segs,nodes,desPartials) calc_partials(segs,nodes,desPartials,whichPrimary);
conname = ['Apse Constraint: '...
    'node ' num2str(whichNode) ', primary ' num2str(whichPrimary)];
con = c_constraint(conname,F,partials,[],whichNode);
end
%% constraint evaluation
function F = conF(segs,nodes,whichPrimary)
    [r,v,r_p,v_p,~] = getrv(segs,nodes,whichPrimary);
    r_rel = r-r_p;
    v_rel = v-v_p;
    % if i want to use this for for than the P1P2 rotating frame, need to
    % account for motion of the primaries as well (location and velocity)
    F = dot(r_rel,v_rel); 
end

%% partials
function partials = calc_partials(segs,nodes,desPartials,whichPrimary)
    partials = nan(1,length(desPartials));
    [r,v,r_p,v_p,a_p] = getrv(segs,nodes,whichPrimary);
    r_rel = r-r_p;
    v_rel = v-v_p;
    for i = 1:length(desPartials)
        desiredPartialCell = desPartials(i,:);
        if ~strcmp(desiredPartialCell(1),'node')
            error('Something went wrong')
        end
        desiredPartial = desiredPartialCell{3};
        switch desiredPartial %string of which partial to take
            case 'x'
                thisPartial = v_rel(1);
            case 'y'
                thisPartial = v_rel(2);
            case 'z'
                thisPartial = v_rel(3);
            case 'xd'
                thisPartial = r_rel(1);
            case 'yd'
                thisPartial = r_rel(2);
            case 'zd'   
                thisPartial = r_rel(3);
            case 't0'
                thisPartial = -dot(v_p,v_rel) - dot(a_p,r_rel);
            case {'mass','dt','SOCalpha0',...
                    'SOCalphadot','SOCbeta0','SOCbetadot'}
                thisPartial = 0;
            otherwise
                error(['Partial "' desiredPartial '" does not exist for Apse constraint'])
        end
        partials(i) = thisPartial;
    end    
end
%% utility
function [r,v,r_p,v_p,a_p] = getrv(segs,nodes,whichPrimary)
    %segs should be empty, should only have 1 node
    if ~isempty(segs) || length(nodes) ~= 1
        error('incorrect segments or nodes input to Apse')
    end
    if strcmp(nodes{1}.system_model.frame,'B1centP1P2rot')
        state = nodes{1}.getStateByIndex(1);
        r = state(1:3);
        v = state(4:6);
        % dot product between relative position and relative velociy vectors
        % F = (nodePos - primaryPos) * (nodeVel - primaryVel)
        mu = nodes{1}.system_model.char.mu;
        switch whichPrimary
            case 1
                r_p = [-mu;0;0];
                v_p = [0;0;0];
                a_p = [0;0;0];
            case 2
                r_p = [1-mu;0;0];
                v_p = [0;0;0];
                a_p = [0;0;0];
            otherwise
                error('Not a valid primary')
        end
    elseif strcmp(nodes{1}.system_model.frame,'J2000')
        if whichPrimary ~= 1
            error('apse constraint in J2000 frame must be with respect to central body')
        end
        state = nodes{1}.getStateByIndex(1);
        r = state(1:3);
        v = state(4:6);
        r_p = [0;0;0];
        v_p = [0;0;0];
        a_p = [0;0;0];
    else
        error('unknown frame')
    end
end