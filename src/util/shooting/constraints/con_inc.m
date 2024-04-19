function con = con_inc(node,desInc)
whichNode = node.node_index;
F = @(segs,nodes) conF(segs,nodes,desInc);
partials = @(segs,nodes,desPartials) calc_partials(segs,nodes,desPartials);
conname = ['inc Constraint: node ' num2str(whichNode) ', desInc ' num2str(desInc)];
con = c_constraint(conname,F,partials,[],whichNode);
end
%% constraint evaluation
function F = conF(segs,nodes,desInc)
    s = getInc(segs,nodes);
    F = s - desInc; 
end

%% partials
function partials = calc_partials(segs,nodes,desPartials)
    partials = nan(1,size(desPartials,1));
end
%% utility
function s = getInc(segs,nodes)
    %segs should be empty, should only have 1 node
    if ~isempty(segs) || length(nodes) ~= 1
        error('incorrect segments or nodes input to State Constraint')
    end
    if strcmp(nodes{1}.system_model.frame(1),'B')
        nodes{1} = nodes{1}.changeFrame('P2centinert');
    end
    mu = nodes{1}.system_model.char.mu;
    primaryLocation = [1-mu;0;0];
    r = nodes{1}.getInitPos.value; v = nodes{1}.getInitVel.value;
    r_rel = r - primaryLocation;
    r_rel = r_rel./norm(r_rel); v = v./norm(v);
    h = cross(r_rel,v); %angular momentum of orbit relative to primary
    inc = acos(dot(h,[0;0;1])); %angle between ang mom and z-axis
    s = inc;
end