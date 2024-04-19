function con = con_elementConstraint(node,elemName,desElem)
whichNode = node.node_index;
F = @(segs,nodes) conF(segs,nodes,elemName,desElem);
partials = @(segs,nodes,desPartials) calc_partials(segs,nodes,desPartials,elemName);
conname = ['State Constraint: node ' num2str(whichNode) ...
    ', ' elemName ', desState ' num2str(desElem)];
con = c_constraint(conname,F,partials,[],whichNode);
end
%% constraint evaluation
function F = conF(segs,nodes,elemName,desElem)
    s = getElem(segs,nodes,elemName);
    F = s - desElem; 
end

%% partials
function partials = calc_partials(segs,nodes,desPartials,elemName)
    partials = nan(1,size(desPartials,1));
end
%% utility
function s = getElem(segs,nodes,elemName)
    %segs should be empty, should only have 1 node
    if ~isempty(segs) || length(nodes) ~= 1
        error('incorrect segments or nodes input to State Constraint')
    end
    if strcmp(nodes{1}.system_model.frame(1),'B')
        nodes{1} = nodes{1}.changeFrame('P2centinert');
    end
    elem = eci_elem(nodes{1});
    s = elem.(elemName);
end