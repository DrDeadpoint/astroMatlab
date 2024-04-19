function con = con_positionDirectionConstraint(node,desDir)
desDir = desDir./norm(desDir); %unit vector
whichNode = node.node_index;
F = @(segs,nodes) conF(segs,nodes,desDir);
partials = @(segs,nodes,desPartials) calc_partials(segs,nodes,desPartials,desDir);
conname = ['Position Unit Vector Constraint: node ' num2str(whichNode) ...
    ', direction [' num2str(desDir(1)) '; ' num2str(desDir(2)) '; ' num2str(desDir(3)) ']'];
con = c_constraint(conname,F,partials,[],whichNode);
end
%% constraint evaluation
function F = conF(segs,nodes,desDir)
    pos = getPos(segs,nodes);
    F = dot(pos,desDir)^2 - (pos(1)^2 + pos(2)^2 + pos(3)^2);
end

%% partials
function partials = calc_partials(segs,nodes,desPartials,desDir)
    partials = nan(1,length(desPartials));
    pos = getPos(segs,nodes);
    xdes = desDir(1);
    ydes = desDir(2);
    zdes = desDir(3);
    for i = 1:length(desPartials)
        desiredPartialCell = desPartials(i,:);
        if ~strcmp(desiredPartialCell(1),'node')
            error('Something went wrong')
        end
        desiredPartial = desiredPartialCell{3};
        switch desiredPartial %string of which partial to take
            case 'x'
                thisPartial = 2*dot(pos,desDir)*xdes - 2*pos(1);
%                 thisPartial = nan;
            case 'y'
                thisPartial = 2*dot(pos,desDir)*ydes - 2*pos(2);
%                 thisPartial = nan;
            case 'z'
                thisPartial = 2*dot(pos,desDir)*zdes - 2*pos(3);
%                 thisPartial = nan;
            case {'t0','xd','yd','zd','mass','dt','SOCalpha0',...
                    'SOCalphadot','SOCbeta0','SOCbetadot'}
                thisPartial = 0;
            otherwise
                error(['Partial "' desiredPartial '" does not exist for Primary Distance constraint'])
        end
        partials(i) = thisPartial;
    end    
end
%% utility
function pos = getPos(segs,nodes)
    %segs should be empty, should only have 1 node
    if ~isempty(segs) || length(nodes) ~= 1
        error('incorrect segments or nodes input to Position Direction Constraint')
    end
    node = nodes{1};
    pos = node.getInitPos().value();
end