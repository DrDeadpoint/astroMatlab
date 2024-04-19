function con = con_SB1angle(node,desAngle)
desAngle = desAngle.change_unit('rad');
desAngle = wrapTo2Pi(desAngle.value);
whichNode = node.node_index;
F = @(segs,nodes) conF(segs,nodes,desAngle);
partials = @(segs,nodes,desPartials) calc_partials(segs,nodes,desPartials);
conname = ['Sun-B1 Angle Constraint: node ' num2str(whichNode) ...
    ', direction ' num2str(rad2deg(desAngle)) ' deg'];
con = c_constraint(conname,F,partials,[],whichNode);
end
%% constraint evaluation
function F = conF(segs,nodes,desAngle)
    [pos,sunangle] = getPos(segs,nodes);
    beta = wrapTo2Pi(atan2(pos(2),pos(1))); %angle in earth-moon frame
    alpha = wrapTo2Pi(beta - (sunangle + pi));
    F = alpha - desAngle;
end

%% partials
function partials = calc_partials(segs,nodes,desPartials)
    partials = nan(1,length(desPartials));
    [pos,~] = getPos(segs,nodes);
    node = nodes{1};
    omega_4 = node.system_model.char.angVelP4;
    x = pos(1);
    y = pos(2);
    % equation: alpha = atan2(y/x) - (omega_4 * t +pi)
    for i = 1:length(desPartials)
        desiredPartialCell = desPartials(i,:);
        if ~strcmp(desiredPartialCell(1),'node')
            error('Something went wrong')
        end
        desiredPartial = desiredPartialCell{3};
        switch desiredPartial %string of which partial to take
            case 'x'
                thisPartial = -y/(x^2+y^2);
            case 'y'
                thisPartial = x/(x^2+y^2);
            case 't0'
                thisPartial = -omega_4;
            case {'z','xd','yd','zd','mass','dt','SOCalpha0',...
                    'SOCalphadot','SOCbeta0','SOCbetadot'}
                thisPartial = 0;
            otherwise
                error(['Partial "' desiredPartial '" does not exist for Primary Distance constraint'])
        end
        partials(i) = thisPartial;
    end    
end
%% utility
function [pos,sunangle] = getPos(segs,nodes)
    %segs should be empty, should only have 1 node
    if ~isempty(segs) || length(nodes) ~= 1
        error('incorrect segments or nodes input to Position Direction Constraint')
    end
    node = nodes{1};
    pos = node.getInitPos().value();
    initSunAngle = node.system_model.char.theta0P4;
    angVelP4 = node.system_model.char.angVelP4;
    time = node.time.value; % EM nondim
    sunangle = wrapTo2Pi(initSunAngle.value + angVelP4.value*time);
end