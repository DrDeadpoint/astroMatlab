function con = con_inequality(con_equality,lessOrGreater)
    conname = ['Inequality: ' lessOrGreater '; ' con_equality.name];
    partials = @(segs,nodes,desPartials) calc_partials(segs,nodes,desPartials,con_equality);
    ows = con_equality.which_segments;
    own = con_equality.which_nodes;
    con_equality.which_segments = 1:length(ows);
    con_equality.which_nodes = 1:length(own);
    F = @(segs,nodes) conF(segs,nodes,con_equality,lessOrGreater);
    con = c_constraint(conname,F,partials,ows,own,con_equality);
end

function partials = calc_partials(segs,nodes,desPartials,con_equality)
    partials = nan(1,length(desPartials));
    for i = 1:length(desPartials)
%         desiredPartialCell = desPartials(i,:);
%         desiredPartial = desiredPartialCell{3};
%         switch desiredPartial
%             case 'slack'
%                 
%             otherwise
%                 partials(i) = con_equality.partials(segs,nodes,desPartials);
%         end
%         partials(i) = 0;
        partials(i) = nan;
    end    
end

function F = conF(segs,nodes,con_equality,lessOrGreater)
    g = con_equality.eval_f(segs,nodes);
    switch lessOrGreater
        case 'less'
            s = -1;
        case 'greater'
            s = +1;
        otherwise
            error('Inequality constraint must use args: ''less'' or ''greater''.')
    end
    if sign(g) == sign(s)
        F = 0;
    else
        F = g; %1
    end
end

