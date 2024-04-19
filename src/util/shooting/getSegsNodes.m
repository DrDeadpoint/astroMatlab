function [segments,nodes] = getSegsNodes(x,param,fcnOptions,wantSTM,varargin)
    if ~iscolumn(x)
        error('input x must be column vector')
    end
    if strcmp(wantSTM,'noSTM')
        wantSTM = false;
    elseif strcmp(wantSTM, 'yesSTM')
        wantSTM = true;
    else
        error('unknown STM argument')
    end
    useKepler = true;
    if ~isempty(varargin)
        useKepler = varargin{1};
    end
    segments = param.segments;
    nodes = param.nodes;
    %go through and replace with free variables for nodes first
    %available free variables:
        % t0
        % x,y,z,xd,yd,zd
        % mass
    for i = 1:length(nodes)
        node = nodes{i};
        node.time.value = node.checkIfFree('t0',x,node.time.value(1));
        pos = node.getInitPos;
        vel = node.getInitVel;
        xpos = node.checkIfFree('x',x,pos.value(1,1));
        ypos = node.checkIfFree('y',x,pos.value(2,1));
        zpos = node.checkIfFree('z',x,pos.value(3,1));
        xvel = node.checkIfFree('xd',x,vel.value(1,1));
        yvel = node.checkIfFree('yd',x,vel.value(2,1));
        zvel = node.checkIfFree('zd',x,vel.value(3,1));
        node.pos.value = [xpos;ypos;zpos];
        node.vel.value = [xvel;yvel;zvel];
        if ~isempty(node.low_thrust.spacecraft)
            node.low_thrust.mass.value = node.checkIfFree('mass',x,...
                node.low_thrust.mass.value(1));
        end
        nodes{i} = node;
    end
    %now do segments
    %available free variables:
        % dt
        % thrust
        % claw coefficients
    for i = 1:length(segments)
        seg = segments{i};
        startNode = nodes{seg.node_indices(1)};
        endNode = nodes{seg.node_indices(2)};
        T0 = startNode.time.value;
%         DT = seg.checkIfFree('dt',x,seg.time.value(end) - seg.time.value(1));
%         seg.time.value = [T0, T0 + DT];
        TF = endNode.time.value;
        seg.time.value = [T0 TF];
        if any(isnan(seg.time.value))
            error('Bad time input')
        end
        seg.pos = startNode.pos;
        seg.vel = startNode.vel;
        seg.low_thrust.mass = startNode.low_thrust.mass;
        if ~isempty(seg.low_thrust.spacecraft)
            SC = seg.low_thrust.spacecraft;
            thrust = SC.Tmax;
            thrust.value = seg.checkIfFree('thrust',x,SC.Tmax.value);
            SC = c_spacecraft(thrust,SC.Isp,SC.M0,seg.system_model,SC.throttle); %might have changed thrust
            seg.low_thrust.spacecraft = SC;
            claw = seg.low_thrust.control_law;
            fnames = fieldnames(claw.coeffs);
            for j = 1:length(fnames)
                coeff_j = seg.checkIfFree(fnames{j},x,claw.coeffs.(fnames{j}));
                claw.coeffs.(fnames{j}) = coeff_j;
            end
            seg.low_thrust.control_law = claw;
        end
        if ~isempty(seg.event_list)
            events = seg.event_list;
        else
            events = 'none';
        end
        seg = seg.prop('prop',fcnOptions.propMethod,...
            'events',events,'options',fcnOptions.integratorOptions,...
            'wantSTM',wantSTM,'useKepler', useKepler); %useKepler
        segments{i} = seg;
    end
end