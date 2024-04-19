function [newSegments, newNodes, newConstraints] = discretizeSetup(segments,nodes,constraints,seg_ind,num_segs,method)
% newProblemSetup = discretizeSetup(segments,nodes,constraints,seg_ind,num_segs,method)
% discretizes a segment, updates all constraints to new references
%   seg_ind: index of segment to be discretized
%   num_segs: how many segments to break it into
%   method: what method to use to discretize
%       'time'
%       'int'egration step
%       'arc'length
if num_segs == 1
    newSegments = segments;
    newNodes = nodes;
    newConstraints = constraints;
else
% extract segment and its nodes
seg = segments{seg_ind};
node_inds = seg.node_indices; %generally [4,5] or [5,4]
node1 = nodes{node_inds(1)};
node2 = nodes{node_inds(2)};

if abs(diff(node_inds)) ~=1
    error('this function relies on nodes''s being adjacent')
end

% create new nodes for segment
new_node_inds = min(node_inds):(min(node_inds) + num_segs); %[4 5 6 7]
if node_inds(2) < node_inds(1)
    new_node_inds = flip(new_node_inds); %[7 6 5 4]
end
seg.pos = node1.pos; seg.vel = node1.vel; %ensure it starts at node1
seg.time.value = [node1.time.value, node2.time.value]; %ensure it has correct tspan
seg = seg.prop;

switch method
    case 'time'
        times = linspace(node1.time.value, node2.time.value, num_segs+1);
    case 'int'
        inds = round(linspace(1,length(seg.time.value),num_segs+1));
        times = seg.time.value(inds);
    case 'arc'
        error('not implemented')
    otherwise
        error('unknown discretization method')
end
new_positions = cell(1,num_segs+1);
new_velocities = cell(1,num_segs+1);
new_positions{1} = node1.pos;
new_velocities{1} = node1.vel;
new_nodes = cell(1,num_segs+1);
node1.node_index = new_node_inds(1);
new_nodes{1} = node1;
for i = 2:length(times)-1
    tempSeg = seg;
    tempSeg.time.value = [node1.time.value, times(i)];
    tempSeg.pos = node1.pos; tempSeg.vel = node1.vel;
    tempSeg = tempSeg.prop;
    node_pos = tempSeg.getFinalPos();
    node_vel = tempSeg.getFinalVel();
    
    new_node = node1;
    new_node.time.value = times(i);
    new_node.pos = node_pos;
    new_node.vel = node_vel;
    new_node.node_index = new_node_inds(i);
    new_node.name = ['Interior node ' num2str(i-1) ' for segment: ' seg.name];
    new_node.free_vars = {'x';'y';'z';'xd';'yd';'zd';'mass'};
    new_nodes{i} = new_node;
    new_positions{i} = node_pos;
    new_velocities{i} = node_vel;
end
node2.node_index = new_node_inds(num_segs+1);
new_positions{num_segs+1} = node2.pos;
new_velocities{num_segs+1} = node2.vel; %should be end of this cell
new_nodes{num_segs+1} = node2;
if node_inds(2) < node_inds(1)
    new_nodes = flip(new_nodes); %{4 5 6 7}
end
beforeInds = 1:(min(node_inds)-1);
afterInds = (max(node_inds)+1):length(nodes);
for i = 1:length(afterInds)
    node = nodes{afterInds(i)};
    node.node_index = node.node_index + num_segs-1;
    nodes{afterInds(i)} = node;
end
newNodes = [nodes(beforeInds), new_nodes, nodes(afterInds)];

% create new segments between nodes
new_seg_inds = seg_ind:(seg_ind + num_segs - 1);
if node_inds(2) < node_inds(1)
    new_seg_inds = flip(new_seg_inds); %[6 5 4]
end
new_segments = cell(1,length(new_seg_inds));
for i = 1:length(new_seg_inds)
    new_seg = seg;
    new_seg.time.value = times(i:i+1);
    new_seg.pos = new_positions{i};
    new_seg.vel = new_velocities{i};
    new_seg.node_indices = new_node_inds(i:i+1);
    new_seg.seg_index = new_seg_inds(i);
    new_seg.name = [seg.name ', section ' num2str(i) '/' num2str(length(new_seg_inds))];
    new_segments{i} = new_seg;
end
if node_inds(2) < node_inds(1)
    new_segments = flip(new_segments); %{4 5 6}
end
beforeInds = 1:seg_ind-1;
afterInds = seg_ind+1:length(segments);
for i = 1:length(afterInds)
    seg = segments{afterInds(i)};
    seg.node_indices = seg.node_indices + num_segs-1;
    seg.seg_index = seg.seg_index + num_segs-1;
    segments{afterInds(i)} = seg;
end
newSegments = [segments(beforeInds), new_segments, segments(afterInds)];

% redefine old constraints, add continuity constraints to interior segs
% watch for seg_ind, node_inds

for i = 1:length(constraints)
    icon = constraints{i};
    %segment constraints are moved to final segment
    if contains(icon.name,'seg')
        if strcmp(icon.name(1:14),'Time of Flight')
            error('need to implement update to TOF constraint')
        end
        if isequal(icon.which_segments, seg_ind)
            %fix name
            
            %fix segment index
            icon.which_segments = max(new_seg_inds);
            %fix node index
            if node_inds(2) < node_inds(1)
                icon.which_nodes = [max(new_node_inds), max(new_node_inds)-1];
            else
                icon.which_nodes = [max(new_node_inds)-1, max(new_node_inds)];
            end
        else
            if seg_ind < icon.which_segments
                icon.which_segments = icon.which_segments + num_segs-1;
                icon.which_nodes = icon.which_nodes + num_segs-1;
            end
        end
    else %node constraint, needs to be updated
        if isequal(icon.which_nodes, max(node_inds))
            %fix name

            %fix node indices
            icon.which_nodes = max(new_node_inds);
        elseif isequal(icon.which_nodes, min(node_inds)) %minimum node index hasn't changed
            %do nothing
        elseif max(node_inds) < icon.which_nodes
            icon.which_nodes = icon.which_nodes + num_segs-1;
        end
    end
    constraints{i} = icon;
end
%add continuity constraints
new_seg_inds = sort(new_seg_inds);
contConstraints = {};
for i = 1:length(new_seg_inds)-1
    for j = 1:7
        contConstraints{end+1} = con_stateContinuity(newSegments{new_seg_inds(i)},j);
    end
    contConstraints{end+1} = con_timeContinuity(newSegments{new_seg_inds(i)});
end
newConstraints = [constraints(:); contConstraints(:)];
end
end