clear
close all
Earth = c_body('earth');
Moon = c_body('Moon');
Sun = c_body('Sun');


descNodeP4 = c_dim_quant(1,'rad');
incP4 = c_dim_quant(0,'rad');
theta0P4 = c_dim_quant(1,'rad');
sysModel4BP = c_system_model('BCR4BP','B1centP1P2rot',Earth,Moon,Sun,...
    'descNodeP4',descNodeP4,'incP4',incP4,'theta0P4',theta0P4);
T0 = c_dim_quant([0], 'nd_t');
TF = c_dim_quant([1], 'nd_t');
time = [T0 TF];
pos = c_dim_quant([0.5; 0; 0.05], 'nd_l');
vel = c_dim_quant([.2; .4; -.05], 'nd_v');
traj = c_traj('testTraj',time,pos,vel,sysModel4BP);
traj = traj.prop();
FCs = traj.getStateByIndex(-1);
FCs = FCs(1:6);
% traj.plot;
desFCs = FCs - [0.2;0;0;0;0;0];
% plot3(desFCs(1),desFCs(2),desFCs(3),'r*','markers',15)

% close all
%%

seg = c_segment('testSeg',time,pos,vel,sysModel4BP);
seg.free_vars = {};
seg.seg_index = 1;
seg.node_indices = [1, 2];
segments = {seg};

node1 = c_node('testNode1',T0,pos,vel,sysModel4BP);
node1.free_vars = {'xd';'yd';'zd'};
node1.node_index = 1;

pos_des = pos;
pos_des.value = desFCs(1:3);
node2 = c_node('testNode2',TF,pos_des,vel,sysModel4BP);
node2.free_vars = {'x';'y';'z';'xd';'yd';'zd'};
node2.node_index = 2;

nodes = {node1, node2};

constraints = {};
for i = 1:6
    constraints{end+1} = con_stateContinuity(seg,i); %all states continuous
end
for i = 1:3
    constraints{end+1} = con_stateConstraint(node2,i,desFCs(i));
end

param = c_problem_parameter(segments,nodes,constraints);
problemSetup = c_problem_setup(param);
soln = newtonRaphsonMethod(problemSetup);

f = plotProblemIterations(soln,problemSetup);