clear
close all
Earth = c_body('earth');
Moon = c_body('Moon');
Sun = c_body('Sun');

%% three body problem
sysModel3BP = c_system_model('CR3BP','B1centP1P2rot',Earth,Moon);
time = c_dim_quant([0 10], 'day');
pos = c_dim_quant([0.5;0;0.1], 'nd_l');
vel = c_dim_quant([0;0.4;0.1], 'nd_v');
traj = c_traj('testTraj',time,pos,vel,sysModel3BP);
traj = traj.prop;
traj.plot;

%% two body problem
sysModel2BP = c_system_model('2BP','P1centinert',Earth);
period = c_dim_quant(95,'min');
ecc = c_dim_quant(0.001,'');
inc = c_dim_quant(35,'deg');
raan = c_dim_quant(1,'rad');
argP = c_dim_quant(10,'deg');
trueAnom = c_dim_quant(pi,'rad');
elem = c_elem(Earth.mu,'period',period,'ecc',ecc,'inc',inc,...
    'raan',raan,'argP',argP,'trueAnom',trueAnom);
traj = c_elem.conv2traj();
traj = traj.prop;
traj.plot;

%% four body problem
close all
descNodeP4 = c_dim_quant(1,'rad');
incP4 = c_dim_quant(0,'rad');
theta0P4 = c_dim_quant(1,'rad');
sysModel4BP = c_system_model('BCR4BP','B1centP1P2rot',Earth,Moon,Sun,...
    'descNodeP4',descNodeP4,'incP4',incP4,'theta0P4',theta0P4);
time = c_dim_quant([0 5], 'nd_t');
pos = c_dim_quant([0.5;0.2;0.2], 'nd_l');
vel = c_dim_quant([0;0.9;0.5], 'nd_v');
traj = c_traj('testTraj',time,pos,vel,sysModel4BP);
traj = traj.prop;
traj.plot;
trueFC = [-0.361217999451362
           0.758649995090825
           0.289332514214419
          -0.327608246402497
           0.431659912817940
           0.201051314393174];
thisFC = traj.getStateByIndex(traj.getLength);
%trueFC - thisFC(1:6)