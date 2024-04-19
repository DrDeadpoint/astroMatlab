clear
close all

Earth = c_body('earth');
Moon = c_body('Moon');
claw = c_control_law('fixed_dir','B1centP1P2rot',[1 0 0]);

%% three body problem
model = c_system_model('CR3BP','B1centP1P2rot',Earth,Moon);
tspan = c_dim_quant([0 1],'nd_t');
traj = c_traj('test seg',tspan,...
    c_dim_quant([0.7;0;0],'nd_l'),...
    c_dim_quant([-0.2;0.5;0],'nd_v'),model);
traj.low_thrust.control_law = claw;
traj = traj.prop();
ax = traj.plot();
makeGray(gcf);

desInitNode = [traj.getInitPos().value; nan;nan;nan];
desFinalNode = [0.72;0.1;0;nan;nan;nan];


%% run collocation
problem = c_coll_problem(traj,'order',5,'isBallistic',true,...
    'desFinalNode',desFinalNode,'desInitNode',desInitNode,...
    'mesh','2_equalTime');
options = optimoptions('fmincon','MaxIterations',1000,...
                        'Display','iter-detailed',...
                        'OptimalityTolerance',1e-9,...
                        'ConstraintTolerance',1e-9);
[soln,axdebug] = problem.solve('options',options);

%%
% problem.plot(soln,'ax',axdebug);
soln.plot('ax',axdebug);
legend
legend(axdebug)
makeGray(gcf);
plot3(axdebug,desFinalNode(1),desFinalNode(2),desFinalNode(3),'r*','DisplayName','target pos')
% traj.pos.value = soln.soln(1:3);
% traj.vel.value = soln.soln(4:6);
traj.pos.value = soln.converged_opt_vars(1:3);
traj.vel.value = soln.converged_opt_vars(4:6);
traj.time = tspan;
traj = traj.prop;
traj.plotting.lineColor = colour('g');
traj.plotting.lineWidth = 1;
traj.plot('ax',axdebug)

%% Compare to shooting



