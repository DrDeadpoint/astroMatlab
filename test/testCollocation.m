clear
close all
Earth = c_body('earth');
Moon = c_body('Moon');
model = c_system_model('CR3BP','B1centP1P2rot',Earth,Moon);
Tmax = c_dim_quant(10,'mN');
Isp = c_dim_quant(2000,'sec');
M0 = c_dim_quant(500,'kg');
spacecraft = c_spacecraft(Tmax,Isp,M0,model,0.01);
controlLaw = c_control_law('fixed_dir','B1centP1P2rot',[1;0;0]);

%
load('L1L2Lyap.mat');
[L1lyap,L2lyap,L1man,L2man] = deal(plotObjects{:});
% plot Lyapunov orbits
dpsL1 = defaultPlotting;
dpsL1.lineStyle = '--'; dpsL1.lineWidth = 1;
L1traj = c_traj('L1 Lyapunov',c_dim_quant(L1lyap.time([1,end]),'nd_t'),...
    c_dim_quant(L1lyap.state(:,1:3)','nd_l'),...
    c_dim_quant(L1lyap.state(:,4:6)','nd_v'),model,...
    'plotting',dpsL1);
L1traj = L1traj.prop;

dpsL2 = defaultPlotting;
dpsL2.lineStyle = '-.'; dpsL2.lineWidth = 1;
L2traj = c_traj('L2 Lyapunov',c_dim_quant(L2lyap.time([1,end]),'nd_t'),...
    c_dim_quant(L2lyap.state(:,1:3)','nd_l'),...
    c_dim_quant(L2lyap.state(:,4:6)','nd_v'),model,...
    'plotting',dpsL2);
L2traj = L2traj.prop;
ax = L1traj.plot;
L2traj.plot('ax',ax);

% plot manifold arcs
dps2 = defaultPlotting;
dps2.lineColor = colour('red');
L2mantraj = c_traj('L2 Stable Manifold',c_dim_quant(L2man.time([1,end]),'nd_t'),...
    c_dim_quant(L2man.state(1,1:3)','nd_l'),...
    c_dim_quant(L2man.state(1,4:6)','nd_v'),model,...
    'spacecraft',spacecraft,'mass',M0,...
    'plotting',dps2,'control_law',controlLaw);
L2mantraj = L2mantraj.prop;
L2mantraj.plot('ax',ax);

dps1 = defaultPlotting;
dps1.lineColor = colour('blue');
L1mantspan = L1man.time([end,1]) + -L1man.time(end) + L2man.time(end);
L1mantraj = c_traj('L1 Stable Manifold',c_dim_quant(L1mantspan,'nd_t'),...
    c_dim_quant(L1man.state(end,1:3)','nd_l'),...
    c_dim_quant(L1man.state(end,4:6)','nd_v'),model,...
    'spacecraft',spacecraft,'plotting',dps1,...
    'mass',L2mantraj.getFinalMass(),'control_law',controlLaw);
L1mantraj = L1mantraj.prop;
L1mantraj.plot('ax',ax);
legend
makeGray(gcf);

desInitNode = L2man.state(1,:)';
desFinalNode = L1man.state(1,:)';

%% run collocation
collTraj = {L2mantraj, L1mantraj};
problem = c_coll_problem(collTraj,'order',7,...
    'desFinalNode',desFinalNode,'desInitNode',desInitNode,...
    'mesh','4_equalTime');
options = optimoptions('fmincon','MaxIterations',10,...
                        'MaxFunctionEvaluations',1e6,...
                        'Display','iter-detailed',...
                        'OptimalityTolerance',1e-9,...
                        'ConstraintTolerance',1e-9,...
                        "EnableFeasibilityMode",true,...
                        "SubproblemAlgorithm","cg");
[soln,axdebug] = problem.solve('options',options);

%%
close all
[f,newAx] = ffigure;
makeGray(f);
Moon.plot(newAx,model,'nd_l');
axis equal
soln.plot('ax',newAx,'DisplayName','Result');
makeGray(gcf);
plotStates(newAx,desFinalNode,'*','DisplayName','Final Node','Color',colour('r'))
plotStates(newAx,desInitNode,'*','DisplayName','Init Node','Color',colour('g'))
legend(newAx)