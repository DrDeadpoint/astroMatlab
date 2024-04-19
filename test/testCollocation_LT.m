clear
close all
Earth = c_body('earth');
Moon = c_body('Moon');
model = c_system_model('CR3BP','B1centP1P2rot',Earth,Moon);
Tmax = c_dim_quant(100,'mN');
Isp = c_dim_quant(2000,'sec');
M0 = c_dim_quant(500,'kg');
spacecraft = c_spacecraft(Tmax,Isp,M0,model,0.2);
controlLaw = c_control_law('fixed_dir','B1centP1P2rot',[1;0;0]);

dps = defaultPlotting;
dps.lineColor = colour('blue'); dps.LTplotmult = 0.1;
tspan = [0 1];
traj = c_traj('traj',c_dim_quant(tspan,'nd_t'),...
    c_dim_quant([0.7; 0; 0],'nd_l'),...
    c_dim_quant([-0.2; 0.5; 0],'nd_v'),model,...
    'spacecraft',spacecraft,'plotting',dps,...
    'mass',M0,'control_law',controlLaw);
traj = traj.prop;
ax = traj.plot();
legend
makeGray(gcf);

desInitNode = traj.getStateByIndex(1);
desInitNode = desInitNode(1:6);
finalNode = traj.getFinalPos().value;
desFinalNode = [finalNode.*1.01; nan; nan; nan];

%% run collocation
collTraj = traj;
problem = c_coll_problem(collTraj,'order',7,...
    'desFinalNode',desFinalNode,'desInitNode',desInitNode,...
    'mesh','5_equalTime');
options = optimoptions('fmincon','MaxIterations',500,...
                        'MaxFunctionEvaluations',1e6,...
                        'Display','iter-detailed',...
                        'OptimalityTolerance',1e-9,...
                        'ConstraintTolerance',1e-9,...
                        "EnableFeasibilityMode",true);%,...
                        %"SubproblemAlgorithm","cg");
[soln,axdebug] = problem.solve('options',options);

%%
close all
[f,newAx] = ffigure;
makeGray(f);
Moon.plot(newAx,model,'nd_l');
axis equal
soln.plot('ax',newAx,'DisplayName','Result','plotLT',true);
makeGray(gcf);
plotStates(newAx,finalNode,'*','DisplayName','Orig Final Node', 'Color', colour('o'))
plotStates(newAx,desFinalNode,'*','DisplayName','Des Final Node','Color',colour('r'))
plotStates(newAx,desInitNode,'*','DisplayName','Des Init Node','Color',colour('g'))
legend(newAx)