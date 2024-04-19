clear
close all
Earth = c_body('earth');
Moon = c_body('Moon');
claw = c_control_law('fixed_dir','P1centinert',[1 0 0]);

%% two body problem
model = c_system_model('2BP','P1centinert',Earth);
sma = c_dim_quant(22000,'km');
ecc = c_dim_quant(0.7,'');
inc = c_dim_quant(0,'deg');
raan = c_dim_quant(0,'rad');
argP = c_dim_quant(0,'deg');
trueAnom = c_dim_quant(0,'rad');
elem = c_elem(Earth.mu,'sma',sma,'ecc',ecc,'inc',inc,...
    'raan',raan,'argP',argP,'trueAnom',trueAnom);
traj = eci_elem(elem,model);
traj.time.value(2) = traj.time.value(2)/1.05;
tspan = traj.time.value;
traj.low_thrust.control_law = claw;
traj = traj.prop;
ax = traj.plot;
makeGray(ax);

%%
% split trajectory and introduce some errors in position and velocity
numSegs = 1;
nodelocs = floor(linspace(1,length(traj.time.value),numSegs+1));
collTraj = {};
ra = (rand(3,1)*0.1 + 0.95);
%ra =[
%    0.964986544247797
%    1.015960525290831
%    1.001859494251054];
for i = 1:numSegs
    ind1 = nodelocs(i);
    ind2 = nodelocs(i+1);
    tempTraj = traj.extractSubTraj(ind1,ind2);
    tempTraj.time.value = tempTraj.time.value([1,end]);
    if i > 1
        tempTraj.pos.value = tempTraj.pos.value .* ra;
        tempTraj.vel.value = tempTraj.vel.value .* ra;
    end
    tempTraj.plotting.lineColor = colour('g');
    tempTraj.plot('ax',ax);
    collTraj{end+1} = tempTraj;
end


%% run collocation
problem = c_coll_problem(collTraj,'order',7,'isBallistic',true);
options = optimoptions('fmincon','MaxIterations',1000,...
                        'Display','iter-detailed',...
                        'OptimalityTolerance',1e-9,...
                        'ConstraintTolerance',1e-9);
[soln,axdebug] = problem.solve('options',options);

%%
soln.plot('ax',axdebug,'DisplayName','Result');
legend
makeGray(gcf);
traj.pos.value = soln.converged_opt_vars(1:3);
traj.vel.value = soln.converged_opt_vars(4:6);
% soln.plot('ax',axdebug,'DisplayName','Result');
% legend
% makeGray(gcf);
% traj.pos.value = soln.converged_opt_vars(1:3);
% traj.vel.value = soln.converged_opt_vars(4:6);
traj.time.value = tspan;
traj = traj.prop;
traj.plotting.lineColor = colour('g');
traj.plotting.lineWidth = 1;
traj.plot('ax',axdebug);