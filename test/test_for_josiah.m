clear
close all
Earth = c_body('earth');
Moon = c_body('Moon');
Sun = c_body('Sun');

descNodeP4 = c_dim_quant(0,'rad');
incP4 = c_dim_quant(0,'rad');
theta0P4 = c_dim_quant(pi+2.383129163559179,'rad');
sysModel4BP = c_system_model('BCR4BP','B2centP4B1rot',Earth,Moon,Sun,...
    'descNodeP4',descNodeP4,'incP4',incP4,'theta0P4',theta0P4);

tof = 1.310675113987149;
tof_EM = tof*sysModel4BP.char.tstar/sysModel4BP.char.tstarP1P2;

time = c_dim_quant([0], 'nd_t');
pos = c_dim_quant([1.000036370326526; -0.000018814052405; 0], 'nd_l');
vel = c_dim_quant([0.159517740776120; 0.333836678554480; 0], 'nd_v');
traj = c_traj('traj',time,pos,vel,sysModel4BP);
traj = traj.changeFrame('B1centP1P2rot');
traj2 = traj.changeFrame('B2centP4B1rot');
traj.time = c_dim_quant([0 tof_EM], 'nd_t');
traj = traj.prop;
traj.plot;
