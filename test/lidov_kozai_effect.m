clear
close all
Earth = c_body('earth');
Moon = c_body('Moon');
%lidov-kozai effect leads to a periodic exchange between eccentricity and
%inclination for highly inclined orbits

%% 30 RE circular orbit, 90 degree inc
sysModel2BP = c_system_model('2BP','P1centinert',Earth);
rad = 30 * Earth.radius; %km
vel = sqrt(sysModel2BP.char.mu / rad); %km/s

%% three body problem
sysModel3BP = c_system_model('CR3BP','B1centP1P2rot',Earth,Moon);
time = c_dim_quant([0 5], 'year');
inc = 95; %deg
pos = c_dim_quant([0;rad/sysModel3BP.char.lstar*cosd(inc);rad/sysModel3BP.char.lstar*sind(inc)], 'nd_l');
vel = c_dim_quant([vel/sysModel3BP.char.lstar*sysModel3BP.char.tstar;0;0], 'nd_v');
traj = c_traj('LKE',time,pos,vel,sysModel3BP);
traj = traj.prop;
traj = traj.changeFrame('P1centinert');
traj.plotting.lineWidth = 0.1;
traj.plot
view(-18,5)
%%
elems = osculating_elements(traj,'P1');
%%
figure
