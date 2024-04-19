function [traj_out] = prop_ephem(traj_in,varargin)
wantSTM = true;
if nargin == 2
    wantSTM = varargin{1};
end
ephem_options = traj_in.etc.ephemOptions;
% seg = ephem_options{1};
seg = traj_in.etc.seg;
segNo = seg.segNumber;
segNos = ephem_options{2};
propDirs = ephem_options{3};
opts = ephem_options{4};
u_TTLc = ephem_options{5};
if length(ephem_options) > 5
    harm = ephem_options{6};
else
    harm = {false;[]}; %don't use harmonics
end

sysModel = traj_in.system_model;
char = sysModel.char;
EMFrame = char.EMFrame;
SEFrame = char.SEFrame;
EMFrame2 = char.EMFrame2;
Body = char.Body;
mu = char.mu;
lstar = char.lstar.value;
tstar = char.tstar.value;
JDFix = char.JDFix;

traj_in.pos = traj_in.pos.change_unit('nd_l',sysModel);
traj_in.vel = traj_in.vel.change_unit('nd_v',sysModel);

Y0 = [traj_in.pos.value(:,1); traj_in.vel.value(:,1)]';
if length(traj_in.time.value) ~= 2
    error('must use a time vector of length 2')
end
traj_in.time = traj_in.time.change_unit('nd_t',sysModel);
t0 = traj_in.time.value(1);
dt = traj_in.time.value(2) - t0;

if dt == 0
   error('should have skipped this one')
end
segNos(end+1) = seg.segNumber;
if dt > 0
    propDirs = [propDirs 1];
else
    propDirs = [propDirs -1];
end

%get propagator frame
if ~strcmp(seg.prop.frame,'J2000')
   error(['Unprecedented frame on segment ' segNo]) 
end
switch seg.prop.centralBody
    case 'EARTH'
        propFrame = EMFrame;
    case 'MOON'
        propFrame = EMFrame2;
    case 'SUN'
        propFrame = SEFrame;
        error('SE frame not implemented');
    otherwise
        error(['Unknown propagation frame on segment ' segNo])
end

if ~strcmp(seg.IC.stateVars{1},'Rx (km)')
   error(['Segment ' num2str(segNo) ...
       ' does not have state in position and velocity space'])
end
tVec = [t0, t0+dt];

% propagation
m0seg1 = traj_in.low_thrust.spacecraft.M0; %kg
m0 = traj_in.low_thrust.mass;
use_lowthrust = false;
use_harmonics = harm{1};
%initialize
uCoef = [];
SC = [];
manFrame = [];
manSys = [];
harmonics_coefs = harm{2};
%declare values
Y0 = [Y0, m0.value(1)];
if ~isnan(seg.maneuver.frame) % a maneuver takes place
    use_lowthrust = true;
    claw = traj_in.low_thrust.control_law;
    spacecraft = traj_in.low_thrust.spacecraft;
    uCoef = [claw.coeffs.alpha0; claw.coeffs.alphadot;...
        claw.coeffs.beta0; claw.coeffs.betadot];
    thrust = spacecraft.Tmax.change_unit('N');
    thrust = thrust.value; 
    thrust = thrust * spacecraft.throttle;
    isp = spacecraft.Isp.value;
    SC.m0D = m0seg1.value;
    SC.thrustMaxD = thrust;
    SC.ispD = isp;
    SC = setSpacecraft(traj_in.etc.system,SC);
    SC = SC{1};
    % sometimes maneuver is in J2000 frame, sometimes in 2B_rotating
    switch claw.frame{1}
        case 'J2000'
            switch claw.frame{2}
                case 'IJK'
                    manFrame = 'J2000';
                case 'VUW'
                    if strcmp(seg.maneuver.centralBody,propFrame.centralBody)
                        manFrame = 'J2000 VUW';
                    else
                        error('need new prop code')
                    end
            end
        case '2B_ROTATING'
            manFrame = '2B rotating VUW';
            manMainBody = seg.maneuver.mainBody;
            manCentralBody = claw.frame{3};
            switch manMainBody
                case 'SUN'
                    manSys = SEFrame;
                case 'EARTH'
                    switch manCentralBody
                        case 'EARTH'
                            manSys = EMFrame;
                        case 'MOON'
                            manSys = EMFrame2;
                    end
            end
    end
end
if wantSTM
    STM_0 = eye(7);
    STM_0 = reshape(STM_0,49,1);
    ddYdT = zeros(6,1);
    Y0 = [Y0'; STM_0; ddYdT];
end
[T, Y] = ode113(@(t,y) ephemeris(t, y, JDFix, propFrame, Body,...
    'lowthrust',use_lowthrust,'uCoef',uCoef,'t0',t0,'spacecraft',SC,...
    'manFrame',manFrame,'manSystem',manSys,'harmonics',use_harmonics,...
    'harmonics_coefs',harmonics_coefs,'stm',wantSTM), ...
        tVec, Y0, opts);
Y = transpose(Y);
if wantSTM
    STM_f = Y(8:56,end);
    STM_f = reshape(STM_f,7,7);
    ddydt = Y(57:62,end);
    STM_f = [STM_f, [ddydt;0]]; %append to STM (STM-like term)
else
    STM_f = [];
end
Y = Y(1:7,:)';

uIMat = zeros(length(T), 3);
uMat_VUW = uIMat;
uMat_R = uIMat;
for i = 1:length(T)
    [~,u_I,u_VUW,u_R] = ephemeris(T(i), Y(i,:)', JDFix, propFrame, Body,...
        'lowthrust',use_lowthrust,'uCoef',uCoef,'t0',t0,'spacecraft',SC,...
        'manFrame',manFrame,'manSystem',manSys,'harmonics',use_harmonics,...
        'harmonics_coefs',harmonics_coefs);
    uIMat(i, :) = u_I;
    uMat_VUW(i,:) = u_VUW;
    uMat_R(i,:) = u_R;
end
u_TTLc{end+1} = {uIMat; uMat_VUW; uMat_R};

traj_out = traj_in;
traj_out.pos = c_dim_quant(Y(:,1:3)','nd_l');
traj_out.vel = c_dim_quant(Y(:,4:6)','nd_v');
traj_out.time = c_dim_quant(T,'nd_t');
traj_out.low_thrust.mass = c_dim_quant(Y(:,7)','nd_m');
STMkey = {'x', 'y', 'z', 'xd', 'yd', 'zd', 'mass', 't0'};
traj_out.stm = c_stm(STM_f,STMkey);

optionsOut = ephem_options;
optionsOut{2} = segNos;
optionsOut{3} = propDirs;
optionsOut{5} = u_TTLc;
traj_out.etc.ephemOptions = optionsOut;
end