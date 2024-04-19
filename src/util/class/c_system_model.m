classdef c_system_model
    % sysModel = MODEL(systemDynamics,frame,bodies,varargin)
    % builds a model based on body data
    properties
        % standard data
        system_dynamics
        unit
        P1
        P2
        P4
        char
        frame
    end
    
    methods
        function sysModel = c_system_model(systemDynamics,frame,P1,P2,P4,varargin)
                % sysModel = c_system_model(systemDynamics,frame,bodies,varargin)
                % builds a model based on body data
            method = 'atd';
            if ~any(strcmp(systemDynamics,{'2BP','CR3BP','BCR4BP','EPHEMERIS','NBODY'}))
               error(['The input system dynamics **' systemDynamics '** has not yet been implemented']) 
            end
            narg = nargin - 2;
            if narg == 0
                error('Need at least one BODY')
            end
            if narg >= 1 % have a P1
                if isa(P1,'c_body')
                    sysModel.P1 = P1;
                else
                    error('Second variable argument must be of class *c_body*')
                end
                J2 = 0;
            end
            if narg >=2 % have a P2
                if isa(P2,'c_body')
                    sysModel.P2 = P2;
                else
                    error('Third variable argument must be of class *c_body*')
                end
            end
            if narg >=3 % have a P4
                if isa(P4,'c_body')
                    sysModel.P4 = P4;
                else
                    error('Fourth variable argument must be of class *c_body*')
                end
                descNodeP4 = c_dim_quant(0,'rad');
                incP4 = c_dim_quant(0,'rad');
                theta0P4 = c_dim_quant(0,'rad');
            end
            if narg > 3 % have 4BP inputs
                for i = 1:length(varargin)/2
                   switch varargin{i*2-1}
                       case 'descNodeP4'
                           descNodeP4 = varargin{i*2};
                       case 'incP4'
                           incP4 = varargin{i*2};
                       case 'theta0P4'
                           theta0P4 = varargin{i*2};
                       case 'J2'
                           J2 = varargin{i*2};
                       case 'JDFix'
                           JDFix = varargin{i*2};
                       case 'seg'
                           seg = varargin{i*2};
                       otherwise
                           error('Unknown varargin')
                   end
                end
            end
            %%
            switch systemDynamics % build the model
                case '2BP'
                    sysModel.P1 = P1;
                    sysModel.unit = P1.radius.unit; 
                    char.mu = P1.mu;
                    char.J2 = J2;
                    char.lstar = c_dim_quant(1,'km');
                    char.tstar = c_dim_quant(1,'sec');
                    char.mstar = c_dim_quant(1,'kg');
                    sysModel.char = char;
                case 'CR3BP'
                    sysModel.P1 = P1;
                    sysModel.P2 = P2;
                    [mu3, mstar, lstar, tstar] = CR3BP_data(P1,P2,method);
                    char.mu = mu3;
                    char.mstar = mstar;
                    char.lstar = lstar;
                    char.tstar = tstar;
                    sysModel.char = char;
                    sysModel.unit = 'nd';
                case 'BCR4BP'
                    [muP1P2, mstarP1P2, lstarP1P2, tstarP1P2,mP4,aP4,wP4,mstarP4B1,muP4B1,lstarP4B1,tstarP4B1] = BCR4BP_data(P1,P2,P4,theta0P4,method);
                    sysModel.P1 = P1;
                    sysModel.P2 = P2;
                    sysModel.P4 = P4;
                    sysModel.unit = 'nd';
                    switch frame
                        case 'B1centP1P2rot'
                            char.mu = muP1P2;
                            char.mstar = mstarP1P2;
                            char.lstar = lstarP1P2;
                            char.tstar = tstarP1P2;
                            char.muP4B1 = muP4B1;
                            char.mstarP4B1 = mstarP4B1;
                            char.lstarP4B1 = lstarP4B1;
                            char.tstarP4B1 = tstarP4B1;
                        case 'B2centP4B1rot'
                            char.mu = muP4B1;
                            char.mstar = mstarP4B1;
                            char.lstar = lstarP4B1;
                            char.tstar = tstarP4B1;
                            char.muP1P2 = muP1P2;
                            char.mstarP1P2 = mstarP1P2;
                            char.lstarP1P2 = lstarP1P2;
                            char.tstarP1P2 = tstarP1P2;
                    end
                    char.massP4 = mP4; %nondim mass by EM mstar
                    char.distP4 = aP4; %nondim dist by EM lstar
                    char.angVelP4 = wP4;
                    char.descNodeP4 = descNodeP4;
                    char.incP4 = incP4;
                    char.theta0P4 = theta0P4;
                    sysModel.char = char;
                case 'EPHEMERIS'
                    %temporarily use Beom's code
                    global MU_EARTH MU_MOON MU_SUN SMA_EARTHMOON SMA_SUNEARTH
                    setGlobalVariableIM;

                    EMFrame = setFrame('EARTH', 'MOON', MU_EARTH, MU_MOON, ...
                        SMA_EARTHMOON, 'J2000', 'EARTH');
                    
                    SEFrame = setFrame('SUN', 'EARTH', MU_SUN, MU_EARTH, ...
                        SMA_SUNEARTH, 'J2000', 'EARTH');
                    
                    mu = EMFrame.mu;
                    lstar = EMFrame.lstar;
                    tstar = EMFrame.tstar;
                    
                    EMFrame2 = EMFrame;
                    EMFrame2.centralBody = 'MOON';
                    
                    Body.ID = {'EARTH', 'MOON', 'SUN'};
                    Body.GM = [MU_EARTH, MU_MOON, MU_SUN];
                    
                    char.lstar = c_dim_quant(lstar,'km');
                    char.tstar = c_dim_quant(tstar,'sec');
                    char.mstar = c_dim_quant(1,'kg');
                    char.mu = mu;
                    char.EMFrame = EMFrame;
                    char.SEFrame = SEFrame;
                    char.EMFrame2 = EMFrame2;
                    char.Body = Body;
                    char.JDFix = JDFix;

                    sysModel.unit = 'nd Earth';
                    sysModel.char = char;
                case 'NBODY'
                    %need more here
            end
            sysModel.system_dynamics = systemDynamics;
            sysModel.frame = frame;
        end

        function pos = getPrimPos(obj,primaryName)
            %uses combo of system dynamics and frame info to get primary positions
            err = false;
            primaryName = lower(primaryName);
            switch obj.system_dynamics
                case '2BP'
                    pos = [0;0;0];
                case {'CR3BP','BCR4BP'}
                    switch obj.frame
                        case 'B1centP1P2rot'
                            switch primaryName
                                case obj.P1.name
                                    pos = [-obj.char.mu;0;0];
                                case obj.P2.name
                                    pos = [1-obj.char.mu;0;0];
                                otherwise
                                    error(['Body **' primaryName '** does not match the bodies associated with the system.'])
                            end
                        case 'B2centP4B1rot'
                            c = obj.char;
                            switch primaryname
                                case obj.P4.name
                                    pos = [-c.mu;0;0];
                                case obj.P1.name
                                    P1radius = c.muP1P2*c.lstarP1P2/c.lstar;
                                    rads = linspace(0,2*pi,100);
                                    x = cos(rads); y = sin(rads); z = zeros(1,length(rads));
                                    xP1 = P1radius*x; yP1 = P1radius*y;
                                    pos = [xP1;yP1;z];
                                    %rotate into OoP orbits if necessary
                                    rot1 = DCM(-incP4,'y',false);
                                    rot2 = DCM(2*pi-descNodeP4,'z',false);
                                    for i = 1:length(z)
                                       pos(:,i) = rot2*rot1*pos(:,i);
                                    end
                                    pos(1,:) = pos(1,:) + 1-c.mu;
                                case obj.P2.name
                                    P2radius = (1-c.muP1P2)*c.lstarP1P2/c.lstar;
                                    rads = linspace(0,2*pi,100);
                                    x = cos(rads); y = sin(rads); z = zeros(1,length(rads));
                                    xP2 = P2radius*x; yP2 = P2radius*y;
                                    pos = [xP2;yP2;z];
                                    %rotate into OoP orbits if necessary
                                    rot1 = DCM(-incP4,'y',false);
                                    rot2 = DCM(2*pi-descNodeP4,'z',false);
                                    for i = 1:length(z)
                                       pos(:,i) = rot2*rot1*pos(:,i);
                                    end
                                    pos(1,:) = pos(1,:) + 1-c.mu;
                                otherwise
                                    error(['Body **' primaryName '** does not match the bodies associated with the system.'])
                            end
                        case 'P1centinert'
                            switch primaryName
                                case obj.P1.name
                                    pos = [0;0;0];
                                otherwise
                                    error(['Body **' primaryName '** either does not match the bodies associated with the system, or cannot be plotted in this frame'])
                            end
                        case 'P2centinert'
                            switch primaryName
                                case obj.P2.name
                                    pos = [0;0;0];
                                otherwise
                                    error(['Body **' primaryName '** either does not match the bodies associated with the system, or cannot be plotted in this frame'])
                            end
                        otherwise
                            err = true;
                    end
                otherwise
                    err = true;
            end
            if err
                error('Current combo of system dynamics and frame not implemented.')
            end
            pos = c_dim_quant(pos,'nd_l');
        end

        function obj = makeND(obj,lstar,tstar,mstar)
            if ~strcmp(obj.system_dynamics,'2BP')
                error('not implemented yet for anything but 2BP')
            end
            ndchar = [lstar, tstar, mstar];
            obj.P1.radius = obj.P1.radius.change_unit('nd_l',ndchar);
            obj.P1.mass = obj.P1.mass.change_unit('nd_m',ndchar);
            obj.P1.mu = obj.P1.mu.change_unit('nd_mu',ndchar);
            obj.char.mu = obj.P1.mu;
            obj.char.lstar = c_dim_quant(lstar,'nd_l');
            obj.char.tstar = c_dim_quant(tstar,'nd_t');
            obj.char.mstar = c_dim_quant(mstar,'nd_m');
        end
    end
end
%%
function [mu3, mstar, lstar, tstar] = CR3BP_data(P1,P2,method)
switch method
    case 'atd'
        % Use ATD to get all important values
        atd.util.initJava;
        import edu.purdue.mbd.atd.crtbp.*;
        % -- Build system
        sys = SystemData(P1.name,P2.name);
        % -- Get quantities
        mu3 = sys.getMassRatio;
        mstar = sys.getCharacteristicMass;
        lstar = sys.getCharacteristicLength;
        tstar = sys.getCharacteristicTime;
    case 'manual'
        error('Not implemented yet')
end
mstar = c_dim_quant(mstar,'kg');
lstar = c_dim_quant(lstar,'km');
tstar = c_dim_quant(tstar,'sec');
mu3 = c_dim_quant(mu3,'nd_mu');
end
%%
function [muP1P2, mstarP1P2, lstarP1P2, tstarP1P2,mP4,aP4,wP4,mstarP4B1,muP4B1,lstarP4B1,tstarP4B1] = BCR4BP_data(P1,P2,P4,theta0P4,method)
switch method
    case 'atd'
        atd.util.initJava;
        import edu.purdue.mbd.atd.bcfbp.*;
        sys = SystemData(P1.name, P2.name, P4.name, theta0P4.value);
        % -- Get quantities
        muP1P2 = sys.getP1P2MassRatio(); %from ATD
        mstarP1P2 = sys.getCharacteristicMass();
        lstarP1P2 = sys.getCharacteristicLength(); %km from ATD
        tstarP1P2 = sys.getCharacteristicTime(); % from ATD
        aP4 = sys.getP3MeanDistanceRatio();
        mP4 = sys.getP3MassRatio();
        MassP4 = mP4*mstarP1P2;
        wP4 = sys.getP3AngularVelocity(); %should be -0.9251...
        mstarP4B1 = MassP4 + mstarP1P2;
        muP4B1 = 1 - MassP4/(mstarP4B1);
        lstarP4B1 = aP4*lstarP1P2;
        tstarP4B1 = sqrt((lstarP4B1*1000)^3/mstarP4B1/6.67430e-11);
    case 'manual'
        error('Not implemented yet')
end
muP1P2 = c_dim_quant(muP1P2,'nd_mu');
mstarP1P2 = c_dim_quant(mstarP1P2,'kg');
tstarP1P2 = c_dim_quant(tstarP1P2,'sec');
lstarP1P2 = c_dim_quant(lstarP1P2,'km');
muP4B1 = c_dim_quant(muP4B1,'nd_mu');
mstarP4B1 = c_dim_quant(mstarP4B1,'kg');
tstarP4B1 = c_dim_quant(tstarP4B1,'sec');
lstarP4B1 = c_dim_quant(lstarP4B1,'km');
mP4 = c_dim_quant(mP4,'nd_m');
aP4 = c_dim_quant(aP4,'nd_l');
wP4 = c_dim_quant(wP4,'nd_av');
end