function out = eci_elem(varargin)
in = varargin{1};
if isa(in,'c_elem')
    sysModel = varargin{2};
    checkClass(sysModel,'c_system_model')
    elem = in;
    
    sma = elem.sma;
    ecc = elem.ecc;
    raan = elem.raan.change_unit('rad');
    inc = elem.inc.change_unit('rad');
    argP = elem.argP.change_unit('rad');
    trueAnom = elem.trueAnom.change_unit('rad');
    if strcmp(elem.mu.unit,'nd_mu')
        sma = sma.change_unit('nd_l',sysModel);
    elseif strcmp(elem.mu.unit,'km3/s2')
        sma = sma.change_unit('km',sysModel);
    end
    smav = sma.value;
    ecc = ecc.value;
    raan = raan.value;
    inc = inc.value;
    argP = argP.value;
    nu = trueAnom.value;
    mu = elem.mu.value;
    per = elem.period;

    p = smav*(1-ecc^2); %semi-latus rectum
    rP1centinert = p/(1+ecc*cos(nu)); %magnitude of radius
    rPF = [rP1centinert*cos(nu); rP1centinert*sin(nu); 0]; %perifocal radius
    vPF = sqrt(mu/p)*[-sin(nu); ecc+cos(nu); 0]; %perifocal velocity
    rP1centinert = (DCM(raan,'z',true))' * (DCM(inc,'x',true))' * (DCM(argP,'z',true))' * rPF; %ECI radius
    vP1centinert = (DCM(raan,'z',true))' * (DCM(inc,'x',true))' * (DCM(argP,'z',true))' * vPF; %ECI velocity
    rP1centinert = c_dim_quant(rP1centinert,sma.unit);
    switch sma.unit
        case 'km'
            velUnit = 'km/s';
        case 'nd_l'
            velUnit = 'nd_v';
    end
    vP1centinert = c_dim_quant(vP1centinert,velUnit);

%     if ~strcmp(sysModel.system_dynamics, '2BP')
%         warning('Your system model is not a 2-body system model')
%     end
    sysModel.frame = 'P1centinert';
    tspan = c_dim_quant([0 per.value], per.unit);
    out = c_traj('2B orbit', tspan, rP1centinert, vP1centinert, sysModel);
elseif isa(in,'c_traj')
    traj = in;
    frame = traj.system_model.frame;
    switch frame(1:2)
        case 'P1'
            mu = traj.system_model.P1.mu;
        case 'P2'
            mu = traj.system_model.P2.mu;
        case 'P4'
            mu = traj.system_model.P4.mu;
        case {'B1','B2'}
            error('Current frame of reference uses a barycenter as the origin. Rotate into a primary centered frame.')
        case 'J2'
            seg = traj.etc.seg;
            GMind = strcmp(traj.system_model.char.Body.ID, seg.prop.centralBody);
            mu = traj.system_model.char.Body.GM(GMind);
            mu = c_dim_quant(mu,'km3/s2');
        otherwise
            error('This type of frame not yet implemented.')

    end
    switch traj.pos.unit
        case 'km'
            mu = mu.change_unit('km3/s2',traj.system_model);
            traj.vel = traj.vel.change_unit('km/s',traj.system_model);
        case 'nd_l'
            mu = mu.change_unit('nd_mu',traj.system_model);
            traj.vel = traj.vel.change_unit('nd_v',traj.system_model);
    end
    muv = mu.value;
    rECI = traj.pos.value(:,1);
    r = norm(rECI);
    vECI = traj.vel.value(:,1);
    v = norm(vECI);
    x = [1; 0; 0];
    z = [0; 0; 1];
    hECI = cross(rECI,vECI);
    h = norm(hECI);
    nodeECI = cross(z,hECI)/h;
    n = norm(nodeECI);
    eECI = (1/muv)*((v^2-muv/r)*rECI-dot(rECI,vECI)*vECI);
    ecc = norm(eECI);
    p = (h^2)/muv;
    sma = p/(1-ecc^2);
    inc = acos(dot(z,hECI)/h);
    RAAN = acos(dot(x,nodeECI)/n);
    if nodeECI(2) < 0
        RAAN = 2*pi - RAAN;
    end
    if isnan(RAAN) %planar motion
        RAAN = 0; %arbitrary
        nodeECI = x;
        n = norm(nodeECI);
    end
    argP = acos(dot(nodeECI,eECI)/(n*ecc));
    if dot(eECI,z) < 0
        argP = 2*pi - argP;
    end
    nu = acos(dot(eECI,rECI)/(ecc*r));
    if dot(rECI,vECI) < 0
        nu = 2*pi - nu;
    end
    if isnan(argP)
        argP = 0;
    end
    if isnan(nu)
       nu = 0; 
    end
    sma = c_dim_quant(sma,traj.pos.unit);
    ecc = c_dim_quant(ecc,'');
    inc = c_dim_quant(inc,'rad');
    RAAN = c_dim_quant(RAAN,'rad');
    argP = c_dim_quant(argP,'rad');
    nu = c_dim_quant(nu,'rad');
    out = c_elem(mu,sma,ecc,inc,RAAN,argP,nu);
else
    error(['First input must be either c_elem or c_traj. Yours was ' class(in)])
end
end