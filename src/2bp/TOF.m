function [out] = TOF(r0,v0,mu,type_in,val)
%[t,nu] = TOF(r0,v0,mu,type,val)
%   TOF gives the time (sec) and change in angle (rad) since periapsis
%   periapsis is assumed to be t=0 and angle=0;
%   r0 and v0: are vectors in ECI coordinate system of satellite position at
%   time or angle since periapsis
%   mu: is GM of center body (km^3/s^2)
%   type: is either 'time' or 'angle' depending on what you are inputting
%   val: is either the change in time or change in angle you know, TOF will
%   output the corresponding change
%%
%this code does not use universal variable. Instead it calculates e and
%uses the conic section cases to find values.
%equations from Dr. Saleh's 4532 notes on Time
%%
%get useful values
h = norm(cross(r0,v0)); %angular momentum
energy = (norm(v0))^2/2-mu/norm(r0); %specific energy
p = ((norm(h))^2)/mu; %semi-latus rectum
a = -mu/(2*energy); %semi-major axis
ecc = sqrt(1-p/a); %eccentricity

if strcmp(type_in,'time') %given time, find angle
    t = val;
    if ecc==0 %circle
        out = t*(mu^2)/(h^3);
    elseif ecc<1 %ellipse
        Tellipse = 2*pi*(h^3)/((mu^2)*sqrt((1-ecc^2)^3)); %orbital period
        Me = 2*pi*t/Tellipse; %mean anomaly
        E = eccAnomaly(ecc,Me,'newton'); %'laguerre'
        out = 2*atan(sqrt((1+ecc)/(1-ecc))*tan(E/2));
        %challenge: write Newton-Raphson for equation for ellipse given
        %below, but solving for nu
    elseif ecc==1 %parabola
        out = 2*atan((3*(mu^2)/(h^3)*t+sqrt(((3*(mu^2)/(h^3)*t)^2+1)))^(1/3)...
            -(3*(mu^2)/(h^3)*t+sqrt(((3*(mu^2)/(h^3)*t)^2+1)))^(1/3));
    else %hyperbola
        Mh = t*(mu^2)*sqrt((ecc^2-1)^3)/(h^3); %mean anomaly
%         F1 = Mh; %eccentric anomaly
%         con = 1; %converges when con < 10^-6
%         while abs(con)>10^-6
%             con = (e*sinh(F1)-(F1+Mh))/(e*cosh(F1)-1);
%             F2 = F1 - con;
%             F1=F2;
%         end
        convFactor = 10^-6; %converges when con < 10^-6
        func = @(F) ecc*sinh(F);
        dfunc = @(F) (F+Mh)/(ecc*cosh(F)-1);
        F = newtonMethod(func,dfunc,Mh,convFactor,'diff');
        out = 2*atan(sqrt((ecc+1)/(ecc-1))*tanh(F/2));
        %challenge: write Newton-Raphson for equation given below for
        %hyperbolics, solving for nu
    end
    if out < 0
    out = out+2*pi;
    end
elseif strcmp(type_in,'angle') %given angle, find time
    nu = val;
    if ecc==0 %circle
        out = (h^3)/(mu^2)*nu;
    elseif ecc<1 %ellipse
        out = (h^3)/((mu^2)*sqrt((1-ecc^2)^3))*(2*atan(sqrt((1-ecc)/(1+ecc))*tan(nu/2))...
            -ecc*sqrt(1-ecc^2)*sin(nu)/(1+ecc*cos(nu)));
    elseif ecc==1 %parabola
        out = (h^3)/(mu^2)*(0.5*tan(nu/2)+(tan(nu/2)^3)/6);
    else %hyperbola
        out = (h^3)/((mu^2)*sqrt((ecc^2-1)^3))*(sqrt(ecc^2-1)*ecc*sin(nu)/(1+ecc*cos(nu))...
            -log((sqrt(ecc+1)+sqrt(ecc-1)*tan(nu/2))/(sqrt(ecc+1)-sqrt(ecc-1)*tan(nu/2))));
    end
end
end