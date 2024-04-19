function [v1, v2] = gaussProblem(mu,r1,r2,dt,dmot,method)
%[v1, v2] = gaussProblem(mu,r1,r2,dt,dmot,method)
%   mu: center body gm
%   r1: initial position [x;y;z]
%   r2: final position [x;y;z]
%   dt: change in time (sec)
%   dmot: direction of motion ('short' or 'long') for hyperbolic vs elliptical
%   method: either 'univ' or 'p-iter'

r1mag = norm(r1); r2mag = norm(r2);
cosNu = dot(r1,r2)/(r1mag*r2mag); %cosine of change in true anomaly
nu = acos(cosNu);
while nu < 0
   nu = nu + 2*pi; 
end
if strcmpi(method,'univ')
    if strcmpi(dmot,'long')
        DM = -1;
    elseif strcmpi(dmot,'short')
        DM = 1;
    else
        error('dmot must be either ''long'' or ''short''.')
    end
    A = DM*sqrt(r1mag*r2mag*(1+cos(nu))); %BMW pg 237
    z_up = 4*pi^2;
    z_low = -4*pi;
    z = 0;
    i=0;
    while i < 7 %six iterations of Bisection algorithm
        S_i = stumpff(z,'S');
        C_i = stumpff(z,'C');
        y = y_func(z,A,r1mag,r2mag);
        if y < 0
            z = z + pi; %not sure what I should be changing z by
            continue
        end
        x = sqrt(y/C_i);
        t = 1/sqrt(mu)*(x^3*S_i+A*sqrt(y));
        if t < dt
            z_low = z; 
        else
            z_up = z;
        end
        z = z_low + (z_up-z_low)/2;
        i = i+1;
    end
    %then call Newton's method
    slow = 0.001;
    func = @(Z) slow*(dt - 1/sqrt(mu)*((x_func(Z,A,r1mag,r2mag))^3*stumpff(Z,'S')...
        + A * sqrt(y_func(Z,A,r1mag,r2mag)))); %BMW pg 234
    dfunc = @(Z) -1/sqrt(mu)*(x_func(Z,A,r1mag,r2mag)^3*(stumpff(Z,'Sdot')...
        - 3/2 * stumpff(Z,'S')*stumpff(Z,'Cdot')/stumpff(Z,'C'))...
        + A/8 * (2*stumpff(Z,'S')*sqrt(y_func(Z,A,r1mag,r2mag))/stumpff(Z,'C')...
        + A/x_func(Z,A,r1mag,r2mag)));
    z = newtonMethod(func,dfunc,z,10^-6,'diff');
    y = y_func(z,A,r1mag,r2mag);
    f = 1-y/r1mag;
    g = A*sqrt(y/mu);
    gdot = 1 - y/r2mag;
    v1 = (r2 - f*r1)/g;
    v2 = (gdot*r2 - r1)/g;
elseif strcmpi(method,'p-iter') %BMW pg 241-251
    piterMethod = 'newton';
    k = r1mag*r2mag*(1-cos(nu));
    l = r1mag + r2mag;
    m = r1mag*r2mag*(1+cos(nu));
    p_i = k/(l+sqrt(2*m));
    p_ii = k/(l-sqrt(2*m));
    if strcmp(piterMethod,'newton')
        switch dmot
            case 'short' %hyperbola
                p = p_i + (p_ii-p_i)/2;
            case 'long' %ellipse
                p = p_ii + (p_ii-p_i)/2;
        end
        a = a_func(m,k,l,p);
    %     g = g_func(r1mag,r2mag,nu,mu,p);
        if a > 0 % ellipse
    %         E = E_func(r1mag,r2mag,p,nu,mu,a);
    %         t = g + sqrt(a^3/mu)*(E - sin(E));
            slow = 1; %helps Newton iteration from running too fast
            func = @(p) slow*(dt - (g_func(r1mag,r2mag,nu,mu,p) + sqrt((a_func(m,k,l,p))^3/mu)*...
                (E_func(r1mag,r2mag,p,nu,mu,a_func(m,k,l,p))...
                - sin(E_func(r1mag,r2mag,p,nu,mu,a_func(m,k,l,p))))));
            dfunc = @(p) -(-g_func(r1mag,r2mag,nu,mu,p)/(2*p)...
                - 3/2*a_func(m,k,l,p)*(dt - g_func(r1mag,r2mag,nu,mu,p))*...
                ((k^2+(2*m-l^2)*p^2)/(m*k*p^2))...
                + sqrt((a_func(m,k,l,p))^3/mu)*...
                2*k*sin(E_func(r1mag,r2mag,p,nu,mu,a_func(m,k,l,p))) / (p*(k-l*p)));
            p = newtonMethod(func,dfunc,p,10^-6,'diff');
        elseif a < 0 % hyperbola
    %         F = F_func(r1mag,r2mag,p,nu,a);
    %         t = g + sqrt((-a)^3/mu)*(sinh(F) - F);
            func = @(p) dt - (g_func(r1mag,r2mag,nu,mu,p) + sqrt((-a_func(m,k,l,p))^3/mu)*...
                (sin(F_func(r1mag,r2mag,p,nu,a_func(m,k,l,p)))...
                - F_func(r1mag,r2mag,p,nu,a_func(m,k,l,p))));
            dfunc = @(p) -(-g_func(r1mag,r2mag,nu,mu,p)/(2*p)...
                - 3/2*a_func(m,k,l,p)*(dt - g_func(r1mag,r2mag,nu,mu,p))*...
                ((k^2+(2*m-l^2)*p^2)/(m*k*p^2))...
                - sqrt((-a_func(m,k,l,p))^3/mu)*...
                2*k*sinh(F_func(r1mag,r2mag,p,nu,a_func(m,k,l,p))) / (p*(k-l*p)));
            p = newtonMethod(func,dfunc,p,10^-6,'diff');
        end
    elseif strcmpi(piterMethod,'bisection')
        switch dmot
            case 'long' %ellipse
                p_up = p_ii;
                p_low = p_i;
                p = p_ii + (p_ii-p_i)/2;
            case 'short' %hyperbola
                p_low = p_ii;
                p_up = p_low*1000; %arbitrary?
                p = p_i + (p_ii-p_i)/2;
        end
        con = inf;
        while abs(con) > 10^-5
            a = a_func(m,k,l,p);
            g = g_func(r1mag,r2mag,nu,mu,p);
            if a > 0 %long way, ellipse
                E = E_func(r1mag,r2mag,p,nu,mu,a);
                t = g + sqrt(a^3/mu)*(E - sin(E));
                if t < dt
                    p_low = p;
                else
                    p_up = p;
                end
            else %short way, hyperbola
                F = F_func(r1mag,r2mag,p,nu,a);
                t = g + sqrt((-a)^3/mu)*(sinh(F) - F);
                if t > dt
                    p_low = p;
                else
                    p_up = p;
                end
            end
            pn = p_low + (p_up-p_low)/2;
            con = p-pn;
            p = pn;
        end
    end
    g = g_func(r1mag,r2mag,nu,mu,p);
    gdot = 1 - r1mag/p*(1-cos(nu));
    f = f_func(r2mag,p,nu);
    v1 = (r2 - f*r1)/g;
    v2 = (gdot*r2 - r1)/g;
end
end
function a = a_func(m,k,l,p)
a = m*k*p/((2*m-l^2)*p^2 + 2*k*l*p - k^2);
end
function F = F_func(r1mag,r2mag,p,nu,a)
coshF = 1 - r1mag/a*(1-f_func(r2mag,p,nu));
F = acosh(coshF);
end
function E = E_func(r1mag,r2mag,p,nu,mu,a)
cosE = 1 - r1mag/a*(1-f_func(r2mag,p,nu));
sinE = -r1mag*r2mag*fdot_func(mu,p,nu,r1mag,r2mag)/sqrt(mu*a);
E = acos(cosE);
if sinE < 0
    E = -E; 
end
while E < 0
    E = E + 2*pi;
end
end
function f = f_func(r2mag,p,nu)
f = 1 - r2mag/p * (1 - cos(nu));
end
function g = g_func(r1mag,r2mag,nu,mu,p)
g = r1mag*r2mag*sin(nu) / sqrt(mu*p);
end
function fdot = fdot_func(mu,p,nu,r1mag,r2mag)
fdot = sqrt(mu/p) * tan(nu/2) * ((1-cos(nu))/p - 1/r1mag - 1/r2mag);
end
function x = x_func(z,A,r1mag,r2mag)
x = sqrt(y_func(z,A,r1mag,r2mag)/stumpff(z,'C'));
end
function y = y_func(z,A,r1mag,r2mag)
y = r1mag + r2mag - A * (1 - z*stumpff(z,'S'))/sqrt(stumpff(z,'C'));
end