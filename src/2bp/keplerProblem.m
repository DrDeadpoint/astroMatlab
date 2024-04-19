function traj_out = keplerProblem(traj_in)
%[r,v] = keplerProblem(r0,v0,t,t0,mu)
%   r0, v0 are vectors
%   t is the time of interest, calculated with TOF if given nu
%   t0 is likely zero
%   mu is GM of center body

traj_out = traj_in;
sysModel = traj_in.system_model;
tspan = traj_in.getTimeSpan().value;
r0vec = traj_in.getInitPos().value;
v0vec = traj_in.getInitVel().value;
mu = sysModel.char.mu.value;

dt = tspan(2) - tspan(1);
r0 = norm(r0vec);
v0 = norm(v0vec);

method = 'shepperd';
switch method
    case 'stumpff'
        a_inv = 2/r0-v0^2/mu;
        
        % calculate universal anomaly
        if a_inv > 0
            X0 = sqrt(mu) * dt * abs(a_inv);
        else
            X0 = sign(dt)*sqrt(-1/a_inv)*log((-2*mu*dt)/(1/a_inv*(dot(r0vec,v0vec)...
                + sign(dt)*sqrt(-mu*1/a_inv)*(1-r0*a_inv))));
        end
        func = @(X) dt - 1/sqrt(mu)*(dot(r0vec,v0vec)/sqrt(mu)*X^2*stumpff(a_inv*X^2,'C') + ...
            (1-a_inv*r0)*X^3*stumpff(a_inv*X^2,'S') + r0*X); %BMW pg 198
        dfunc = @(X) -1/sqrt(mu)*(X^2*stumpff(a_inv*X^2,'C') ... %negative for newton method?
            + dot(r0vec,v0vec)/sqrt(mu)*X*(1-a_inv*X^2*stumpff(a_inv*X^2,'S')) ...
            + r0*(1-a_inv*X^2*stumpff(a_inv*X^2,'C')));
        % ddfunc: need derivative of dfunc; look at Saleh's notes on Time for help
        X = newtonMethod(func,dfunc,X0,10^-12,'diff');
        
        % calculate f and g functions
        f = 1 - X^2/r0*stumpff(a_inv*X^2,'C');
        g = dt - X^3/sqrt(mu)*stumpff(a_inv*X^2,'S');
        rvec = f*r0vec + g*v0vec;
        rmag = norm(rvec);
        fdot = X*sqrt(mu)/(rmag*r0)*(a_inv*X^2*stumpff(a_inv*X^2,'S')-1);
        gdot = 1 - X^2/rmag*stumpff(a_inv*X^2,'C');
        vvec = fdot*r0vec + gdot*v0vec;
        stm = c_stm(nan(8,8),{'x','y','z','xd','yd','zd', 'mass', 't0'});
    case 'shepperd'
        T = dt;
        nu0 = dot(r0vec,v0vec);
        beta = 2*mu/r0 - dot(v0,v0);
        u = 0; %initial guess
        dU = 0;
        if beta > 0 % elliptical orbits only
            P = 2*pi*mu*beta^(-3/2);
            n = floor(1/P * (T + P/2 - 2*nu0/beta));
            dU = 2*n*pi*beta^(-5/2);
        end
        % Kepler iteration loop
        t = inf;
        while abs(t-T) > 10^-10
            q = beta * u^2 / (1+beta*u^2);
            U0w2 = 1-2*q;
            U1w2 = 2*(1-q)*u;
            U = 16/15 * U1w2^5 * GaussFraction(-9,3,15,0,1,1,1,q) + dU;
            U0 = 2*U0w2^2 - 1;
            U1 = 2*U0w2*U1w2;
            U2 = 2*U1w2^2;
            U3 = beta*U + 1/3*U1*U2;
            r = r0*U0 + nu0*U1 + mu*U2;
            t = r0*U1 + nu0*U2 + mu*U3;
            u = u - (t-T)/(4*(1-q)*r);
        end
        % solution
        f = 1-(mu/r0)*U2;
        g = r0*U1 + nu0*U2;
        F = -mu*U1/(r*r0);
        G = 1-(mu/r)*U2;
        rvec = f*r0vec + g*v0vec;
        vvec = F*r0vec + G*v0vec;
        % stm
        W = g*U2 + 3*mu*U;
        M11 = (U0/(r*r0) + 1/r0^2 + 1/r^2)*F - mu/r^3*mu/r0^3*W;
        M12 = F*U1/r + (G-1)/r^2;
        M13 = (G-1)*U1/r - mu/r^3*W;
        M21 = -F*U1/r0 - (f-1)/r0^2;
        M22 = -F*U2;
        M23 = -(G-1)*U2;
        M31 = (f-1)*U1/r0 - mu/r0^3*W;
        M32 = (f-1)*U2;
        M33 = g*U2 - W;
        I = eye(3);
        p11 = f*I + [rvec vvec] * [M21 M22; M31 M32] * [r0vec v0vec]';
        p12 = g*I + [rvec vvec] * [M22 M23; M32 M33] * [r0vec v0vec]';
        p21 = F*I - [rvec vvec] * [M11 M12; M21 M22] * [r0vec v0vec]';
        p22 = G*I - [rvec vvec] * [M12 M13; M22 M23] * [r0vec v0vec]';
        stm = [p11 p12;
               p21 p22];
        stm = [stm, zeros(6,2);
               zeros(2,6), [1 0; 0 0]];
        stm = c_stm(stm, {'x','y','z','xd','yd','zd', 'mass', 't0'});
    case 'der'
        a = 2/r0-v0^2/mu;
        sigma0 = dot(r0vec,v0vec)/sqrt(mu);
        F = @(x) r0*fU1(x,a) + sigma0*fU2(x,a) + fU3(x,a) - sqrt(mu)*dt;
        Fp = @(x) r0*fU0(x,a) + sigma0*fU1(x,a) + fU2(x,a);
        Fpp = @(x) sigma0*fU0(x,a) + (1-a*r0)*fU1(x,a);
        xn = a*sqrt(mu)*dt; %x0
        xnp = xn*2; %init
        iter = 0;
        while abs(xnp - xn) > 10^-10
            xnp = xn - 5*F(xn) / (Fp(xn) + Fp(xn)/abs(Fp(xn)) * sqrt(...
                16*Fp(xn)^2 - 20*F(xn)*Fpp(xn)));
            iter = iter + 1;
            if iter > 15
                error('took too many iterations')
            end
            xn = xnp;
        end
        x = xnp;
        % calculate values with newfound solution for x
        U0 = fU0(x,a);
        U1 = fU1(x,a);
        U2 = fU2(x,a);
        U3 = fU3(x,a);
        r = r0*U0 + sigma0*U1 + U2;
        sigma = sigma0*U0 + (1-a*r0)*U1;
%         y = fy(x,a);
%         C = fC(x,a);
%         S = fS(x,a);
        f = 1-U2/r0;
        g = r0*U1/sqrt(mu) + sigma0*U2/sqrt(mu);
        fd = -sqrt(mu)/r/r0*U1;
        gd = 1-U2/r;
        
        % get kepler solution
        rvec = f*r0vec + g*v0vec;
        vvec = fd*r0vec + gd*v0vec;

        % get stm
        I = eye(3);
        c11 = 1/a/r/r0^2*(3*U1*U3 + (a*r0-2)*U2^2) + U1^2/r + U2/r0;
        c12 = v0*U1*U2/r/sqrt(mu);
        c13 = v0/a/r/r0^2/sqrt(mu)*(r0*U1*U2 + 2*sigma0*U2^2 + 3*U2*U3 - 3*r*U3 + a*r0^2*U1*U2);
        c14 = v0^2*U2^2/r/mu;
        c21 = r0*U1*U2/r/sqrt(mu);
        c22 = v0/a/r/mu * 3*U1*U3 + (a*r0-2)*U2^2;
        c23 = r0*v0*U2^2/r/mu;
        c24 = v0^2/a/r/mu/sqrt(mu)*r0*U1*U2 + 2*sigma0*U2^2 + 3*U2*U3 - 3*r*U3;
        c31 = sqrt(mu)/a/r^3/r0^2 * (r*(3*U0*U3 - U1*U2 + 2*a*r0*U1*U2) - sigma*(3*U1*U3 - 2*U2^2 + a*r0*U2^2)) +...
            sqrt(mu)/r^3/r0 * (2*r*r0*U0*U1 + r^2*U1 -  sigma*r0*U1^2);
        c32 = v0/r^3*(r*(U0*U2+U1^2) - sigma*U1*U2);
        c33 = -3*v0*U2/a/r/r0^2 + v0/a/r^2/r0^2 * (4*sigma0*U1*U2 + r0*(U0*U2+U1^2) + 3*(U1*U3+U2^2))...
            - sigma*v0*U2/a/r^3/r0^2*(r0*U1 + 2*sigma0*U2 + 3*U3) + v0/r^3*(r*(U0*U2 + U1^2) - sigma*U1*U2);
        c34 = v0^2/r^3/sqrt(mu)*(2*r*U1*U2 - sigma*U2^2);
        c41 = r0/r^3*(r*(U0*U2 + U1^2) - sigma*U1*U2);
        c42 = v0/a/r^3/sqrt(mu)*(r*(3*U0*U3 - U1*U2 + 2*a*r0*U1*U2) - sigma*(3*U1*U3 - 2*U2^2 + a*r0*U2^2));
        c43 = r0*v0/r^3/sqrt(mu)*(2*r*U1*U2 - sigma*U2^2);
        c44 = -3*v0^2*U2/a/r/mu + v0^2/a/r^2/mu * (4*sigma0*U1*U2 + r0*(U0*U2 + U1^2) + 3*(U1*U3 + U2^2)) ...
            - sigma*v0^2*U2/a/r^3/mu*(r0*U1 + 2*sigma0*U2 + 3*U3);
        M1 = (r0vec * r0vec')/r0^2;
        M2 = (r0vec * v0vec')/r0/v0;
        M3 = (v0vec * r0vec')/r0/v0;
        M4 = (v0vec * v0vec')/v0^2;
        p11 = f*I + c11*M1 + c12*M2 + c13*M3 + c14*M4;
        p12 = g*I + c21*M1 + c22*M2 + c23*M3 + c24*M4;
        p21 = fd*I + c31*M1 + c32*M2 + c33*M3 + c34*M4;
        p22 = gd*I + c41*M1 + c42*M2 + c43*M3 + c44*M4;
        stm = [p11 p12; 
               p21 p22];
        stm = [stm, zeros(6,2);
               zeros(2,6), [1 0; 0 0]];
        stm = c_stm(stm, {'x','y','z','xd','yd','zd', 'mass', 't0'});
end

traj_out.pos = c_dim_quant([r0vec rvec],'nd_l');
traj_out.vel = c_dim_quant([v0vec vvec],'nd_v');
traj_out.low_thrust.mass.value = [traj_out.low_thrust.mass.value(1) traj_out.low_thrust.mass.value(1)];
traj_out.stm = stm;
end
%%
function G = GaussFraction(k,l,d,n,A,B,G,x)
    if abs(B) > 10^-12
        k = -k;
        l = l+2;
        d = d+4*l;
        n = n+(1+k)*l;
        A = d/(d-n*A*x);
        B = (A-1)*B;
        G = G+B;
        G = GaussFraction(k,l,d,n,A,B,G,x);
    end
end
%% 
function U0 = fU0(x,a)
    U0 = 1 - a*x^2*fC(x,a);
end
function U1 = fU1(x,a)
    U1 = x*(1-a*x^2*fS(x,a));
end
function U2 = fU2(x,a)
    U2 = x^2*fC(x,a);
end
function U3 = fU3(x,a)
    U3 = x^3*fS(x,a);
end
function C = fC(x,a)
    y = fy(x,a);
    if y > 0
        C = (1-cos(sqrt(y)))/y;
    elseif y == 0
        C = 1/2;
    else
        C = (cosh(sqrt(-y))-1)/(-y);
    end
end
function S = fS(x,a)
    y = fy(x,a);
    if y > 0
        S = (sqrt(y) - sin(sqrt(y)))/(sqrt(y)^3);
    elseif y == 0
        S = 1/6;
    else
        S = (sinh(sqrt(-y)) - sqrt(-y))/(sqrt(-y)^3);
    end
end
function y = fy(x,a)
    y = a*x^2;
end