function [x,iter] = laguerreMethod(n,func,dfunc,ddfunc,x0,convFactor,convMethod)
% x = laguerreMethod(n,func,dfunc,ddfunc,x0,convFactor,convMethod)
% Laguerre is usually used for determining roots of polynomials of degree n
%
% n is the degree of the polynomial you are trying to solve
% func is the function you are trying to find the root of
% dfunc is the derivative of func
% ddfunc is the derivative of dfunc
% x0 is the initial guess
% convFactor is to stop the loop when convergence is achieved (10^-6)
% convMethod is how to test for convergence, currently either:
%   "diff" checks difference between x_n+1 and x_n
%   "eval" evaluates func(x_n+1)

con = inf; %start large, then replace
iter = 0;
while abs(con)>convFactor
    F = func(x0);
    dF = dfunc(x0);
    ddF = ddfunc(x0);
    if dF > 0
        x = x0 - n*F/(dF + abs(sqrt(((n-1)*dF)^2-n*(n-1)*F*ddF)));
    else
        x = x0 - n*F/(dF - abs(sqrt(((n-1)*dF)^2-n*(n-1)*F*ddF)));
    end
    switch convMethod
        case 'diff'
            con = x-x0;
        case 'eval'
            con = func(x);
    end
    x0 = x;
    iter = iter + 1;
end
end