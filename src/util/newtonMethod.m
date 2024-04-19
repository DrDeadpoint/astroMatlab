function [x,iter] = newtonMethod(func,dfunc,x0,convFactor,convMethod)
% x = newtonMethod(func,dfunc,x0,convFact,convMethod)
%   func is the function you are trying to find the root of
%   dfunc is the derivative of func
%   x0 is the initial guess
%   convFactor is to stop the loop when convergence is achieved (10^-6)
%   convMethod is how to test for convergence, currently either:
%       "diff" checks difference between x_n+1 and x_n
%       "eval" evaluates func(x_n+1)

con = inf; %start large, then replace
iter = 0;
while abs(con)>convFactor
    F = func(x0);
    dF = dfunc(x0);
    x = x0 - F/dF;
    switch convMethod
        case 'diff'
            con = F/dF;
        case 'eval'
            con = func(x);
    end
    x0 = x;
    iter = iter+1;
end
end