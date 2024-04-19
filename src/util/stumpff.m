function y = stumpff(x,C_or_S)
%y = stumpff(x,C_or_S)
%   x:
%   C_or_S:
series = false;
if series
    y = 0;
    n = 15;
end
if strcmpi(C_or_S,'C')
    if series
        for i = 0:n
            y = y + (-x)^i/factorial(2*i+2);
        end
    else
        if x > 0
            y = (1-cos(sqrt(x))) / x;
        elseif x == 0
            y = 1/2;
        else
            y = (cosh(sqrt(-x))-1) / -x;
        end
    end
elseif strcmpi(C_or_S,'S')
    if series
        for i = 0:n
            y = y + (-x)^i/factorial(2*i+3);
        end
    else
        if x > 0
            y = (sqrt(x)-sin(sqrt(x))) / x^(3/2);
        elseif x == 0
            y = 1/6;
        else
            y = (sinh(sqrt(-x))-sqrt(-x)) / (-x)^(3/2);
        end
    end
elseif strcmpi(C_or_S,'Cdot')
    if series
        for i = 0:n
            y = y + ((-1)^i*i*x^(i-1))/factorial(2*i+2);
        end
    else
        y = 1/(2*x) * (1-x*stumpff(x,'S')-2*stumpff(x,'C'));
    end
elseif strcmpi(C_or_S,'Sdot')
    if series
        for i = 0:n
            y = y + ((-1)^i*i*x^(i-1))/factorial(2*i+3);
        end
    else
        y = 1/(2*x)*(stumpff(x,'C') - 3*stumpff(x,'S'));
    end 
else
    error('No valid C_or_S input. Must be ''C'', ''Cdot'', ''S'', or ''Sdot''.')
end
end