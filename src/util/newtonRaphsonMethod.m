function soln = newtonRaphsonMethod(problemSetup)
% soln = newtonRaphsonMethod(problemSetup)
% 
% 'problemSetup' is created using the function: buildProblem
% 'soln' is a structure containing:
%   x: feasible x found by newton raphson method
%   xvec: list of all intermediate x's at each iteration
%   normF: vector showing what the normF value was at each iteration
%   exitflag: tells how well the solution did
%       1: good convergence
%       0: hit max iterations
%       -1: error went above max allowable error
%       -2: stayed above max unconverged error for too long
%       -3: violated some facet of the problem, produced NaNs in the F vector
%       -4: some other error occurred, check the exitmessage
%   iter: how many iterations to reach solution
%   exitmessage: details on the exitflag
%
% Alex Hoffman
% 11/08/2021
options = problemSetup.options.solverOptions; %solver options only
param = problemSetup.param;
maxIter = options.maxIter;
fcnTol = options.fcnTol;
verbose = options.verbose;
attenuationFactor = options.attenuationFactor;
normType = options.normType;
invMethod = options.invMethod;
W = options.weights;
if isempty(W)
%     W = eye(length(param.x0));
    W = ones(size(param.x0));
end
conVec = {};
for i = 1:length(param.constraints)
    conVec = [conVec; param.constraints{i}.name];
end
if length(attenuationFactor) == 2 %stop attentuating after some number of iterations
    attenuationLength = attenuationFactor(2);
    attenuationFactor = attenuationFactor(1);
else
    attenuationLength = inf;
end
maxError = options.maxError; %how high can normF go?
notConvError = options.notConvError; %at what error level do we say it isn't converging after three iterations
numNotConv = 0;
x0 = param.x0; %initial guess
if ~iscolumn(x0)
    error('x0 must be a column vector of states')
end
objFunc = problemSetup.objFunc; %objective function
iter = 0;
exitflag = 1;
exitmessage = 'Converged well!';
if verbose
    fprintf('=== Solving problem...\n')
end
[F] = objFunc(x0); %for printing quickly
normF = norm(F,normType);
normFvec = normF;
if verbose
    fprintf('Newton Raphson Method, input. Error norm: %e\n', normF)
end
xvec = x0;
Fvec = F;
% try
if normF > fcnTol
    while 1
        iter = iter + 1;
        % check if max iteration has been reached
        if iter > maxIter
           iter = iter - 1;
           exitflag = 0;
           exitmessage = ['Max number of iterations reached, maxIter = ' num2str(maxIter)];
           break
        end
        % check for print
        if verbose
            fprintf('NRM, step %d. ', iter);
        end
        % update x
        if attenuationLength > iter
            aF = attenuationFactor;
        else
            aF = 1;
        end
%         dxs2try = [];
        [F,DF] = objFunc(x0); %now get DF for next iteration
        if size(DF,1) == size(DF,2) %solve
            % am i in a local basin, and pinv is just too good at taking
            % small steps
            % is pinv ever used in software for other problems?
            % use condition number to inform choice of method
            % check which constraints/free variables are not stepping, and
            % why the pseudoinverse is keeping that from happening (basin?)
            switch invMethod
                case 'qr'
                    [Q,R] = qr(DF);
                    dx = R\(Q\F);
                case 'svd'
                    [U,S,V] = svd(DF);
                    dx = V * inv(S) * U * F;
                case 'ppinv'
                    dx = DF' * pinv(DF*DF')*F;
                case 'inv'
                    dx = inv(DF) * F;
                case 'pinv'
                    dx = pinv(DF) * F;
                otherwise
                    error('unknown inverse method')
            end
% 
%             VarNames = {'inv', 'pinv', 'ppinv', 'qr', 'svd'};
% 
%             dxtable = table(dxinv, dxpinv, dxppinv, dxqr, dxsvd, 'VariableNames',VarNames);
%             dxs2try = [dxinv,dxpinv,dxppinv,dxqr,dxsvd];

        elseif size(DF,1) < size(DF,2)
%             [Q,R]  = qr(DF');
%             dxqr = Q*(inv(R')*F);
            switch invMethod
                case 'ppinv'
                    dx = DF' * pinv(DF*DF')*F;
                case 'qr'
                    [Q,R] = qr(DF);
                    dx = R\(Q\F);
                case 'wppinv'
%                     invW = inv(W);
%                     dx = invW * DF' * pinv(DF*invW*DF')*F;
                    dx = DF' * pinv(DF*DF')*F;
                    dx = W.*dx;
            end
        %     dx = aF*DF'*((DF*DF')\F);
%             dx = dxqr;
        else
            switch invMethod
                case 'qr'
                    [Q,R]  = qr(DF'*DF);
                    dx = R\(Q\(DF'*F));
                case 'ppinv'
                    dx = DF' * pinv(DF*DF')*F;
                case 'best'
                    [Q,R]  = qr(DF'*DF);
                    dxqr = R\(Q\(DF'*F));
                    dxppinv = DF' * pinv(DF*DF')*F;
                    dxs2try = [dxqr, dxppinv];
                otherwise
                    error('unknown inverse method')
            end
        end
        if strcmp(invMethod,'best')
            warning('this tends to get stuck')
            dxs2try = aF * dxs2try;
            Fii = inf(1,size(dxs2try,2));
            for ii = 1:size(dxs2try,2)
                if ~any(isnan(dxs2try))
                    thisx = x0 - dxs2try(:,ii);
                    Fvecii = objFunc(thisx);
                    Fii(ii) = norm(Fvecii);
                end
            end
            [~,minFind] = min(Fii);
            dx = dxs2try(:,minFind);
        end
        if any(isnan(dx))
            exitflag = -3;
            exitmessage = 'Something went wrong when evaluating the update; produced NaNs in the dx vector';
            break
        end
        dx = aF * dx;
        x = x0 - dx;
        % updates for next iteration
        x0 = x;
        xvec = [xvec, x];
        % evaluate objective function, check normF
        [F] = objFunc(x0); %first only get F to check for convergence
        if any(isnan(F))
            exitflag = -3;
            exitmessage = 'Something went wrong when evaluating the obj function; produced NaNs in the F vector';
            break
        end
        normF = norm(F,normType);
        normFvec = [normFvec, normF];
        Fvec = [Fvec, F];
        if verbose
            fprintf('Current error norm: %e\n', normF);
        end
        % check for convergence
        if normF < fcnTol % met the tolerance
           break 
        elseif normF > maxError % bigger than allowable tolerance
           exitflag = -1;
           exitmessage = ['Error grew too large: reached ' num2str(maxError)];
           break 
        elseif normF > notConvError % bigger than impatient tolerance
           if numNotConv == 2 %there were already two iterations directly previously which had this error
               exitflag = -2;
               exitmessage = ['Error stayed above ' num2str(notConvError) ' for ' num2str(numNotConv+1) ' iterations'];
               break
           else
               numNotConv = numNotConv + 1;
           end
        else
            numNotConv = 0; %set this equal to zero and move on
        end
    end
end
% catch ME
%     exitflag = -4;
%     exitmessage = ME.message;
%     warning(ME.message);
% end
soln.x = x0;
indx = 1:length(x0);
soln.xtable = table(indx',param.x0key,param.x0,x0); 
soln.xvec = xvec;
soln.normF = normFvec;
soln.Fvec = table(conVec,Fvec);
soln.exitflag = exitflag;
soln.iter = iter;
soln.exitmessage = exitmessage;
if verbose
    fprintf('%s\n',exitmessage)
    fprintf('Done.\n')
end
end