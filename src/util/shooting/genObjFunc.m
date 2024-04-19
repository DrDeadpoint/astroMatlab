function [F,varargout] = genObjFunc(x,param,fcnOptions)
% [F,DF] = genObjFunc(x,param,fcnOptions)
% 'x' is a column vector of free variables
% 'param' contains segments and constraints
% 'fcnOptions' is passed in from the bigger options struct
%
% Alex Hoffman
% 11/15/21
    
switch fcnOptions.DFmethod
    case 'numeric'
        [segs,nodes] = getSegsNodes(x,param,fcnOptions,'noSTM');
        if nargout == 1
            [F] = all_numeric_partials(x,param,fcnOptions,segs,nodes);
        else
            [F,DF] = all_numeric_partials(x,param,fcnOptions,segs,nodes);
            varargout{1} = DF;
        end
    case 'analytic'
        [segs,nodes] = getSegsNodes(x,param,fcnOptions,'yesSTM');
        if nargout == 1
            [F] = analytic_partials(x,param,fcnOptions,segs,nodes);
        else
            [F,DF] = analytic_partials(x,param,fcnOptions,segs,nodes);
            varargout{1} = DF;
        end
end
end
%%
function [F,varargout] = analytic_partials(x,param,fcnOptions,segs,nodes)
    checkAnalyticPartials = false; %for debugging only
    con = param.constraints;
    m = length(con);
    F = obj(param,segs,nodes);
    if nargout == 2
        DF = nan(m,length(x));
        for i = 1:m
            thisCon = con{i};
            DF(i,:) = thisCon.eval_df(segs,nodes,param.x0key);
        end
        nanLoc = isnan(DF);
        nanCol = any(nanLoc,1); %find columns with nan (unknown partials)
        if any(nanCol) % replace values where I don't have an analytic solution
            %calculate only the partials I don't know
            [DFnum] = some_numeric_partials(x,param,fcnOptions,nanCol);
            DF(nanLoc) = DFnum(nanLoc);
        end
        if checkAnalyticPartials 
            [~,DFnum] = all_numeric_partials(x,param,fcnOptions,segs,nodes);
%             Ferr = (F - Fnum) ./ Fnum;
%             fprintf('Max relerr in F: %.3e\n',max(Ferr));
            DFrelerr = (DF - DFnum) ./ DFnum;
            DFabserr = (DF - DFnum);
            [~,ir] = max(abs(DFrelerr),[],'all');
            [rr,cr] = ind2sub(size(DFrelerr),ir);
            [~,ia] = max(abs(DFabserr),[],'all');
            [ra,ca] = ind2sub(size(DFabserr),ia);
            fprintf('Max relerr in DF: %.3e\n',DFrelerr(rr,cr));
            fprintf('   Corresponds to: %s\n',con{rr}.name)
            fprintf('   For parameter: %s\n',param.x0key{cr})
            fprintf('   Actual error: %.3e\n',DFabserr(rr,cr))
            fprintf('Max abserr in DF: %.3e\n',DFabserr(ra,ca));
            fprintf('   Corresponds to: %s\n',con{ra}.name)
            fprintf('   For parameter: %s\n',param.x0key{ca})
        end
        varargout{1} = DF;
    end
end
%%
function [F,varargout] = all_numeric_partials(x,param,fcnOptions,segs,nodes)
    F = obj(param,segs,nodes);
    if nargout == 2
        DF = zeros(length(param.constraints),length(x));
        h = fcnOptions.step_centralDifferencing;
        for i = 1:length(x)
            hVec = zeros(1, length(x));
            hVec(i) = h;
            [pSTs,pNodes] = getSegsNodes(x+hVec', param, fcnOptions, 'noSTM'); %propogate using slight changes in x
            [nSTs,nNodes] = getSegsNodes(x-hVec', param, fcnOptions, 'noSTM');
            Fp = obj(param, pSTs, pNodes);
            Fm = obj(param, nSTs, nNodes);
            DF(:,i) = (Fp - Fm)/2/h;
        end
        varargout{1} = DF;
    end
end

function [DF] = some_numeric_partials(x,param,fcnOptions,nanCol)
    DF = zeros(length(param.constraints),length(x));
    h = fcnOptions.step_centralDifferencing;
    for i = find(nanCol)
        hVec = zeros(1, length(x));
        hVec(i) = h;
        [pSegs,pNodes] = getSegsNodes(x+hVec', param, fcnOptions, 'noSTM'); %propogate using slight changes in x
        [nSegs,nNodes] = getSegsNodes(x-hVec', param, fcnOptions, 'noSTM');
        Fp = obj(param, pSegs, pNodes);
        Fm = obj(param, nSegs, nNodes);
        DF(:,i) = (Fp - Fm)/2/h;
    end
end
%%
function F = obj(param,segs,nodes)
    con = param.constraints;
    m = length(con);
    F = nan(m,1); %mx1 constraint vector
    for i = 1:m
        thisCon = con{i};
        F(i) = thisCon.eval_f(segs,nodes);
    end
end