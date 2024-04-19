function [fs, varargout] = plotProblemIterations(soln,problemSetup,varargin)

irange = 1:soln.iter+1;
param = problemSetup.param;
fcnOptions = problemSetup.options.fcnOptions;

massplot = false;
finalIter = false;
initIter = false;
axin = {};
useDim = true;
useColor = false;
SB1plot = false;
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'massplot'
                massplot = true;
            case 'finalIter'
                finalIter = true;
            case 'noDim'
                useDim = false;
            case 'initIter'
                initIter = true;
            case 'useColor'
                useColor = true;
            case 'SB1plot'
                SB1plot = true;
        end
    else
        axin = varargin{i};
    end
end
if finalIter
    irange = irange(end);
end
if initIter
    irange = irange(1);
end

% plotting figure
if isempty(axin)
    [f,ax] = ffigure;
    hold(ax,'on');
    axis(ax,'equal');
else
    f = axin{1};
    ax = axin{2};
end

if massplot
    if isempty(axin)
        % mass figure
        [fmass, axmass] = ffigure;
        hold(axmass,'on');
    else
        f = axin{3};
        ax = axin{4};
    end
end
if SB1plot
    [fSB1, axSB1] = ffigure;
end

for i = irange
    x = soln.xvec(:,i);
    [segs, nodes] = getSegsNodes(x,param,fcnOptions,'noSTM',false);
    segs_o = segs; nodes_o = nodes;
    % plot on main figure
    [segs, nodes] = plotSegsNodes(segs,nodes,ax,i,irange(end),useDim,useColor);
    % change frame if desired
    if SB1plot
        for j = 1:length(segs_o)
            segs_o{j} = segs_o{j}.changeFrame('B2centP4B1rot');
            segs_o{j}.plotting.LTplotmult = 0.001;
        end
        for j = 1:length(nodes_o)
            nodes_o{j} = nodes_o{j}.changeFrame('B2centP4B1rot');
        end
        plotSegsNodes(segs_o,nodes_o,axSB1,i,irange(end),useDim,useColor);

    end
    % energy plot

    % mass plot
    if massplot
        for j = 1:length(segs)
            seg = segs{j};
            seg.time = seg.time.change_unit('day',seg.system_model);
            seg.low_thrust.mass = seg.low_thrust.mass.change_unit('kg',seg.low_thrust.spacecraft);
            if i == irange(end)
                plot(axmass,seg.time.value,seg.low_thrust.mass.value,...
                'Color',seg.plotting.lineColor,'LineWidth',seg.plotting.lineWidth,'DisplayName',seg.name);
            else
                plot(axmass,seg.time.value,seg.low_thrust.mass.value,...
                    'Color',seg.plotting.lineColor,'LineWidth',seg.plotting.lineWidth,'HandleVisibility','off');
            end
        end
        for j = 1:length(nodes)
            node = nodes{j};
            node.time = node.time.change_unit('day',node.system_model);
            node.low_thrust.mass = node.low_thrust.mass.change_unit('kg',node.low_thrust.spacecraft);
            if i == irange(end)
                plot(axmass,node.time.value, node.low_thrust.mass.value,node.plotting.markerStyle,...
                    'Color',node.plotting.lineColor,'DisplayName',node.name);
            else
                plot(axmass,node.time.value, node.low_thrust.mass.value,node.plotting.markerStyle,...
                    'Color',node.plotting.lineColor,'HandleVisibility','off');
            end
        end
        xlabel(axmass,'Time [days]')
        ylabel(axmass,'Mass [kg]')
        legend(axmass)
    end
end

varargout{1} = ax;
fs = {f};
if massplot
    fs{end+1} = fmass;
    varargout{2} = axmass;
end
if SB1plot
    fs{end+1} = fSB1;
    varargout{3} = axSB1;
end

end

function [segs, nodes] = plotSegsNodes(segs,nodes,ax,iter,maxiter,useDim,useColor)
    lstar = segs{1}.system_model.char.lstar;
    styleList = {':','--','-.'};
    for i = 1:length(segs)
        seg = segs{i};
        if iter < maxiter
            if useColor
                seg.plotting.lineColor = colour(iter);
                if iter == 1
                    seg.name = ['Initial guess'];
                else
                    seg.name = ['Iteration ' num2str(iter-1)];
                end
                seg.plotting.lineStyle = styleList{mod(iter,length(styleList))+1};
            else
                seg.plotting.lineColor = colour('gray');
                seg.plotting.HandleVisibility = 'off';
            end
            seg.plotting.lineWidth = 2;
        else
            if useColor
                seg.name = ['Final iteration'];
                seg.plotting.lineColor = colour(iter);
            else
                %do nothing
            end
        end
        if strcmp(seg.system_model.system_dynamics,'EPHEMERIS')
            plotStates(ax,seg.getAllStates.*lstar,'color',...
                seg.plotting.lineColor,'lineWidth',seg.plotting.lineWidth);
            plotStates(ax,seg.getStateByIndex(-1).*lstar,'.','color',...
                seg.plotting.lineColor);
            if ~isnan(seg.etc.seg.maneuver.frame) % a maneuver takes place
                u_TTLc = seg.etc.ephemOptions{5};
                if iter < maxiter
                    arrowColor = colour('gray'); %low thrust arcs get blue arrows
                    mag = 0.1*lstar;
                else
                    arrowColor = colour('red');
                    mag = 0.1*lstar;
                end
                nPlot = 20;
                T = seg.time.value;
                Y = seg.pos.value';
                nQ = floor(length(T)/nPlot);
                uIMat = u_TTLc{end}{1};
                for ii = 1:length(T)
                    if mod(ii, nQ) == 0
                        quiver3(ax,Y(ii,1)*lstar, Y(ii,2)*lstar, Y(ii,3)*lstar, ...
                            uIMat(ii,1)*mag, uIMat(ii,2)*mag, uIMat(ii,3)*mag, 'Color',arrowColor);
                    end
                end
            end
        else
            seg.plot('ax',ax,'useDim',useDim);
            if useDim
                plotStates(ax,seg.getStateByIndex(-1).*lstar,'.','color',...
                    seg.plotting.lineColor,'HandleVisibility','off');
            else
                plotStates(ax,seg.getStateByIndex(-1),'.','color',...
                    seg.plotting.lineColor,'HandleVisibility','off');
            end
        end
        segs{i} = seg;
    end
    for i = 1:length(nodes)
        node = nodes{i};
        if iter == maxiter
            if i == 1
                node.plotting.lineColor = colour('green');
                node.plotting.markerStyle = 'o';
            elseif i < length(nodes)
                node.plotting.lineColor = colour('blue');
                node.plotting.markerStyle = 's';
            else
                node.plotting.lineColor = colour('red');
                node.plotting.markerStyle = 'o';
            end
        else
            if i == 1
                node.plotting.lineColor = colour('green');
                node.plotting.markerStyle = '.';
            elseif i < length(nodes)
                node.plotting.lineColor = colour('blue');
                node.plotting.markerStyle = '.';
            else
                node.plotting.lineColor = colour('red');
                node.plotting.markerStyle = '.';
            end
            node.plotting.HandleVisibility = 'off';
        end
%         if iter < maxiter
%             node.plotting.markerSize = node.plotting.markerSize/2;
%         end
        if strcmp(node.system_model.system_dynamics,'EPHEMERIS')
            plotStates(ax,node.getAllStates.*lstar,node.plotting.markerStyle,'color',...
                node.plotting.lineColor,'markerSize',node.plotting.markerSize);
        else
            node.plot('ax',ax,'useDim',useDim);
        end
        nodes{i} = node;
    end
end