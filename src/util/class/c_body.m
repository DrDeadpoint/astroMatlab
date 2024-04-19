classdef c_body
    %obj = c_body(bodyString)
    
    properties
        name char
        radius c_dim_quant %km
        mass c_dim_quant %kg
        mu c_dim_quant
    end
    
    methods
        function obj = c_body(bodyString)
            %obj = c_body(bodyString)
            [radius,mass,mu] = int_getBodyData(bodyString);
            obj.name = lower(bodyString);
            obj.radius = c_dim_quant(radius,'km');
            obj.mass = c_dim_quant(mass,'kg');
            obj.mu = c_dim_quant(mu, 'km3/s2');
        end
        
        function H = plot(obj,ax,sysModel,varargin)
            %H = plot(obj,ax,sysModel)
            %Adds this body to the axes defined by ax
            plotObjs = findall(ax);
            alreadyPlotted = false;
            for i = 1:length(plotObjs) %check to see if the body is already there
               if isa(plotObjs(i),'matlab.graphics.chart.primitive.Surface') ...
                       || isa(plotObjs(i),'matlab.graphics.chart.primitive.Line')
                   if strcmp(plotObjs(i).DisplayName,obj.name)
                      alreadyPlotted = true; %don't plot body again
                      H = plotObjs(i);
                   end
               end
            end
            if ~alreadyPlotted
                rad = obj.radius;
                pos = sysModel.getPrimPos(obj.name);
                if nargin > 3 %change units to user defined
                    desUnit = varargin{1};
                    rad = rad.change_unit(desUnit,sysModel);
                    pos = pos.change_unit(desUnit,sysModel);
                else %default is nd
                    rad = rad.change_unit(sysModel.char.lstar.unit,sysModel);
                    pos = pos.change_unit(sysModel.char.lstar.unit,sysModel);
                end
                pos = pos.value;
                if size(pos,1) == 3 %body is fixed
                    H = plot3DBody(ax,obj.name, rad, pos);
                    set(H,'HandleVisibility','off','DisplayName',obj.name) 
                else
                    col = sysModel.getPrimCol(obj.name);
                    H = plotStates(ax,pos,'-',...
                        'color',col,...
                        'LineWidth',linewidth/2,...
                        'DisplayName',obj.name,...
                        'HandleVisibility','off');
                end
            end
        end
    end
end

function [radius,mass,mu] = int_getBodyData(bodyString)
switch lower(bodyString)
    case 'sun'
        radius = 695508;
        mass = 1.9885e+30;
    case 'mercury'
        radius = 2439.7;
        mass = 3.3011e+23;
    case 'venus'
        radius = 6051.8;
        mass = 4.8675e+24;
    case 'earth'
        radius = 6371;
        mass = 5.97236526e+24; %to get desired mu
    case 'moon'
        radius = 1737.4;
        mass = 7.34603131e+22;
    case 'mars'
        radius = 3389.5;
        mass = 6.41706422e+23;
    case 'phobos'
        radius = 11.2667;
        mass = 1.0659e+16;
    case 'jupiter'
        radius = 69911;
        mass = 1.8982e+27;
    case 'saturn'
        radius = 58232;
        mass = 5.6834e+26;
    case 'uranus'
        radius = 25362;
        mass = 8.6813e+25;
    case 'neptune'
        radius = 24622;
        mass = 1.0241e+26;
    case 'pluto'
        radius = 1187;
        mass = 1.3030e+22;
    otherwise
        error(['Input body, ' bodyString ', not known.'])
end
G = 6.67408 * 10^-20; %km^3/kg/s^2
mu = mass * G; %km^3/s^2
end