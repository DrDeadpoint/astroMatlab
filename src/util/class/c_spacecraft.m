classdef c_spacecraft
    %SPACECRAFT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Tmax c_dim_quant
        Isp c_dim_quant
        g0 c_dim_quant
        M0 c_dim_quant
        m0ND c_dim_quant
        TmaxND c_dim_quant
        IspND c_dim_quant
        g0ND c_dim_quant
        throttle = 1 %{mustBePositive,mustBeLessThan(throttle,1)} = 1
        mdot c_dim_quant
        mdotND c_dim_quant
    end
    
    methods
        function obj = c_spacecraft(Tmax,Isp,M0,sysModel,varargin)
            obj.Tmax = Tmax.change_unit('kN',sysModel); %kilonewtons, kg*km/s^2
            obj.Isp = Isp.change_unit('sec',sysModel); %seconds
            obj.M0 = M0.change_unit('kg',sysModel); %kilograms
            g0 = c_dim_quant(9.80665/1000,'km/s2'); %km/s^2 (don't change, used two lines down)
            obj.g0 = g0;
            mdot = obj.Tmax.value/obj.Isp.value/obj.g0.value;
            obj.mdot = c_dim_quant(mdot,'kg/s');
            
%             sysModel.char.mstar = obj.M0; %briefly change how mstar is defined, kg
            obj.m0ND = obj.M0.change_unit('nd_m',sysModel,obj); %M0/M0; %kg/kg
            obj.TmaxND = obj.Tmax.change_unit('nd_f',sysModel,obj); %Tmax/M0*(tstar^2/lstar); %kN * sec^2/kgkm
            obj.IspND = obj.Isp.change_unit('nd_t',sysModel,obj); %Isp/tstar; %sec/sec
            obj.g0ND = obj.g0.change_unit('nd_a',sysModel,obj); %g0/(lstar/tstar^2); %km/sec^2 / km/sec^2
            mdotND = obj.TmaxND.value/obj.IspND.value/obj.g0ND.value;
            obj.mdotND = c_dim_quant(mdotND,'nd_md');
            
            if ~isempty(varargin)
                obj.throttle = varargin{1};
            end
        end
    end
end

