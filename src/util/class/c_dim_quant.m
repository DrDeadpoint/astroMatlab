classdef c_dim_quant
    %C_DIM_QUANT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        value {mustBeNumeric}
        unit char
    end
    
    methods
        function obj = c_dim_quant(value,unit)
            %C_DIM_QUANT Construct an instance of this class
            %   Detailed explanation goes here
            obj.value = value;
            obj.unit = unit;
        end

        function obj_out = index(obj,varargin)
            obj_out = obj;

            if varargin{1} == -1 %should make this work regardless of 2D or 3D array, etc
                varargin{1} = length(obj_out.value);
            end
            obj_out.value = obj_out.value(varargin{:});
        end
        
        function obj_out = change_unit(obj,desUnit,varargin)
            % obj_out = change_unit(obj,desDim,varargin)
            % obj: input c_dim_quant
            % desUnit: string containing desired unit
            % varargin:
            %   (...,sysModel) contains a c_system_model for nondimensionalization
            %   (...,spacecraft) allows only for mass nondimensionalization
            %   (...,[lstar tstar mstar]) allows for manual nd
            %   (...,sysModel,spacecraft) allows for all nd, with mass of
            %       spacecraft instead of system mstar
            obj_out = obj;
            if strcmp(desUnit,obj.unit)
            %do nothing, already in correct units
            else
            unitTable = {
              %unit     type        to km,s,kg    
              'mm',     'length',   1e-6;
              'cm',     'length',   1e-5;
              'm',      'length',   1e-3;
              'km',     'length',   1;
              'sec',    'time',     1; 
              'min',    'time',     60;
              'hour',   'time',     60*60;
              'day',    'time',     60*60*24;
              'year',   'time',     60*60*24*365;
              'km/s',   'speed',    1;
              'm/s',    'speed',    1e-3;
              'km/s2',  'accel',    1;
              'm/s2',   'accel',    1e-3;
              'mN',     'force',    1e-6;
              'N',      'force',    1e-3;
              'kN',     'force',    1; %kg*km/s^2
              'g',      'mass',     1e-3;
              'kg',     'mass',     1;
              'deg',    'angle',    1;
              'rad',    'angle',    180/pi;
              'kg/s',   'mdot',     1;
              'km3/s2', 'mu',       1; 
            };
            switch nargin
                case 2
                case 3
                    argin = varargin{1};
                    if isa(argin,'c_system_model')
                        sysModel = argin;
                        lstar = sysModel.char.lstar.value;
                        tstar = sysModel.char.tstar.value;
                        mstar = sysModel.char.mstar.value;
                        unitTableND = {
                            %unit       %type       %to km,s,kg
                            'nd_l',     'length',   lstar;
                            'nd_t',     'time',     tstar;
                            'nd_m',     'mass',     mstar;
                            'nd_a',     'accel',    lstar/tstar^2;
                            'nd_v',     'speed',    lstar/tstar;
                            'nd_av',    'angVel',   1/tstar;
                            'nd_f',     'force',    mstar*lstar/tstar^2;
                            'nd_md',    'mdot',     mstar/tstar;
                            'nd_mu',    'mu',       lstar^3/tstar^2;
                        };
                    elseif isa(argin,'c_spacecraft')
                        mstar = argin.M0.value;
                        unitTableND = {
                            %unit       %type       %to km,s,kg
                            'nd_m',     'mass',     mstar;
                        };
                    elseif length(argin) == 3
                        lstar = argin(1);
                        tstar = argin(2);
                        mstar = argin(3);
                        unitTableND = {
                            %unit       %type       %to km,s,kg
                            'nd_l',     'length',   lstar;
                            'nd_t',     'time',     tstar;
                            'nd_m',     'mass',     mstar;
                            'nd_a',     'accel',    lstar/tstar^2;
                            'nd_v',     'speed',    lstar/tstar;
                            'nd_av',    'angVel',   1/tstar;
                            'nd_f',     'force',    mstar*lstar/tstar^2;
                            'nd_md',    'mdot',     mstar/tstar;
                            'nd_mu',    'mu',       lstar^3/tstar^2;
                        };
                    else
                        error('unknown system characteristics')
                    end
                    unitTable = [unitTable; unitTableND]; %append nd
                case 4
                    argin1 = varargin{1};
                    argin2 = varargin{2};
                    if isa(argin1,'c_system_model') && isa(argin2,'c_spacecraft')
                        sysModel = argin1;
                        lstar = sysModel.char.lstar.value;
                        tstar = sysModel.char.tstar.value;
                        mstar = argin2.M0.value;
                        unitTableND = {
                            %unit       %type       %to km,s,kg
                            'nd_l',     'length',   lstar;
                            'nd_t',     'time',     tstar;
                            'nd_m',     'mass',     mstar;
                            'nd_a',     'accel',    lstar/tstar^2;
                            'nd_v',     'speed',    lstar/tstar;
                            'nd_av',    'angVel',   1/tstar;
                            'nd_f',     'force',    mstar*lstar/tstar^2;
                            'nd_md',    'mdot',     mstar/tstar;
                            'nd_mu',    'mu',       lstar^3/tstar^2;
                        };
                    end
                    unitTable = [unitTable; unitTableND]; %append nd
                otherwise
                    error('Too many input arguments')
            end
            inUnit = obj.unit;
            ind_in = strcmp(unitTable(:,1),inUnit);
            if ~any(ind_in)
                if strcmp(inUnit,'s')
                    warning('Did you mean "sec"?')
                end
                error(['The input unit, ' inUnit ', does not exist in the lookup table'])
            end
            ind_out = strcmp(unitTable(:,1),desUnit);
            if ~any(ind_out)
                if strcmp(desUnit,'s')
                    warning('Did you mean "sec"?')
                end
                error(['The output unit, ' desUnit ', does not exist in the lookup table'])
            end
            type_in = unitTable{ind_in,2};
            type_out = unitTable{ind_out,2};
            if ~strcmp(type_in,type_out)
                error(['The requested units do not match. Input unit is of type: '...
                    type_in ', while output unit is of type: ' type_out])
            end
            %convert from input to SI
            stUnit = unitTable{ind_in,3} .* obj.value;
            %convert from SI to desUnit
            obj_out.value = stUnit ./ unitTable{ind_out,3};
            obj_out.unit = desUnit;
            end
        end

        %% ======================= overload ==================================
        %operations performed on dimensional quantities will turn them into doubles
        function out = plus(a,b) 
            if isa(b,'c_dim_quant')
                b = b.value;
            end
            if isa(a,'c_dim_quant')
                a = a.value;
            end
            out = a + b;
        end

        function out = add(a,b)
            if ~isa(a,'c_dim_quant') || ~isa(b,'c_dim_quant')
                error('both inputs must be c_dim_quant')
            end
            if ~strcmp(a.unit,b.unit)
                error('units must match between inputs')
            end
            out = a;
            out.value = a.value + b.value;
        end

        function out = minus(a,b) 
            if isa(b,'c_dim_quant')
                b = b.value;
            end
            if isa(a,'c_dim_quant')
                a = a.value;
            end
            out = a - b;
        end

        function out = sub(a,b)
            if ~isa(a,'c_dim_quant') || ~isa(b,'c_dim_quant')
                error('both inputs must be c_dim_quant')
            end
            if ~strcmp(a.unit,b.unit)
                error('units must match between inputs')
            end
            out = a;
            out.value = a.value - b.value;
        end

        function out = uplus(a) 
            out = a.value;
        end

        function out = uminus(a) 
            out = -a.value;
        end

        function out = times(a,b) 
            if isa(b,'c_dim_quant')
                b = b.value;
            end
            if isa(a,'c_dim_quant')
                a = a.value;
            end
            out = a .* b;
        end

        function out = mult(a,b)
            if ~isa(a,'c_dim_quant') || ~isa(b,'c_dim_quant')
                error('both inputs must be c_dim_quant')
            end
            if ~strcmp(a.unit,b.unit)
                error('units must match between inputs')
            end
            out = a;
            out.value = a.value .* b.value;
        end

        function out = mtimes(a,b) 
            if isa(b,'c_dim_quant')
                b = b.value;
            end
            if isa(a,'c_dim_quant')
                a = a.value;
            end
            out = a * b;
        end

        function out = rdivide(a,b) 
            if isa(b,'c_dim_quant')
                b = b.value;
            end
            if isa(a,'c_dim_quant')
                a = a.value;
            end
            out = a ./ b;
        end

        function out = ldivide(a,b) 
            if isa(b,'c_dim_quant')
                b = b.value;
            end
            if isa(a,'c_dim_quant')
                a = a.value;
            end
            out = a .\ b;
        end

        function out = mrdivide(a,b) 
            if isa(b,'c_dim_quant')
                b = b.value;
            end
            if isa(a,'c_dim_quant')
                a = a.value;
            end
            out = a / b;
        end

        function out = mldivide(a,b) 
            if isa(b,'c_dim_quant')
                b = b.value;
            end
            if isa(a,'c_dim_quant')
                a = a.value;
            end
            out = a \ b;
        end

        function out = power(a,b) 
            if isa(b,'c_dim_quant')
                b = b.value;
            end
            if isa(a,'c_dim_quant')
                a = a.value;
            end
            out = a .^ b;
        end

        function out = mpower(a,b) 
            if isa(b,'c_dim_quant')
                b = b.value;
            end
            if isa(a,'c_dim_quant')
                a = a.value;
            end
            out = a ^ b;
        end

        function out = lt(a,b) 
            if isa(b,'c_dim_quant')
                b = b.value;
            end
            if isa(a,'c_dim_quant')
                a = a.value;
            end
            out = a < b;
        end

        function out = gt(a,b) 
            if isa(b,'c_dim_quant')
                b = b.value;
            end
            if isa(a,'c_dim_quant')
                a = a.value;
            end
            out = a > b;
        end

        function out = le(a,b) 
            if isa(b,'c_dim_quant')
                b = b.value;
            end
            if isa(a,'c_dim_quant')
                a = a.value;
            end
            out = a <= b;
        end

        function out = ge(a,b) 
            if isa(b,'c_dim_quant')
                b = b.value;
            end
            if isa(a,'c_dim_quant')
                a = a.value;
            end
            out = a >= b;
        end

        function out = ne(a,b) 
            if isa(b,'c_dim_quant')
                b = b.value;
            end
            if isa(a,'c_dim_quant')
                a = a.value;
            end
            out = a ~= b;
        end

        function out = eq(a,b) 
            if isa(b,'c_dim_quant')
                b = b.value;
            end
            if isa(a,'c_dim_quant')
                a = a.value;
            end
            out = a == b;
        end

        function out = ctranspose(a) 
            out = a.value';
        end

        function out = transpose(a) 
            out = a.value.';
        end

        function out = horzcat(a,varargin)
            vararg = cell(1,length(varargin));
            for i = 1:length(varargin)
                if isa(varargin{i},'c_dim_quant')
                    b = varargin{i};
                    b = b.change_unit(a.unit);
                else
                    error('cannot horzcat a c_dim_quant and something else')
                end
                vararg{i} = b.value;
            end
            outval = horzcat(a.value,vararg{:});
            out = c_dim_quant(outval,a.unit);
        end
        
        function out = vertcat(a,varargin)
            vararg = cell(1,length(varargin));
            for i = 1:length(varargin)
                if isa(varargin{i},'c_dim_quant')
                    b = varargin{i};
                    b = b.change_unit(a.unit);
                else
                    error('cannot horzcat a c_dim_quant and something else')
                end
                vararg{i} = b.value;
            end
            outval = vertcat(a.value,vararg{:});
            out = c_dim_quant(outval,a.unit);
        end
    end
end

