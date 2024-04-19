classdef c_elem
    %   c_elem(mu,sma,ecc,inc,raan,argP,trueAnom)
            %
            %   if you don't have exactly the 6 standard elements, then you
            %   must input any 6 valid elements as varargin
            %   c_elem(mu,varargin)
            %       where vararign:
            %       [must include any 2 from first 4, inc, and any 3 after]
            %   (...,'period',period,...) is the period of the orbit
            %   (...,'sma',sma,...) is the semi-major axis
            %   (...,'ecc',ecc,...) is the eccentricity
            %   (...,'slr',p,...) is the semi-latus rectum
            %   (...,'inc',inc,...) is the inclination
            %   (...,'raan',raan,...) is the right ascension of the ascending node
            %   (...,'argP',argP,...) is the argument of periapis
            %   (...,'trueAnom',trueAnom,...) is the true anomaly
            %   (...,'meanAnom',meanAnom,...) is the mean anomaly
            %   (...,'eccAnom',eccAnom,...) is the eccentric anomaly
            %   (...,'meanLong',meanLong,...) is the mean longitude
            %   (...,'argLat',argLat,...) is the argument of latitude
            %   (...,'trueLong',trueLong,...) is the true longitude
            %   (...,'longPeri',longPeri,...) is the longitude of periapsis
    
    properties
        mu c_dim_quant
        period c_dim_quant
        sma c_dim_quant
        ecc c_dim_quant
        semilatusrectum c_dim_quant
        inc c_dim_quant
        raan c_dim_quant
        argP c_dim_quant
        trueAnom c_dim_quant
        meanAnom c_dim_quant
        eccAnom c_dim_quant
        meanLong c_dim_quant %argP + raan + meanAnom
        argLat c_dim_quant %argP + tAnom
        trueLong c_dim_quant %raan + argP + tAnom
        longPeri c_dim_quant %raan + argP
    end
    
    methods
        function obj = c_elem(mu,varargin)
            %   c_elem(mu,sma,ecc,inc,raan,argP,trueAnom)
            %
            %   if you don't have exactly the 6 standard elements, then you
            %   must input any 6 valid elements as varargin
            %   c_elem(mu,varargin)
            %       where vararign:
            %       [must include any 2 from first 4, inc, and any 3 after]
            %   (...,'period',period,...) is the period of the orbit
            %   (...,'sma',sma,...) is the semi-major axis
            %   (...,'ecc',ecc,...) is the eccentricity
            %   (...,'slr',p,...) is the semi-latus rectum
            %   (...,'inc',inc,...) is the inclination
            %   (...,'raan',raan,...) is the right ascension of the ascending node
            %   (...,'argP',argP,...) is the argument of periapis
            %   (...,'trueAnom',trueAnom,...) is the true anomaly
            %   (...,'meanAnom',meanAnom,...) is the mean anomaly
            %   (...,'eccAnom',eccAnom,...) is the eccentric anomaly
            %   (...,'meanLong',meanLong,...) is the mean longitude
            %   (...,'argLat',argLat,...) is the argument of latitude
            %   (...,'trueLong',trueLong,...) is the true longitude
            %   (...,'longPeri',longPeri,...) is the longitude of periapsis
            if strcmp(mu.unit,'km3/s2')
                smaUnit = 'km';
            elseif strcmp(mu.unit,'nd_mu')
                smaUnit = 'nd_l';
            else
                error('unknown unit for mu')
            end
            obj.mu = mu;
            switch nargin
                case 7
                    [sma, ecc, inc, raan, argP, trueAnom] ...
                        = deal(varargin{:});
                    obj.sma = sma.change_unit(smaUnit);
                    obj.ecc = ecc;
                    obj.inc = inc.change_unit('rad');
                    obj.raan = raan.change_unit('rad');
                    obj.argP = argP.change_unit('rad');
                    obj.trueAnom = trueAnom.change_unit('rad');
                    period = 2*pi*sqrt(sma.value^3/mu);
                    obj.period = c_dim_quant(period,'nd_t');
                case 13
                    for i = 1:length(varargin)/2
                        switch varargin{i*2-1}
                            case 'sma'
                                obj.sma = varargin{i*2}.change_unit(smaUnit);
                            case 'period'
                                obj.period = varargin{i*2}.change_unit('sec');
                            case 'ecc'
                                obj.ecc = varargin{i*2};
                            case 'slr'
                                obj.semilatusrectum = varargin{i*2}.change_unit(smaUnit);
                            case 'inc'
                                obj.inc = varargin{i*2}.change_unit('deg');
                            case 'raan'
                                obj.raan = varargin{i*2}.change_unit('deg');
                            case 'argP'
                                obj.argP = varargin{i*2}.change_unit('deg');
                            case 'trueAnom'
                                obj.trueAnom = varargin{i*2}.change_unit('deg');
                            case 'meanAnom'
                                obj.meanAnom = varargin{i*2}.change_unit('deg');
                            case 'eccAnom'
                                obj.eccAnom = varargin{i*2}.change_unit('deg');
                            case 'meanLong'
                                obj.meanLong = varargin{i*2}.change_unit('deg');
                            case 'argLat'
                                obj.argLat = varargin{i*2}.change_unit('deg');
                            case 'trueLong'
                                obj.trueLong = varargin{i*2}.change_unit('deg');
                            case 'longPeri'
                                obj.longPeri = varargin{i*2}.change_unit('deg');
                        end
                    end
                    mu = obj.mu.value;
                    % get all shape parameters
                        if isempty(obj.sma) && isempty(obj.period)
                            p = obj.semilatusrectum.value;
                            e = obj.ecc.value;
                            sma = p/(1-e^2);
                            obj.sma = c_dim_quant(sma,smaUnit);
                            period = 2*pi*sqrt(sma^3/mu);
                            obj.period = c_dim_quant(period,'sec');
                        else
                            if isempty(obj.sma)
                                period = obj.period.value;
                                sma = (mu*(period/2/pi)^2)^(1/3);
                                obj.sma = c_dim_quant(sma,smaUnit);
                            else
                                sma = obj.sma.value;
                                period = 2*pi*sqrt(sma^3/mu);
                                obj.period = c_dim_quant(period,'sec');
                            end
                            if isempty(obj.ecc)
                                slr = obj.semilatusrectum.value;
                                ecc = sqrt(1-slr/sma);
                                obj.ecc = c_dim_quant(ecc,'');
                            else
                                ecc = obj.ecc.value;
                                slr = sma*(1-ecc^2);
                                obj.semilatusrectum = c_dim_quant(slr,smaUnit);
                            end
                        end
                            
                    % get all angle parameters

                otherwise
                    error('must have either 6 standard elements or a 12-length list of other elements')
            end
        end
    end
end


function M = meanAnomaly(meanLong,argP,RAAN)
M = wrapTo2Pi(meanLong - argP - RAAN);
end
