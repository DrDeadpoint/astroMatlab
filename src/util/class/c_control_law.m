classdef c_control_law
    %c_control_law : gives details on low-thrust spacecraft control laws
    %law could be SOC or OCT or SAR or fixed_dir
    %frame is what reference frame the law is defined in
    %coeffs are any necessary coefficients for the given LAW
    
    properties
        law char
        frame
        coeffs 
    end
    
    methods
        function obj = c_control_law(law,frame,coeffs)
            %obj = c_control_law(law,frame,coeffs)
            %ensure that rates of angle changes are nondimensionalized correctly
            obj.law = law;
            obj.frame = frame;
            switch law
                case 'SOC'
                    my_fields = {'alpha0','beta0','alphadot','betadot','alphaddot','betaddot',...
                                'alphaAmp','betaAmp','alphaFreq','betaFreq','alphaPhase','betaPhase'};
                    my_coeffs = struct;
                    for i = 1:length(my_fields)
                        fi = my_fields{i};
                        if isfield(coeffs,fi)
                            my_coeffs.(fi) = coeffs.(fi);
                        else
                            my_coeffs.(fi) = 0; %remove that portion of the equation
                        end
                    end
                case 'fixed_dir'
                    if abs(norm(coeffs) - 1) > 1e-5
                       error('Thrust direction must be unit vector')
                    end
                    my_coeffs = coeffs;
                otherwise
                    error(['Your designated control law, ' law ', has not been implemented.'])
            end
            obj.coeffs = my_coeffs;
        end
    end
end

