
function fxns = reactor_fxns() 
    fxns.get_supercritical_c02_density = @get_supercritical_c02_density;

end

function rho = get_supercritical_c02_density(condition, opt)
    % Input: condition = T or P. Depending on option
    %   P [=] bar
    %   T [=] celcius   
    % Assumptions:
    %   Isobaric model is at 140 C
    %   Isothermal model is at 150 bar
    % Ranges of input
    %   P = [150 bar, 200 bar] ???
    %   T = [30 C, 150 C] ???
    
    if opt == 'isothermal'
        P = condition;
        rho.kg_m3 = 1.8853 * P - 31.755;
    else if opt == 'isobaric'
        T = condition;
        rho.kg_m3 = 0.037 * T^2 - 10.46 * T + 1060.8;
    else
        disp("SUPERCRICIAL C02 DENSITY FUNCTION ERROR: invalid opt")
        rho = NaN;
    end

end