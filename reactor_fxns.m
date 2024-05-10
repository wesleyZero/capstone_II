


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
        rho = 
    else if opt == 'isobaric'

    else
        disp("SUPERCRICIAL C02 DENSITY FUNCTION ERROR: invalid opt")
    end

end