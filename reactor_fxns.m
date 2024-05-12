
function fxns = reactor_fxns() 
    % fxns.get_supercritical_c02_density = @get_supercritical_c02_density;

end

function r = get_reaction_rate(reaction, condition, opt, F)
    % input:
    %   condition = T or P
    %   opt = 'isothermal' or 'isobaric'

    C = get_concentrations();
    k = get_rate_constant(reaction, condition, opt);
    switch reaction
        case '2f'
            r = k * (C.ethylene_carbonate)^0.8;
        case '2r'
            r = k * C.dimethyl_carbonate * C.ethylene_carbonate;
        case '3'
            r = k * C.ethylene_carbonate;
        otherwise
            r = NaN;
            disp("ERROR: get_reaction_rate(): invalid reaction option")
    end
end

function C = get_concentrations()
    C.ethylene_carbonate = 1; % ?? 
    C.dimethyl_carbonate = 1; % ?? 
end 


function V = get_reactor_volume()
    V = 1; % ?? 
end

function k = get_rate_constant(reaction, condition, opt)

    if opt == 'isothermal'
        P = condition;
        k = get_isothermal_rate_constant(reaction, P);
    elseif opt == 'isobaric'
        T = condition;
        k = get_isobaric_rate_constant(reaction, T);
    else
        k = NaN;
        disp("ERROR: get_rate_constant(): invalid opt")
    end
end

function k = get_isobaric_rate_constant(reaction, T)
    % input: 
    %   T [ C ]
    % output:
    %   k [mol / L s ]

    const = get_constants();
    thermo = const.thermo;
    T = const.units.temperature.c_to_k(T); % [ K ]

    switch reaction
        case '2f'
            k = 6.69 * 10^2 * exp(-37200 / (thermo.R * T));
        case '2r'
            k = 1.19 * 10^4 * exp(-53700 / (thermo.R * T));
        case '3'
            k = 1.89 * 10^6 * exp(-82400 / (thermo.R * T));
        otherwise
            k = NaN;
            disp("ERROR: get_isobaric_rate_constant(): invalid reaction option");
    end
end

function k = get_isothermal_rate_constant(reaction, P)
    % input:
    %   P [ bar ]

    if reaction == '2f'
        rho = get_supercritical_c02_density(P, 'isothermal')
        if rho > 246.82 % [g / L] 
            k = 0.02486 - 4.943 * (10^(-5)) * rho; 
        else
            k = 0.013; 
        end
    elseif reaction == '2r'
        k = 0.01486 * rho^(-0.873);
    elseif reaction == '3'
        k = 3.014 * (10^(-4)) * exp(-5.99 * (10^(-3)) * rho);
    else 
        k = NaN;
        disp("ERROR: get_isothermal_rate_constant(): invalid reaction option")
    end

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
    elseif opt == 'isobaric'
        T = condition;
        T = get_constants().units.temperature.c_to_k(T); % [ K ]
        rho.kg_m3 = 0.037 * T^2 - 10.46 * T + 1060.8;
    else
        disp("SUPERCRICIAL C02 DENSITY FUNCTION ERROR: invalid opt")
        rho = NaN;
    end

end



