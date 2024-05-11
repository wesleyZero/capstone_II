
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
    elseif opt == 'isobaric'
        T = condition;
        rho.kg_m3 = 0.037 * T^2 - 10.46 * T + 1060.8;
    else
        disp("SUPERCRICIAL C02 DENSITY FUNCTION ERROR: invalid opt")
        rho = NaN;
    end

end


function k = get_rate_constant(reaction, condition, opt)

    if opt == 'isothermal'
        P = condition;
        k = get_isothermal_rate_constant(reaction, P);
    elseif opt == 'isobaric'
        T = condition;
        k = get_isobaric_rate_constant(reaction, T);
    else
        disp("ERROR: get_rate_constant(): invalid opt")
    end
end

function k = get_isobaric_rate_constant(reaction, T)
    % input: 
    %   T [ K ]
    % output:
    %   k [mol / L s ]

    thermo = get_constants().thermo;

    switch reaction
        case '2f'
            k = 6.69 * 10^2 * exp(-37200 / (thermo.R * T));
        case '2r'
            k = 1.19 * 10^4 * exp(-53700 / (thermo.R * T));
        case '3'
            k = 1.89 * 10^6 * exp(-82400 / (thermo.R * T));
        otherwise
            disp("ERROR: get_isobaric_rate_constant(): invalid reaction option");
    end

end






























% function isoT = get_isothermal_rate_constants() 

%     isoT.unitt = "";
%     isoT.k2f_lowRho = 0.013; 
%     isoT.k2f_highRho = @(rho) 0.02486 - 4.943 * 10^-5 * rho;
%     isoT.k3 = @(rho) 3.014 * 10^-4 * exp(-5.99 * rho);
% end

% function rate_const = get_rate_constants() 

%     rate_const.isoP.units = "";
%     rate_const.isoP.k2f = @(T_K) 6.59 * 10^2 * exp( -37000 / thermo.R / T_K);
%     rate_const.isoP.k2r = @(T_K) 1.19 * 10^4 * exp( -53700 / thermo.R / T_K);
%     rate_const.isoP.k3 = @(T_K) 1.89 * 10^6 * exp( -82400 / thermo.R / T_K);
% end

% function rate = get_reaction_rates()

%     rate.units = "";
%     rate.r2f = @(T_K,C_EC) thermo.rate_const.k2f(T_K) * C_EC^0.8;
%     rate.r2r = @(T_K,C_DMC, C_EG) thermo.rate_const.k2r(T_K) * C_DMC * C_EG;
%     rate.r3 = @(T_K, C_EC) thermo.rate_const.k3(T_K) * C_EC;
% end
