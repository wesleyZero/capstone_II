
function fxns = reactor_fxns() 
    fxns.get_reaction_rate = @get_reaction_rate;
    fxns.get_reactor_flows = @get_reactor_flows;
end

function [F, P, R] = get_reactor_flows(F, tau, T, P, opt)
    % basis calculations  
    q = get_volumetric_flowrates(F, T, P, opt);
    q_total = get_total_volumetric_flowrate(q).value;
    V_rxtr.basis = q_total * tau;
    
    C_init = get_initial_concentrations(F, V_rxtr.basis, tau);


    F = NaN; P = NaN; R = NaN;
end

function C = get_initial_concentrations(F, V_rxtr, tau);
    % F.mol = mol /s 
    % V = L
    % tau = s

    fieldNames = fieldnames(F);
    for i = 1:length(fieldNames)
        C.(fieldNames{i}) = F.(fieldNames{i}).mol * tau / V_rxtr;
    end
end 

function condition = get_condition(T, P, opt)
    switch opt
        case 'isobaric'
            condition = T;
        case 'isothermal'
            condition = P; 
        otherwise
            condition = NaN; 
            disp("ERROR : get_volumetric_flowrates : opt not valid");
    end
end

function rho = get_all_densities(T, P, opt)
    % rho = const.densities; 
    % rho.methanol = get_methanol_density(T, P);
    % rho.carbon_dioxide = get_supercritical_c02_density(condition, opt, P);

    condition = get_condition(T, P, opt);
    rho = get_constants().densities; 
    rho.methanol = get_methanol_density(T, P);
    rho.carbon_dioxide = get_supercritical_c02_density(condition, opt, P);
    
end

function q_total = get_total_volumetric_flowrate(q)
    q_total.value = 0;
    fieldNames = fieldnames(q);
    for i = 1:length(fieldNames)
        if fieldNames{i} == "units"
            q_total.units = q.units; 
            continue
        end
        q_total.value = q_total.value +  q.(fieldNames{i});
    end
end

function q = get_volumetric_flowrates(F, T, P, opt)
    const = get_constants();
    % condition = get_condition(T, P, opt);
    % rho = const.densities; 
    % rho.methanol = get_methanol_density(T, P);
    % rho.carbon_dioxide = get_supercritical_c02_density(condition, opt, P);
    rho = get_all_densities(T, P, opt);    

    q.units = "m^3 / s";
    % Carbon dioxide 
    % (L / s)        = (mol / s)            * (g / mol)                      
    q.carbon_dioxide = F.carbon_dioxide.mol * const.molar_mass.carbon_dioxide * ... 
        ...% (kg / g)        * (m^3 / kg)             * (L / m^3) 
        const.units.mass.kg_per_g * (1/rho.carbon_dioxide) * const.units.volume.l_per_m3;
    % Methanol
    q.methanol = F.methanol.mol * const.molar_mass.methanol * ...
        const.units.mass.kg_per_g * (1 / rho.methanol) * const.units.volume.l_per_m3;
    % Ethylene Carbonate
    q.ethylene_carbonate = F.ethylene_carbonate.mol * const.molar_mass.ethylene_carbonate *...
         const.units.mass.kg_per_g * (1 / rho.ethylene_carbonate) * const.units.volume.l_per_m3;
    
    q.carbon_dioxide = q.carbon_dioxide * const.units.volume.l_per_m3;
    q.methanol = q.methanol * const.units.volume.l_per_m3;
    q.ethylene_carbonate = q.ethylene_carbonate * const.units.volume.l_per_m3;
    q.units = "L / s" ;
    
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




% function V = get_reactor_volume()
%     V = 1; % ?? 
% end

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

    switch reaction
        case '2f'
            rho = get_supercritical_c02_density(P, 'isothermal');
            if rho > 246.82 % [g / L]
                k = 0.02486 - 4.943 * (10^(-5)) * rho;
            else
                k = 0.013;
            end
        case '2r'
            k = 0.01486 * rho^(-0.873);
        case '3'
            k = 3.014 * (10^(-4)) * exp(-5.99 * (10^(-3)) * rho);
        otherwise
            k = NaN;
            disp("ERROR: get_isothermal_rate_constant(): invalid reaction option")
    end


end

function rho = get_supercritical_c02_density(condition, opt, pressure)
    % Input: condition = T or P. Depending on option
    %   P [=] bar
    %   T [=] celcius   
    % Assumptions:
    %   Isobaric model is at 140 C
    %   Isothermal model is at 150 bar
    % Ranges of input
    %   P = [50 bar, 150 bar]
    %   T = [80 C, 140 C]
    % Output: 
    %   rho [=] kg / m^3
    withinTempRange = @(T) T >= 80 && T <= 140;
    withinPressureRange = @(P) P >= 50 && P <= 150;
    
    rho.units = "kg / m^3";
    switch opt
        case 'isothermal'
            P = condition;
            if withinPressureRange(P)
                rho = 1.6746 * P - 12.592;
                    % NIST / Excel Regression
            else
                rho = NaN;
                disp("get_supercritical_c02_density : ERROR : P out of range")
            end
        case 'isobaric'
            T = condition;
            if withinTempRange(T)
                if pressure < 125 % [C]
                    rho = 356.08 * exp(-0.006 * T);
                        % NIST Data at 100 bar
                else
                    rho = 838.87 * exp(-0.009 * T);
                        % NIST Data at 150 bar
                end
            else
                rho = NaN;
                disp("get_supercritical_c02_density : ERROR : T out of range")
            end
        otherwise
            disp("SUPERCRITICAL C02 DENSITY FUNCTION ERROR: invalid opt")
            rho = NaN;
    end
end

function rho = get_methanol_density(T, P)
    % Input:
    %   P [=] bar
    %   T [=] celcius   
    % Assumptions:
    %   Isobaric model is at 140 C
    %   Isothermal model is at 150 bar
    % Ranges of input
    %   P = [50 bar, 150 bar]
    %   T = [80 C, 140 C]
    % Output: 
    %   rho [=] kg / m^3
    withinTempRange = @(T) T >= 80 && T <= 140;
    % withinPressureRange = @(P) P >= 50 && P <= 150;
    rho.units = "kg / m^3";  
    if withinTempRange(T)%  && withinPressureRange(P)
        if P < 125; % [C]
            rho = -1.0866 * T + 833.31;
                % NIST Data at 100 bar
        else
            rho = -1.0354 * T + 834.79;
                % NIST Data at 150 bar
        end
    else
        rho = NaN;
        disp("get_methanol_density : ERROR : T or P out of range")
    end

end



