
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
    condition = get_condition(T, P, opt);
    C = get_reactor_effluent_concentrations(C_init, tau, T, P, opt);


    F = NaN; P = NaN; R = NaN;
end

function C_vector = get_conc_vector(C)
    C_vector(1) = C.ethylene_carbonate;
    C_vector(2) = C.ethylene_glycol;
    C_vector(3) = C.methanol;
    C_vector(4) = C.carbon_dioxide;
    C_vector(5) = C.dimethyl_carbonate;
    C_vector(6) = C.methoxy_ethanol;
end

function k = get_all_rate_constants(condition, opt)
    k.k2f = get_rate_constant('2f', condition, opt);
    k.k2r = get_rate_constant('2r', condition, opt);
    k.k3 = get_rate_constant('3', condition, opt);
end

function b = get_sys_eqns_constants(T, P, opt)
    % This b vector are what the system of non-linear equations are equal to.
    % I used 'b' because it's like a non-linear version of Ax=b

    user = get_user_inputs(); 
    rho = get_all_molar_densities(T, P, opt); 

    sum = 0;
    sum = user.level3.molar_ratio_carbon_dioxide_EO / rho.carbon_dioxide;
    sum = sum + user.level3.molar_ratio_methanol_EO / rho.methanol;
    sum = sum + 1 / rho.carbon_dioxide;
    sum_recip = 1 / sum;

    % Note: this calculation is kind of weird, because of the way that we are 
    % modeling this 'virtual reactor' before the first 'real' reactor. Reminder 
    % for myself in the future and others is that the kinetics of the first 
    % reaction are so fast that we are assuming that first reaction goes 
    % to completion. This first reaction is the 'virtual reactor'. When 
    % calculating the 'real' reactor feed we are using the molar ratios to 
    % ethylene oxide, however we don't actually have any ethylene oxide
    % because it ran to completion with C02 in massive stiochiometric excess. 
    % So we are using the ethylene carbonotate instead. 

    b(1) = sum_recip;
    b(2) = -user.level3.molar_ratio_methanol_EO * sum_recip;
    b(3) = user.level3.molar_ratio_carbon_dioxide_EO * sum_recip;
    b(4) = 0;
    b(5) = 0;
    b(6) = 0;

end

function C = get_reactor_effluent_concentrations(C_init, tau, T, P, opt)
    C_init_vector = get_conc_vector(C_init);
    condition = get_condition(T, P, opt);
    k = get_all_rate_constants(condition, opt);
    b = get_sys_eqns_constants(T, P, opt);
    eqns = @(C) sys_of_eqns(C, k, b, tau);

    % C = fsolve() 

    C = NaN;
end

function F = sys_of_eqns(C, k, b, tau)
    F(1) = -C(1) - tau * C(1)^0.8 - tau * k.k3 * C(1) + tau * k.k2r * C(5) * C(2) - b(1); 
    F(2) = -tau * k.k2f * C(1)^0.8 - tau * k.k3 * C(1) - C(3) + tau * k.k2r * C(5) * C(2) + b(2);
    F(3) = -C(4) + tau * k.k3 * C(1) - b(3);
    F(4) = -C(5) + tau * k.k2f * C(1)^0.8 - tau * k.k2r * C(5) * C(2) - b(4);
    F(5) = C(2) + tau * k.k2f * C(1)^0.8

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
            disp("ERROR : get_condition : opt not valid");
    end
end

function rho = get_all_molar_densities(T, P, opt)
    rho = get_all_densities(T, P, opt);
    const = get_constants();

    fieldNames = fieldnames(rho);
    for i = 1:length(fieldNames)
        if strcmp(fieldNames{i},'units')
            rho.(fieldNames{i}) = 'mol / L';
            continue
        end

        % mol / L              = (kg / m^3)              * (g / kg)
        rho.(fieldNames{i}) = rho.(fieldNames{i}) * const.units.mass.g_per_kg * ... 
        ...% (mol / g)                            * (m^3 / L) 
        (1 / const.molar_mass.(fieldNames{i})) * const.units.volume.m3_per_l;
    end
end

function rho = get_all_densities(T, P, opt)
    % out is in kg / m^3

    condition = get_condition(T, P, opt);
    rho = get_constants().densities; 
    rho.methanol = get_methanol_density(T, P);
    rho.carbon_dioxide = get_supercritical_c02_density(condition, opt, P);
    
end

function q_total = get_total_volumetric_flowrate(q)
    % output [L / s]
    q_total.value = 0;
    fieldNames = fieldnames(q);
    for i = 1:length(fieldNames)
		if strcmp(fieldNames{i},'units')
            q_total.units = q.units; 
            continue
        end
        q_total.value = q_total.value +  q.(fieldNames{i});
    end
end

function q = get_volumetric_flowrates(F, T, P, opt)
    const = get_constants();
    rho = get_all_densities(T, P, opt);    

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

    q.units = 'L / s' ;
    
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

	if strcmp(opt, 'isothermal')
    	P = condition;
    	k = get_isothermal_rate_constant(reaction, P);
	elseif strcmp(opt, 'isobaric')
    	T = condition;
    	k = get_isobaric_rate_constant(reaction, T);
	else
    	k = NaN;
    	disp("ERROR: get_rate_constant(): invalid opt");
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
    
    rho.units = 'kg / m^3';
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



