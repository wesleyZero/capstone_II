
function fxns = reactor_fxns() 
    fxns.get_reaction_rate = @get_reaction_rate;
    fxns.get_reactor_flows = @get_reactor_flows;
end

function [F_in, F_out, R] = get_reactor_flows(F_basis, T, P, opt, tau)
    % basis calculations  
    q_tot = get_total_volumetric_flowrate(F_basis, T, P, opt);
    V_rxtr.basis = q_tot * tau;
    V_rxtr = get_reactor_volume(F_basis, T, P, opt, tau);
    
    C = get_reactor_effluent_concentrations(F_basis, T, P, opt, tau);

    if any(imag(C) ~= 0)
        disp('ERROR : Complex valued concentrations');
%         C
    elseif any(real(C) < 0)
        disp('ERROR : Negative valued concentrations');
%         C
    else
        disp('Valid solution!!!!!!!!!!')
    end
    F_in = NaN; F_out = NaN; R = NaN;
end

function C_vector = get_conc_vector(F, T, P, opt, tau)
    C = get_concentrations(F, T, P, opt, tau);
    C_vector(1) = C.ethylene_carbonate;
    C_vector(2) = C.ethylene_glycol;
    C_vector(3) = C.methanol;
    C_vector(4) = C.carbon_dioxide;
    C_vector(5) = C.dimethyl_carbonate;
    C_vector(6) = C.methoxy_ethanol;
end

function k = get_all_rate_constants(T, P, opt)
    k.k2f = get_rate_constant('2f', T, P, opt);
    k.k2r = get_rate_constant('2r', T, P, opt);
    k.k3 = get_rate_constant('3', T, P, opt);
end

function C = get_reactor_effluent_concentrations(F, T, P, opt, tau)
    user = get_user_inputs();
    % Ci0 = get_conc_vector(F, T, P, opt, tau);
    Ci0 = get_concentrations(F, T, P, opt, tau);
    Ci_init = get_conc_vector(F, T, P, opt, tau);
    params.Ci0 = Ci0;
    params.tau = tau;
    params.T = T;
    params.P = P;
    params.opt = opt;

    eqns = @(C) sys_of_eqns(C, params);
    C = fsolve(eqns, Ci_init, user.level3.fsolveOpt);
end

function C_struct = get_concentration_struct(C_vector)
    C_struct.ethylene_carbonate = C_vector(1);
    C_struct.ethylene_glycol = C_vector(2);
    C_struct.methanol = C_vector(3);
    C_struct.carbon_dioxide = C_vector(4);
    C_struct.dimethyl_carbonate = C_vector(5);
    C_struct.methoxy_ethanol = C_vector(6);
end

function eqn = sys_of_eqns(C, params)
    T = params.T;
    P = params.P;
    opt = params.opt;
    Ci0 = params.Ci0;
    tau = params.tau;
    C_struct = get_concentration_struct(C);
    r = get_all_reaction_rates(C_struct, T, P, opt);

    r.ec =  r.r2r - r.r2f - r.r3; 
    r.meoh = (2 * r.r2r) - (2 * r.r2f) - r.r3;
    r.co2 = r.r3;
    r.dmc = r.r2f - r.r2r;
    r.eg = r.r2f - r.r2r; 
    r.me = r.r3; 

    eqn(1) = Ci0.ethylene_carbonate - C(1) + (tau * r.ec);
    eqn(2) = Ci0.methanol - C(3) + (tau * r.meoh);
    eqn(3) = Ci0.carbon_dioxide - C(4) + (tau * r.co2);
    eqn(4) = (-C(5)) + (tau * r.dmc); 
    eqn(5) = (-C(2)) + (tau * r.eg);
    eqn(6) = (-C(6)) + (tau * r.me);
end

function r = get_all_reaction_rates(C, T, P, opt)
    r.r2f = get_reaction_rate(C, '2f', T, P, opt);
    r.r2r = get_reaction_rate(C, '2r', T, P, opt);
    r.r3 = get_reaction_rate(C, '3', T, P, opt);
end
% function F = sys_of_eqns(C, k, b, tau)
%     F(1) = -C(1) - tau * k.k2f * C(1)^0.8 - tau * k.k3 * C(1) + tau * k.k2r * C(5) * C(2) - b(1); 
%     F(2) = -tau * k.k2f * C(1)^0.8 - tau * k.k3 * C(1) - C(3) + tau * k.k2r * C(5) * C(2) - b(2);
%     F(3) = -C(4) + tau * k.k3 * C(1) - b(3);
%     F(4) = -C(5) + tau * k.k2f * C(1)^0.8 - tau * k.k2r * C(5) * C(2) - b(4);
%     F(5) = -C(2) + tau * k.k2f * C(1)^0.8 - tau * k.k2r * C(5) * C(2) - b(5); 
%     F(6) = -C(6) + tau * k.k3 * C(1) - b(6);
% end

function r = get_reaction_rate(C, reaction, T, P, opt)
    % input:
    %   opt = 'isothermal' or 'isobaric'

    k = get_rate_constant(reaction, T, P, opt);
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

function C = get_concentrations(F, T, P, opt, tau)
    V = get_reactor_volume(F, T, P, opt, tau);
    fieldNames = fieldnames(F);
    for i = 1:length(fieldNames)
        C.(fieldNames{i}) = F.(fieldNames{i}).mol * tau / V;
    end
end

function V = get_reactor_volume(F, T, P, opt, tau)
    q_tot = get_total_volumetric_flowrate(F, T, P, opt);
    V = q_tot * tau;
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

    rho = get_constants().densities; 
    rho.methanol = get_methanol_density(T, P);
    rho.carbon_dioxide = get_supercritical_c02_density(T, P, opt);
end

function q_total = get_total_volumetric_flowrate(F, T, P, opt)
    q = get_volumetric_flowrates(F, T, P, opt);
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
	q_total = q_total.value;
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

function k = get_rate_constant(reaction, T, P, opt)
	if strcmp(opt, 'isothermal')
    	k = get_isothermal_rate_constant(reaction, P);
	elseif strcmp(opt, 'isobaric')
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

function rho = get_supercritical_c02_density(T, P, opt)
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
            if withinPressureRange(P)
                rho = 1.6746 * P - 12.592;
                    % NIST / Excel Regression
            else
                rho = NaN;
                disp("get_supercritical_c02_density : ERROR : P out of range")
            end
        case 'isobaric'
            if withinTempRange(T)
                if P < 125 % [C]
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
