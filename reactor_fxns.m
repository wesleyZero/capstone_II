% Reactor Simulation Functions 

function fxns = reactor_fxns() 
    % fxns.get_reaction_rate = @get_reaction_rate;
    fxns.get_reactor_flows = @get_reactor_flows;
    fxns.get_conversion = @get_conversion;
end

function [F_fresh, F_real_feed, F_real_eff, R, V_rxtr] = get_reactor_flows(F_real_feed_basis, T, P, opt, tau)
    F_fresh = NaN; F_real_feed = NaN; F_real_eff = NaN; R = NaN;
    % basis calculations for the real reactor 
    q_tot.basis = get_total_volumetric_flowrate(F_real_feed_basis, T, P, opt);
    V_rxtr.basis = q_tot.basis * tau;
    C_out = get_reactor_effluent_concentrations(F_real_feed_basis, T, P, opt, tau);

    if any(imag(C_out) ~= 0)
        % disp('ERROR : Complex valued concentrations');
        return 
    elseif any(real(C_out) < 0)
        % disp('ERROR : Negative valued concentrations');
        return
    else
        % disp('Valid solution!!!!!!!!!!')
        C_out = get_concentration_struct(C_out);
    end

    F_real_eff_basis = conc_to_flowrate(C_out, q_tot.basis);
        % assumption : liquid flow has no change in vol in effluent 

    % Plant Scale Calculations 
    [F_fresh, F_real_feed, F_real_eff, R] = get_plant_flowrates(F_real_feed_basis, F_real_eff_basis);
    scale_factor = get_scale_factor(F_real_eff_basis);
    V_rxtr.plant = V_rxtr.basis * scale_factor ; 

    V_rxtr = V_rxtr.plant;

end

function chi = get_conversion(F_real_feed, F_real_eff)
    if ~isstruct(F_real_feed) || ~isstruct(F_real_eff), chi = NaN;, return, end;
    if ~isreal(F_real_feed.ethylene_carbonate.mol) || ~isreal(F_real_eff.ethylene_carbonate.mol), chi = NaN;, end
    chi = ( F_real_feed.ethylene_carbonate.mol - F_real_eff.ethylene_carbonate.mol)  / ...
            F_real_feed.ethylene_carbonate.mol;
end


function [F_fresh, F_real_feed, F_real_effluent, R] = get_plant_flowrates(F_real_feed_basis, F_real_eff_basis)
    scale_factor = get_scale_factor(F_real_eff_basis);

    F_real_effluent = get_scaled_flowrate(F_real_eff_basis, scale_factor);
    F_real_feed = get_scaled_flowrate(F_real_feed_basis, scale_factor);

    F_virt_feed = get_virtual_reactor_feed(F_real_feed);
    [F_fresh, R] = get_recycle_and_fresh_flowrates(F_virt_feed, F_real_effluent);

end

function F_scaled  = get_scaled_flowrate(F_basis, scale_factor)
    % Conserves mass by scaling the mass flowrates 

    flow_fxns = flowrate_fxns();
    F_scaled = flow_fxns.get_blank_flowstream();
    fieldNames = fieldnames(F_basis);
    for i = 1:length(fieldNames)
        if strcmp(fieldNames{i},'units') , continue, end ;
        F_scaled.(fieldNames{i}).kta = scale_factor * F_basis.(fieldNames{i}).kta;
	end
    F_scaled = flow_fxns.set_F_kta(F_scaled); 
end

function F_virt_feed = get_virtual_reactor_feed(F_virtual_effluent)
    user = get_user_inputs();
    flow_fxns = flowrate_fxns();
    F_virt_feed = flow_fxns.get_blank_flowstream();

    % Assume that C02 is in excess
    % Get the flow into the virtual reactor
    F_virt_feed.ethylene_oxide.mol = F_virtual_effluent.ethylene_carbonate.mol;
        % Assume that EO -> EC Complete conversion in virtual reactor
    F_virt_feed.methanol.mol = F_virt_feed.ethylene_oxide.mol * user.level3.molar_ratio_methanol_EO;
    F_virt_feed.carbon_dioxide.mol = F_virt_feed.ethylene_oxide.mol * user.level3.molar_ratio_carbon_dioxide_EO;
    F_virt_feed = flow_fxns.set_F_mol(F_virt_feed);
end

function [F_fresh, R] = get_recycle_and_fresh_flowrates(F_virt_feed, F_real_effluent)
    flow_fxns = flowrate_fxns();

    % Initialize 
    R = flow_fxns.get_blank_flowstream();
    F_fresh = flow_fxns.get_blank_flowstream();

    % Recycle flow
    R.ethylene_carbonate.mol = F_real_effluent.ethylene_carbonate.mol;
    R.methanol.mol = F_real_effluent.methanol.mol;
    R.carbon_dioxide.mol = F_real_effluent.carbon_dioxide.mol;
    R = flow_fxns.set_F_mol(R);

    % Fresh feed flow 
    F_fresh.ethylene_oxide.mol = F_virt_feed.ethylene_oxide.mol - R.ethylene_carbonate.mol;
    F_fresh.methanol.mol = F_virt_feed.methanol.mol - R.methanol.mol;
    F_fresh.carbon_dioxide.mol = F_virt_feed.carbon_dioxide.mol - R.carbon_dioxide.mol;
    F_fresh = flow_fxns.set_F_mol(F_fresh);

    % F_fresh.ethylene_oxide.mol = F_fresh.ethylene_carbonate.mol;
    % F_fresh.ethylene_carbonate.mol = 0;
		% EC should be turned back into EO because we need the feed into the virtual reactor
		% ?? Look into more detail of the EO / EC and recycle stream because the
		% VR really complicates things


end


function scale_factor = get_scale_factor(F_out)
    scale_factor = (get_user_inputs().dmc_production_rate) / F_out.dimethyl_carbonate.kta; 

end

function F = conc_to_flowrate(C, q_tot)
    % ?? Assumption : densities don't change after reaction (which will change
    % the T and P, thus densities should change) modify q_tot(T, P) to get a more
    % accurate result
    flow_fx = flowrate_fxns();
    F = flow_fx.get_blank_flowstream();

    fieldNames = fieldnames(C);
    for i = 1:length(fieldNames)
        if strcmp(fieldNames{i},'units') , continue, end ;

        % mol / s         = ( L / s )         * (mol / L)
        F.(fieldNames{i}).mol = q_tot * C.(fieldNames{i}); 
    end

    F = flow_fx.set_F_mol(F); 
end

function C_vector = get_concentration_vector(F, T, P, opt, tau)
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
    tau = tau / 10;
    user = get_user_inputs();
    Ci0 = get_concentrations(F, T, P, opt, tau);
    Ci_init = get_concentration_vector(F, T, P, opt, tau);
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

    % % Carbon dioxide 
    % % (L / s)        = (mol / s)            * (g / mol)                      
    % q.carbon_dioxide = F.carbon_dioxide.mol * const.molar_mass.carbon_dioxide * ... 
    %     ...% (kg / g)        * (m^3 / kg)             * (L / m^3) 
    %     const.units.mass.kg_per_g * (1/rho.carbon_dioxide) * const.units.volume.l_per_m3;
    % % Methanol
    % q.methanol = F.methanol.mol * const.molar_mass.methanol * ...
    %     const.units.mass.kg_per_g * (1 / rho.methanol) * const.units.volume.l_per_m3;
    % % Ethylene Carbonate
    % q.ethylene_carbonate = F.ethylene_carbonate.mol * const.molar_mass.ethylene_carbonate *...
    %      const.units.mass.kg_per_g * (1 / rho.ethylene_carbonate) * const.units.volume.l_per_m3;


    fieldNames = fieldnames(F);
    for i = 1:length(fieldNames)
        if strcmp(fieldNames{i},'units')
            % rho.(fieldNames{i}) = 'mol / L';
            continue
        end
        % (L / s)        = (mol / s)            * (g / mol)                      
        q.(fieldNames{i}) = F.(fieldNames{i}).mol * const.molar_mass.(fieldNames{i}) * ... 
            ...% (kg / g)             * (m^3 / kg)             * (L / m^3) 
            const.units.mass.kg_per_g * (1/rho.(fieldNames{i})) * const.units.volume.l_per_m3;
    end
    q.units = 'L / s' ;
    
end

function k = get_rate_constant(reaction, T, P, opt)
	if strcmp(opt, 'isothermal')
    	k = get_isothermal_rate_constant(reaction, T, P);
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

function k = get_isothermal_rate_constant(reaction, T, P)
    % input:
    %   P [ bar ]
    % These functions are from the research paper
    rho = get_supercritical_c02_density(T, P, 'isothermal');
    switch reaction
        case '2f'
            if rho > 246.82 % [g / L]
                k = (2.486 * 10^(-2)) - (4.943 * (10^(-5)) * rho);
            else
                k = (1.362 * 10^(-2)) - (1.569 * (10^(-6)) * rho);
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
    %   Isobaric model is at 150 bar
    %   Isothermal model is at 140 C
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
                if P < 125 % [ Bar ]
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
