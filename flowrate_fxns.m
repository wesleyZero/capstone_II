
function fxns = flowrate_fxns()
    fxns.set_F_kta = @set_F_kta;
    fxns.set_F_mol = @set_F_mol;
    fxns.get_blank_flowstream = @get_blank_flowstream;
    fxns.get_basis_feed_flowrates = @get_basis_feed_flowrates;
end

function F = get_basis_feed_flowrates()
    user = get_user_inputs();

    F = get_blank_flowstream();

    % Define the basis ratios
    F.ethylene_carbonate.mol = 1; 
    F.methanol.mol = F.ethylene_carbonate.mol * user.level3.molar_ratio_methanol_EO;
    F.carbon_dioxide.mol = F.ethylene_carbonate.mol * user.level3.molar_ratio_carbon_dioxide_EO;

    % set the flowrates
    F = set_F_mol(F);

    % update the mol fractions
    F = set_mol_fractions(F);
end 

function F = get_blank_flowstream()
    F.carbon_dioxide.mol = 0;
    F.ethylene_oxide.mol = 0; 
    F.methanol.mol = 0; 
    F.ethylene_carbonate.mol = 0; 
    F.ethylene_glycol.mol = 0; 
    F.methoxy_ethanol.mol = 0; 
	F.dimethyl_carbonate.mol = 0; 

    F = set_F_mol(F);

    % update the mol fractions
    F = set_mol_fractions(F);
end

function F = set_F_kta(F) 
    % Calculates the molar flow rates in mol / s for a given flow rates in kta
    const = get_constants();

    % Calculate the mol / s flowrate for given kta flowrate inputs
    fieldNames = fieldnames(F);
    for i = 1:length(fieldNames)
        % mol / s            = (kt / yr )           * (g / kt) 
        F.(fieldNames{i}).mol = F.(fieldNames{i}).kta * const.units.mass.g_per_kt * ...
            ... %  (yr / sec)      * (mol / g)
            const.units.time.yr_per_sec * (1 / const.molar_mass.(fieldNames{i}));
    end

    % update the mol fractions with new calculated molar flowrates
    F = set_mol_fractions(F);
end

function F = set_F_mol(F)
    % Calculates the the molar flow rates in kta for given flow rates in mol / s
    const = get_constants();

    % Set the kta flowrates, for given flowrates in mol / s 
    fieldNames = fieldnames(F);
    for i = 1:length(fieldNames)
        % kt / yr            = (mol / s)            * (g / mol)
        F.(fieldNames{i}).kta = F.(fieldNames{i}).mol * const.molar_mass.(fieldNames{i}) * ...
            ... % (kt / g)            * (s / yr) 
            const.units.mass.kt_per_g * const.units.time.sec_per_yr;
    end

   % update the mol fractions
    F = set_mol_fractions(F);
end

function F = set_mol_fractions(F)

    % Find total molar flow rate
    fieldNames = fieldnames(F);
    F_total = 0;
    for i = 1:length(fieldNames)
		species = fieldNames{i};
        F_total = F_total + F.(species).mol;
    end

    % Calculate the mol fractions
    for i = 1:length(fieldNames)
        if F_total
            F.(fieldNames{i}).x = F.(fieldNames{i}).mol / F_total;
        else % For blank flowstreams
            F.(fieldNames{i}).x = 0;
        end
    end

end