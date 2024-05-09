

function fxns = flowrate_fxns()

    fxns.set_F_kta = @set_F_kta;
    fxns.set_F_mol = @set_F_mol;

end

function F = set_F_kta(F) 
    % Calculates the molar flow rates in mol / s for a given flow rates in kta

    const = get_constants();

    % mol / s            = (kt / yr )           * (g / kt) 
    F.carbon_dioxide.mol = F.carbon_dioxide.kta * const.units.mass.g_per_kt * ...
        ... %  (yr / sec)      * (mol / g)
        const.units.time.yr_per_sec * (1 / const.molar_mass.carbon_dioxide);

    F.ethylene_oxide.mol = F.ethylene_oxide.kta * const.units.mass.g_per_kt * ...
        const.units.time.yr_per_sec * (1 / const.molar_mass.ethylene_oxide);

    F.methanol.mol = F.methanol.kta * const.units.mass.g_per_kt * ...
        const.units.time.yr_per_sec * (1 / const.molar_mass.methanol);
    
    F.ethylene_carbonate.mol = F.ethylene_carbonate.kta * const.units.mass.g_per_kt * ...
        const.units.time.yr_per_sec * (1 / const.molar_mass.ethylene_carbonate);
    
    F.ethylene_glycol.mol = F.ethylene_glycol.kta * const.units.mass.g_per_kt * ...
        const.units.time.yr_per_sec * (1 / const.molar_mass.ethylene_glycol);
    
    F.methoxy_ethanol.mol = F.methoxy_ethanol.kta * const.units.mass.g_per_kt * ...
        const.units.time.yr_per_sec * (1 / const.molar_mass.methoxy_ethanol);


end

function F = set_F_mol(F)
    const = get_constants();

    % kt / yr            = (mol / s)            * (g / mol)
    F.carbon_dioxide.kta = F.carbon_dioxide.mol * const.molar_mass.carbon_dioxide * ...
        ... % (kt / g)            * (s / yr) 
        const.units.mass.kt_per_g * const.units.time.sec_per_yr;

    F.ethylene_oxide.kta = F.ethylene_oxide.mol * const.molar_mass.ethylene_oxide * ...
        const.units.mass.kt_per_g * const.units.time.sec_per_yr;

    F.methanol.kta = F.methanol.mol * const.molar_mass.methanol * ...
        const.units.mass.kt_per_g * const.units.time.sec_per_yr;

    F.ethylene_carbonate.kta = F.ethylene_carbonate.mol * const.molar_mass.ethylene_carbonate * ...
        const.units.mass.kt_per_g * const.units.time.sec_per_yr;

    F.ethylene_glycol.kta = F.ethylene_glycol.mol * const.molar_mass.ethylene_glycol * ...
        const.units.mass.kt_per_g * const.units.time.sec_per_yr;

    F.methoxy_ethanol.kta = F.methoxy_ethanol.mol * const.molar_mass.methoxy_ethanol * ...
        const.units.mass.kt_per_g * const.units.time.sec_per_yr;

end
