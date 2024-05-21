
function const = get_constants()
    const.units = get_unit_conversions();
    const.molar_mass = get_molar_masses();
    const.stoich = get_stoichiometric_coeff();
    const.thermo = get_thermodynamic_constants();
    const.econ = get_economic_constants();
    const.densities = get_species_densities();
    % const.heat_cap = get_heat_capacities();

end 


function rho = get_species_densities();
    rho.units = "kg / m^3 "
    rho.ethylene_carbonate = 1.3214;
        % https://pubchem.ncbi.nlm.nih.gov/compound/Ethylene-carbonate#section=Melting-Point
        % g / ml = kg / m^3 right? 

end 


function econ = get_economic_constants()

    econ.value.units = "$ / MT";
    econ.value.dimethyl_carbonate = 1100;
    econ.value.ethylene_glycol = 500;
    econ.value.methanol = 600;
    econ.value.ethylene_oxide = 1250;
    econ.value.carbon_dioxide_feedstock = 45;

    econ.value.fuel.units = "$ / GJ";
    econ.value.fuel.natural_gas = 3;

    % econ.value
end

function thermo = get_thermodynamic_constants()
    thermo.R = 8.314;					% [ J / mol K ]
    % thermo.heat_cap = get_heat_capacities();
    % thermo.rate_const = get_rate_constants();
    % thermo.rate = get_reaction_rates();
    % thermo.rate_const = get_isothermal_rate_constants();
    thermo.enthalpy = get_reaction_enthalpies();
end 

function enthalpy = get_reaction_enthalpies()
    enthalpy.e1.kj_per_mol = -60.8;
    enthalpy.e2.kj_per_mol = -53.5;
    enthalpy.e3.kj_per_mol = -55.3;

end

function stoich = get_stoichiometric_coeff()
    % WARNING : DOES NOT GIVE NEGATIVE VALUES

    stoich.r1.ethylene_oxide = 1;
    stoich.r1.carbon_dioxide = 1; 
    stoich.r1.ethylene_carbonate = 1;

    stoich.r2.ethylene_carbonate = 1;
    stoich.r2.methanol = 2; 
    stoich.r2.dimethyl_carbonate = 1;
    stoich.r2.ethylene_glycol = 1;

    stoich.r3.ethylene_carbonate = 1;
    stoich.r3.methanol = 1;
    stoich.r3.methoxy_ethanol = 1; 
    stoic.r3.carbon_dioxide = 1;
end

function units = get_unit_conversions() 
    
    % Mass
    units.mass.mt_per_kt = 10^3;     % [ MT / kt ]
    units.mass.g_per_kt = 10^9;       % [ g / kt ]
    units.mass.kt_per_g = 10^-9;       % [ kt / g ]  
    units.mass.kg_per_kt =  10^6;     % [ kg / MT ]
    units.mass.mt_per_g = 10^-6;       % [ MT / g ] 

    % Energy
    units.energy.gj_per_kj = 10^-6;        % [ GJ / kJ ]
    units.energy.kj_per_gj = 10^6;        % [ kJ / GJ ]

    % Temperature      
	units.temperature.c_to_k = @(T_C) T_C + 273.15;        % [ C -> K ]
	units.temperature.k_to_c = @(T_K) T_K - 273.15;        % [ K -> C ]

    % Value
    units.value.mmdolla_per_dolla = 10^-6;    % [ $ MM / $]
    units.value.dolla_per_mmdolla = 10^6;    % [ $ / $ MM ]

    % Pressure
    units.pressure.bar_per_psia = 0.0689476;    % [ Bar / Psia ]

    % Time
    units.time.yr_per_sec = 1 / (3.154 * 10^7);    % [ yr / s ]
    units.time.sec_per_yr = 3.154 * 10^7;       % [ s / yr ]
    units.time.yr_per_hr = (1/8760 );       % [ yr / hr ]
    units.time.hr_per_yr = 8760;       % [ hr / yr ]

    % Volumes 
    units.volume.m3_per_l = 0.001;

    % heat 
    units.heat.millionbtu_per_gj = 1.0551;       % [ ]
end

function molar_mass = get_molar_masses()
    molar_mass.units = "g / mol";

    molar_mass.ethylene_oxide = 44.0526;            % [ g / mol ]
        % source: https://webbook.nist.gov/cgi/cbook.cgi?ID=C75218&Mask=80 
    molar_mass.carbon_dioxide =  44.0095;           % [ g / mol ]
        % source: https://webbook.nist.gov/cgi/cbook.cgi?ID=124-38-9 
    molar_mass.ethylene_carbonate =  88.0621;       % [ g / mol ]
        % source: https://webbook.nist.gov/cgi/cbook.cgi?ID=C96491&Mask=200 
    molar_mass.methanol = 32.0419;                  % [ g / mol ]
        % source: https://webbook.nist.gov/cgi/cbook.cgi?ID=67-56-1
    molar_mass.dimethyl_carbonate = 90.0779;        % [ g / mol ]
        % source: https://webbook.nist.gov/cgi/cbook.cgi?ID=C616386&Mask=200
    molar_mass.ethylene_glycol =62.0678;            % [ g / mol ]
        % source: https://webbook.nist.gov/cgi/cbook.cgi?ID=107-21-1 
    molar_mass.methoxy_ethanol = 76.10;              % [ g / mol ]
        % source: https://webbook.nist.gov/cgi/cbook.cgi?ID=109-86-4
end


