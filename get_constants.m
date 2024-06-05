
function const = get_constants()
    const.units = get_unit_conversions();
    const.molar_mass = get_molar_masses();
    const.stoich = get_stoichiometric_coeff();
    const.thermo = get_thermodynamic_constants();
    const.econ = get_economic_constants();
    const.densities = get_species_densities();
    % const.heat_cap = get_heat_capacities();
    const.pi = 3.14159;

end 


function rho = get_species_densities();
    rho.units = "kg / m^3 ";

    rho.ethylene_carbonate = 1.3214; % g / ml 
        % https://pubchem.ncbi.nlm.nih.gov/compound/Ethylene-carbonate#section=Melting-Point
        % kg / m3 = (g/ml) * (1000ml / L) * (kg / 1000 g) * (1000L / m^3)
    rho.ethylene_carbonate = rho.ethylene_carbonate * 1000;

    rho.dimethyl_carbonate = 1.069; % g / ml
        % https://www.sigmaaldrich.com/US/en/product/sial/517127
    rho.dimethyl_carbonate = rho.dimethyl_carbonate * 1000;

    rho.ethylene_glycol = 1.113; % g / ml  
        % https://www.sigmaaldrich.com/US/en/product/sial/324558
    rho.ethylene_glycol = rho.ethylene_glycol * 1000;

    rho.methoxy_ethanol = 0.965; % g / ml 
        % https://www.sigmaaldrich.com/US/en/product/sial/284467
    rho.methoxy_ethanol = rho.methoxy_ethanol * 1000;

    rho.ethylene_oxide = 0.882; % g / ml 
        % https://www.sigmaaldrich.com/US/en/product/aldrich/387614
    rho.ethylene_oxide = rho.ethylene_oxide * 1000;
    
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
    units.mass.mt_per_kt = 10^3;     
    units.mass.g_per_kt = 10^9;       
    units.mass.kt_per_g = 10^-9;      
    units.mass.kg_per_kt =  10^6;    
    units.mass.mt_per_g = 10^-6;       
    units.mass.kg_per_g = 10^-3;
    units.mass.g_per_kg = 10^3;

    % Energy
    units.energy.gj_per_kj = 10^-6;       
    units.energy.kj_per_gj = 10^6;      

    % Temperature      
	units.temperature.c_to_k = @(T_C) T_C + 273.15;       
	units.temperature.k_to_c = @(T_K) T_K - 273.15;      

    % Value
    units.value.mmdolla_per_dolla = 10^-6;    % [ $ MM / $]
    units.value.dolla_per_mmdolla = 10^6;    % [ $ / $ MM ]

    % Pressure
    units.pressure.bar_per_psia = 0.0689476;   
    units.pressure.psia_per_bar = 1 / units.pressure.bar_per_psia ;

    % Time
    units.time.yr_per_sec = 1 / (3.154 * 10^7);   
    units.time.sec_per_yr = 3.154 * 10^7;      
    units.time.yr_per_hr = (1/8760 );      
    units.time.hr_per_yr = 8760;    

    % Volumes 
    units.volume.m3_per_l = 0.001;
    units.volume.l_per_m3 = 1000;

    % Length
    units.length.ft_per_m =  3.28084;
    units.length.m_per_ft = 1 / units.length.ft_per_m;

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

    molar_mass.aniline = 93.13;
    molar_mass.water = 18;
end


