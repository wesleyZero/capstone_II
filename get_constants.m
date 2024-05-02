
function const = get_constants()

    const.units = get_unit_conversions();

    const.molar_mass = get_molar_masses();
    
    const.stoich = get_stoichiometric_coeff();
    
    const.thermo = get_thermodynamic_constants();


end 

function thermo = get_thermodynamic_constants()

    thermo.R = 8.314;					% [ J / mol K ]

    % thermo.heat_cap = get_heat_capacities();

    % themrmo.rate_const.isoP.units = "";
    % thermo.rate_const.isoP.k2f = @(T_K) 6.59 * 10^2 * exp( -37000 / thermo.R / T_K);
    % thermo.rate_const.isoP.k2r = @(T_K) 1.19 * 10^4 * exp( -53700 / thermo.R / T_K);
    % thermo.rate_const.isoP.k3 = @(T_K) 1.89 * 10^6 * exp( -82400 / thermo.R / T_K);
    thermo.rate_const = get_rate_constants();

    % thermo.rate.units = "";
    % thermo.rate.r2f = @(T_K,C_EC) thermo.rate_const.k2f(T_K) * C_EC^0.8;
    % thermo.rate.r2r = @(T_K,C_DMC, C_EG) thermo.rate_const.k2r(T_K) * C_DMC * C_EG;
    % thermo.rate.r3 = @(T_K, C_EC) thermo.rate_const.k3(T_K) * C_EC;

    thermo.rate = get_reaction_rates();

    % thermo.rate_const.isoT.unitt = "";
    % thermo.rate_const.isoT.k2f_lowRho = 0.013; 
    % thermo.rate_const.isoT.k2f_highRho = @(rho) 0.02486 - 4.943 * 10^-5 * rho;
    % thermo.rate_const.isoT.k3 = @(rho) 3.014 * 10^-4 * exp(-5.99 * rho);

    thermo.rate_const = get_isothermal_rate_constants();

    thermo.enthalpy.heat_formation.units = "BS VALUES";
    thermo.enthalpy.heat_formation.ethylene_oxide = 44.0526;            
        % source:  
    thermo.enthalpy.heat_formation.carbon_dioxide =  44.0095;          
        % source:  
    thermo.enthalpy.heat_formation.ethylene_carbonate =  88.0621;      
        % source:  
    thermo.enthalpy.heat_formation.methanol = 32.0419;                 
        % source: 
    thermo.enthalpy.heat_formation.dimethyl_carbonate = 90.0779;       
        % source: 
    thermo.enthalpy.heat_formation.ethylene_glycol =62.0678;          
        % source: 
    thermo.enthalpy.heat_formation.methoxy_ethanol = 76.10;           
        % source: 

    stoich = get_stoichiometric_coeff()
    thermo.enthalpy.e1 = 

end 

function rate_const = get_rate_constants() 

    rate_const.isoP.units = "";
    rate_const.isoP.k2f = @(T_K) 6.59 * 10^2 * exp( -37000 / thermo.R / T_K);
    rate_const.isoP.k2r = @(T_K) 1.19 * 10^4 * exp( -53700 / thermo.R / T_K);
    rate_const.isoP.k3 = @(T_K) 1.89 * 10^6 * exp( -82400 / thermo.R / T_K);
end

function rate = get_reaction_rates()

    rate.units = "";
    rate.r2f = @(T_K,C_EC) thermo.rate_const.k2f(T_K) * C_EC^0.8;
    rate.r2r = @(T_K,C_DMC, C_EG) thermo.rate_const.k2r(T_K) * C_DMC * C_EG;
    rate.r3 = @(T_K, C_EC) thermo.rate_const.k3(T_K) * C_EC;
end


function isoT = get_isothermal_rate_constants() 

    isoT.unitt = "";
    isoT.k2f_lowRho = 0.013; 
    isoT.k2f_highRho = @(rho) 0.02486 - 4.943 * 10^-5 * rho;
    isoT.k3 = @(rho) 3.014 * 10^-4 * exp(-5.99 * rho);
end

% function heat_cap = get_heat_capacities()

%     heat_cap.units = 'kJ / mol K';

%     heat_cap.ethylene_oxide =             % [ kJ / mol K ]
%         % source:
%     heat_cap.carbon_dioxide =             % [ kJ / mol K ]
%         % source:
%     heat_cap.ethylene_carbonate =        % [ kJ / mol K ]
%         % source:
%     heat_cap.methanol =                 % [ kJ / mol K ]
%         % source:
%     heat_cap.dimethyl_carbonate =         % [ kJ / mol K ]
%         % source:
%     heat_cap.ethylene_glycol =        % [ kJ / mol K ]
%         % source:
%     heat_cap.methoxy_ethanol = ;              % [ kJ / mol K ]
%         % source:

% end

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
	units.temperature.c_to_k = @(T_K) T_K - 273.15;        % [ C -> K ]

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
