clc; clear; close all ;

test_level2_feedstream()
% test_constants()

function void = test_units(units)
	disp("test units")
	units

	units.mass
	units.energy
	units.temperature
	units.value
	units.pressure
	units.time
	units.volume
	units.heat
end

function void = test_constants()
	
	disp("const")
	const = get_constants()

	test_units(const.units);
	% test_molar_mass(const.molar_mass);
	% test_stoich(const.stoich);
	% test_thermo(const.thermo);
	% test_econ(const.econ);

end

function void = test_level2_feedstream()

    const = get_constants(); 

    user = get_user_inputs();

    F_fxns = flowrate_fxns();

	console = get_console_constants();

    F = user.level2.feed_stream;

	fprintf("LEVEL 2 FEEDSTREAM FLOWRATES%s\n",  console.divider);
	fprintf("Carbon Dioxide \n\t%4.3f kta\t %4.3f mol /s\n",F.carbon_dioxide.kta, F.carbon_dioxide.mol)
	fprintf("Ethylene Oxide \n\t%4.3f kta\t %4.3f mol /s\n", F.ethylene_oxide.kta, F.ethylene_oxide.mol);
	fprintf("Methanol \n\t%4.3f kta\t %4.3f mol /s\n", F.methanol.kta, F.methanol.mol);
	fprintf("Ethylene Carbonate \n\t%4.3f kta\t %4.3f mol /s\n", F.ethylene_carbonate.kta, F.ethylene_carbonate.mol);
	fprintf("Ethylene Glycol \n\t%4.3f kta\t %4.3f mol /s\n", F.ethylene_glycol.kta, F.ethylene_glycol.mol);
	fprintf("Methoxy Ethanol \n\t%4.3f kta\t %4.3f mol /s\n", F.methoxy_ethanol.kta, F.methoxy_ethanol.mol);
	
end