clc; clear; close all ;

test_level2_feedstream()
test_units()
% test_constants()

function void = test_units()
	const = get_constants();
	console = get_console_constants();

	fprintf("TESTING UNIT CONVERSIONS%s\n",  console.divider);
	fprintf("\tmass %s\n",  console.divider);
	const.units.mass
	
	fprintf("\tenergy %s\n",  console.divider);
	const.units.energy
	
	fprintf("\ttemperature %s\n",  console.divider);
	const.units.temperature
	
	fprintf("\tvalue %s\n",  console.divider);
	const.units.value
	
	fprintf("\tpressure %s\n",  console.divider);
	const.units.pressure
	
	fprintf("\ttime %s\n",  console.divider);
	const.units.time
	
	fprintf("\tvolume %s\n",  console.divider);
	const.units.volume
	
	fprintf("\theat %s\n",  console.divider);
	const.units.heat
	
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