clc; clear; close all ;

test_constants()

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

