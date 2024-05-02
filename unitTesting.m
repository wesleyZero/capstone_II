clc; clear; close all ;


test_constants()


function void = test_constants()

    const = get_constants()

    const.units.mass
	const.units.energy
	const.units.temperature
	const.units.time
	const.units.pressure
	const.units.volume
	const.units.heat

end