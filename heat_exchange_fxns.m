
function Cp_avg = avg_heat_capacity(F)
	global HEAT_CAPACITY_HYDROGEN HEAT_CAPACITY_METHANE HEAT_CAPACITY_ETHANE ...
		HEAT_CAPACITY_ETHYLENE HEAT_CAPACITY_PROPANE HEAT_CAPACITY_BUTANE ...
		HEAT_CAPACITY_WATER

	% weighted average of Cp's 
	F_tot = total_molar_flowrate(F);
	
	Cp_avg = (molar_flowrate(F, 'hydrogen') / F_tot) * HEAT_CAPACITY_HYDROGEN + ...
		(molar_flowrate(F, 'methane') / F_tot) * HEAT_CAPACITY_METHANE + ...
		(molar_flowrate(F, 'ethane') / F_tot) * HEAT_CAPACITY_ETHANE + ...
		(molar_flowrate(F, 'ethylene') / F_tot) * HEAT_CAPACITY_ETHYLENE + ...
		(molar_flowrate(F, 'propane') / F_tot) * HEAT_CAPACITY_PROPANE + ...
		(molar_flowrate(F, 'butane') / F_tot) * HEAT_CAPACITY_BUTANE + ...
		(molar_flowrate(F, 'water') / F_tot) * HEAT_CAPACITY_WATER;

end

function sep = hex(sep, T_out)
	% Temperatures must be in kelvin
	global GJ_PER_KJ

	% GJ/yr 	=  (mol / yr) 				* (kJ / mol K )				  * ( K    -  K   ) * (GJ / kJ)
	sep.heat = total_molar_flowrate(sep.F) * avg_heat_capacity(sep.F) * (T_out - sep.T) * GJ_PER_KJ;
	sep.T = T_out;	
end

function sep = hex_e101(sep)
	global GJ_PER_KJ

	% user inputs
	T_out = 25 + 273.15;	% [ K ]	

	% GJ/yr 	=  (mol / yr) 				* (kJ / mol K )				  * ( K    -  K   ) * (GJ / kJ)
	sep.heat = total_molar_flowrate(sep.F) * avg_heat_capacity(sep.F) * (T_out - sep.T) * GJ_PER_KJ;
	sep.T = T_out;	
end

function [combusted_fuel_flowrates, heatflux_left] = fuel_combustion(heat_flux, flowrates)
	global HYDROGEN METHANE ETHYLENE PROPANE BUTANE;
	global ENTHALPY_METHANE ENTHALPY_PROPANE ENTHALPY_BUTANE HEAT_CAPACITY_ETHANE;
	global MT_PER_KT G_PER_KT GJ_PER_KJ KJ_PER_GJ MOLMASS_METHANE KT_PER_G MOLMASS_BUTANE ...
			MOLMASS_PROPANE PSA_TOGGLE ENTHALPY_HYDROGEN MOLMASS_HYDROGEN

	% Note! : Longest Chain Hydrocarbons are cheapest to combust

	% initialize all values in the array to be zero 
	combusted_fuel_flowrates = flowrates * 0;

	% LOGIC : Goes through each heat source in order, returns if the heat flux supplied is sufficient.
	heatflux_left = heat_flux; 

	% (GJ / yr)           = (kt / yr)          * (g / kt) * (kJ / g)        * (GJ / kJ)
	Q_combust_all_hydrogen = flowrates(HYDROGEN) * G_PER_KT * ENTHALPY_HYDROGEN * GJ_PER_KJ;

	if (~PSA_TOGGLE)
		% Hydrogen
		if (heatflux_left > Q_combust_all_hydrogen)
			combusted_fuel_flowrates(HYDROGEN) = flowrates(HYDROGEN);
			heatflux_left = heatflux_left - Q_combust_all_hydrogen;
		else
			% (kt / yr)                       = ((GJ)                 ) * (KJ / GJ) *
			combusted_fuel_flowrates(HYDROGEN) = (heatflux_left) * KJ_PER_GJ * ...
				... % (mol / KJ)        * (g / mol)       * (kt / g)
				( 1 / ENTHALPY_HYDROGEN) * MOLMASS_HYDROGEN * KT_PER_G;
			heatflux_left = 0;
			return
		end
	end

	% (GJ / yr) 		  = (kt / yr)          * (g / kt) * (kJ / g)		* (GJ / kJ)
	Q_combust_all_methane = flowrates(METHANE) * G_PER_KT * ENTHALPY_METHANE * GJ_PER_KJ;
	
	% Methane
	if (heatflux_left > Q_combust_all_methane)
		combusted_fuel_flowrates(METHANE) = flowrates(METHANE);
		heatflux_left = heatflux_left - Q_combust_all_methane;
	else
		% (kt / yr) 					  = ((GJ)                 ) * (KJ / GJ) *
		combusted_fuel_flowrates(METHANE) = (heatflux_left) * KJ_PER_GJ * ...
			... % (mol / KJ) 		* (g / mol) 	  * (kt / g)
			( 1 / ENTHALPY_METHANE) * MOLMASS_METHANE * KT_PER_G;
		heatflux_left = 0;
		return
	end
end
