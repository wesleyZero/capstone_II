%% MASTER FUNCTION

function fxns = get_economic_functions()
    fxns.get_npv = @get_npv;
	fxns.get_work_min_npv = @get_work_min_npv;

end

function energy = get_energy_plant()
	
	energy = 0;
end

function Fp = get_closest_Fp(P)
    % Define the table values
    pressures = [50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000];
    Fp_values = [1.00, 1.05, 1.15, 1.20, 1.35, 1.45, 1.60, 1.80, 1.90, 2.30, 2.50];
    
    % Find the closest pressure in the table
    [~, idx] = min(abs(pressures - P));
    
    % Get the corresponding Fp value
    Fp = Fp_values(idx);
end


function installed_cost = get_cost_reactor(V, P)
	% V = volume of reactor in m^3
	% P = bar
	user = get_user_inputs();
	const = get_constants();
	coeff = (user.marshall_and_swift_index / 280) * 101.9;
	F_m = 1.0;
		% Assumption : Carbon Steel is the material ??
	P = P * const.units.pressure.psia_per_bar;
	F_p = get_closest_Fp(P );
	F_c = F_m * F_p;
	k = 2.942 / 4.413; % D / H ratio is assumed constant 
	D.m = (4 * k * V / const.pi)^(1/3);
	H.m = D.m / k;
	D.ft = D.m * const.units.length.ft_per_m;
	H.ft = H.m * const.units.length.ft_per_m; 
	installed_cost = coeff * (D.ft)^1.066 * H.ft^0.82 * (2.18 + F_c);
	installed_cost = installed_cost * 1.5;
	V;

end

function isbl_capital_cost = get_ISBL_aspen_data(V_rxtr, P, opt)
	if strcmp(opt, 'aspen')
		isbl_capital_cost = get_cost_reactor(V_rxtr, P);
		disp("aspen")
	elseif strcmp(opt, 'matlab')
		isbl_capital_cost = get_cost_reactor(V_rxtr, P);	
		disp("matlab")
	else
		isbl = NaN;
		disp("ERROR | get_ISBL_aspen_data | opt not valid");
	end
end

		

function charge = get_CO2_sustainability_charge()
	charge = 0;
end

function cost = get_separation_system_cost()

	cost = 0;
end

function lifetime_npv = get_work_min_npv(F, T, P, V_rxtr, conversion, opt)
	sep_fxns = separation_fxns();
	w_min = sep_fxns.get_work_min(F, T);

	npv_params.mainProductRevenue = get_main_product_revenue(F);
	npv_params.byProductRevenue = get_byproduct_revenue(F);
	npv_params.rawMaterialsCost = get_raw_material_cost(F);
	energy = get_energy_plant();
	npv_params.utilitiesCost = get_utilities_cost(energy);
	npv_params.conversion = conversion;
	% npv_params.ISBLcapitalCost = get_cost_reactor(V_rxtr, P);
	npv_params.ISBLcapitalCost = get_ISBL_aspen_data(V_rxtr, P, opt);

	npv_params.CO2sustainabilityCharge = get_CO2_sustainability_charge();

	lifetime_npv = get_npv(npv_params) / 10^6 ;
end

function value = chemical_value(F, species)

	switch true
		case strcmp(species, 'dimethyl_carbonate')
			value = 1100 * 10^3 * F.dimethyl_carbonate.kta;
		case strcmp(species, 'ethylene_glycol')
			value = 500 * 10^3 * F.ethylene_glycol.kta;
		case strcmp(species, 'methanol')
			value = 600 * 10^3 * F.methanol.kta;
		case strcmp(species, 'ethylene_oxide')
			value = 1250 * 10^3 * F.ethylene_glycol.kta; 
		case strcmp(species, 'carbon_dioxide')
			value = 45 * 10^3 * F.carbon_dioxide.kta;
		otherwise
			value = NaN;
			fprintf("ERROR | chemical_value | invalid case | %s\n", species);
	end
end

function value = get_main_product_revenue(F)
	% Products are the value of DMC
	value = chemical_value(F, 'dimethyl_carbonate');
end

function value = get_byproduct_revenue(F)
	% byproducts are EG
	value = chemical_value(F, 'ethylene_glycol');

end

function cost = get_raw_material_cost(F)
	% raw material costs are EO, MeOH, C02, water
	% ask TJ how to cost the water for the sep unit

	cost = 0; 
	cost = cost + chemical_value(F, 'ethylene_oxide');
	cost = cost + chemical_value(F, 'methanol');
	cost = cost + chemical_value(F, 'carbon_dioxide');
	% cost = cost + 
		% water ??
end

function cost = get_utilities_cost(energy)
	% utilities are electricity, fuel??
	rate.dollar_per_GJ = 3;
	cost = rate.dollar_per_GJ * energy;
end

function value = get_CO2_sustainability_value()
	% 45 / MT, confirm with TJ

	value = 0;

end

% function cost = cost_reactor(V_plant_input)
% 	global FT_PER_METER STEAM_TO_FEED_RATIO
% 	FT_PER_METER = 3.28084;
% 	% ??? WHAT ARE THE UNITS OF TIME
	
% 	pi = 3.14159;
% 	D = 0.05;								% [m]
% 	V_plant_max = pi * (0.025)^2 * 20; 		%[m^3]
	
% 	% Reactors have a max length, so calculate the number of full size reactors
% 	% and add it to the cost of the one non-max length reactor

% 	cost = 0;

% 	% Find the Cost of the max-sized reactors
% 	num_of_additional_reactors = int64(V_plant_input / V_plant_max);	
% 	num_of_additional_reactors = double(num_of_additional_reactors);
	
% 	V_plant = V_plant_max;
% 	factor_1 = 4.18;
% 	factor_2 = (V_plant / (pi * (D/2)^2) * FT_PER_METER)^0.82;
% 	factor_3 = (101.9 * D * FT_PER_METER)^1.066;
% 	factor_4 = (1800 / 280);
% 	cost_max_reactor = factor_1 * factor_2 * factor_3 * factor_4; 
% 	cost = cost + num_of_additional_reactors * cost_max_reactor;
	
% 	% Find the cost of the non-max size reactor 
% 	V_plant = V_plant_input - V_plant_max * num_of_additional_reactors;
% 	if V_plant < 0
% 		V_plant = 0;
% 	end
% 	factor_1 = 4.18;
% 	factor_2 = (V_plant / (pi * (D/2)^2) * FT_PER_METER)^0.82;
% 	factor_3 = (101.9 * D * FT_PER_METER)^1.066;
% 	factor_4 = (1800 / 280);
% 	cost = cost + factor_1 * factor_2 * factor_3 * factor_4;

% end


function lifetime_npv = get_npv(npv)
	% USER_INPUTS | All inputs are in units of $MM
		% npv.mainProductRevenue = value_ethylene(P_ethylene);
		% npv.byProductRevenue = value_h2_chem(P_hydrogen - combusted_hydrogen); 
		% npv.rawMaterialsCost = value_ethane(F_fresh_ethane);
		% npv.utilitiesCost = cost_steam(F_steam, COST_RATES_STEAM(STEAM_CHOICE, STEAM_COST_COL)); 
		% npv.CO2sustainabilityCharge = tax_C02(combusted_fuel_flow_rates, F_natural_gas); 
		% npv.conversion = conversion(i);
		% npv.isbl = cost_rxt_vec + cost_separation_system(P_flowrates, F_steam, R_ethane);

		% look at photos screenshot, cost of electricity, some shit 
	user = get_user_inputs();

	YEARS_IN_OPERATION = user.npv.period_plant_operation;

	% % Economic Assumptions 
	% npv.discountRate = 0.15;		% [ % in decimal ]
	% npv.taxRate = 0.27;				% [ % in decimal ]
	% npv.salvageValue = 0.05;		% [ % in decimal ]	

	npv.discountRate = user.npv.enterprise_rate;	% [ % in decimal ]
	npv.taxRate = user.npv.total_tax_rate;		% [ % in decimal ]
	npv.salvageValue = user.salvage_value;		% [ % in decimal ]	

	WORKING_CAP_PERCENT_OF_FCI = 0.15; 		% [ % in decimal ]
	STARTUP_COST_PERCENT_OF_FCI = 0.10;		% [ % in decimal ]
	LENGTH_CONSTRUCTION_TABLE = 6;
	LAST_ROW_CONSTRUCTION = LENGTH_CONSTRUCTION_TABLE; 
	YEARS_OF_CONSTUCTION = user.npv.construction_period;

	% Revenues & Production Costs	
	npv.consummablesCost = 0;
	npv.VCOP = npv.rawMaterialsCost + npv.utilitiesCost + ...
				npv.consummablesCost + npv.CO2sustainabilityCharge - ...
														npv.byProductRevenue;
	npv.salaryAndOverhead = 0;
	npv.maintenenace = 0;
	% npv.interest = 15;
	npv.interest = user.npv.enterprise_rate;

	% npv.AGS = (npv.mainProductRevenue + npv.byProductRevenue)*0.05;		% ~5% revenue
	npv.AGS = (npv.mainProductRevenue + npv.byProductRevenue) * ...
		user.npv.admin_and_general_services;	

	npv.FCOP = npv.salaryAndOverhead + npv.maintenenace +...
						 npv.AGS + npv.interest;

	% Capital Costs 
	npv.OSBLcapitalCost = npv.ISBLcapitalCost * 0.40;
		% ??  npv.ISBLcapitalCost, npv.OSBLcapitalCost)
	npv.contingency = (npv.ISBLcapitalCost + npv.OSBLcapitalCost) * 0.25;
	npv.indirectCost = (npv.ISBLcapitalCost + npv.OSBLcapitalCost + ...
													npv.contingency) * 0.30;
	npv.totalFixedCapitalCost = npv.ISBLcapitalCost + ...
								npv.OSBLcapitalCost + ...
								npv.indirectCost + ... 
								npv.contingency;

	npv.workingCapital = npv.totalFixedCapitalCost * WORKING_CAP_PERCENT_OF_FCI;
	npv.startupCost = npv.totalFixedCapitalCost * STARTUP_COST_PERCENT_OF_FCI;
	npv.land = 10;
	npv.totalCapitalInvestment = npv.totalFixedCapitalCost + ...
									npv.workingCapital + ...
									npv.startupCost + ...
									npv.land;


	% CONSTRUCTION SCHEDULE INDICIES 
	YEAR = 1;
	FC = 2;
	WC = 3;
	SU = 4;
	FCOP = 5;
	VCOP = 6;
	construction_matrix = zeros(LENGTH_CONSTRUCTION_TABLE + 1, VCOP);
	
	% Generate the construction schedule matrix
	for yr = 0:LENGTH_CONSTRUCTION_TABLE
		row = yr + 1;
		if yr > 0 && yr < 4
			construction_matrix(row, FC) = 0.33;
		end
		if yr == 3
			construction_matrix(row, WC) = 1.00;
			construction_matrix(row, SU) = 1.00;
		end
		if yr > 3 && yr <= 6
			construction_matrix(row, FCOP) = 1.00;
			construction_matrix(row, VCOP) = 1.00;
		end
	end

	% NPV COLUMN INDICIES 
	YEAR = 1;
	CAPITAL_EXPENSE = 2;
	REVENUE = 3;
	COM = 4;
	GROSS_PROFIT = 5;
	DEPRECIATION = 6;
	TAXABLE_INC = 7;
	TAXES_PAID = 8;
	CASH_FLOW = 9;
	CUM_CASH_FLOW = 10;
	PV_OF_CF = 11;
	CUM_PV_OF_CF = 12;
	NPV = 13;
	cash_flow_matrix = zeros(YEARS_IN_OPERATION + 1, NPV);
	LAST_ROW_CASHFLOW = YEARS_IN_OPERATION + 1; 


	for yr = 0:YEARS_IN_OPERATION
		row = yr + 1;
		cash_flow_matrix(row, YEAR) = yr;

		% Capital Expenses Column
		if yr == 0
			cash_flow_matrix(row, CAPITAL_EXPENSE) = npv.land;
		elseif yr >= 1 && yr <= 5
			cash_flow_matrix(row, CAPITAL_EXPENSE) ...
				= npv.totalFixedCapitalCost * construction_matrix(row,FC) + ...
				  npv.workingCapital * construction_matrix(row, WC) + ... 
				  npv.startupCost * construction_matrix(row, SU) ;
		elseif yr == YEARS_IN_OPERATION
			cash_flow_matrix(row, CAPITAL_EXPENSE) = - npv.salvageValue * npv.totalFixedCapitalCost;
		else
			cash_flow_matrix(row, CAPITAL_EXPENSE) ...
				= npv.totalFixedCapitalCost * construction_matrix(LAST_ROW_CONSTRUCTION,FC) + ...
				  npv.workingCapital * construction_matrix(LAST_ROW_CONSTRUCTION, WC) + ... 
				  npv.startupCost * construction_matrix(LAST_ROW_CONSTRUCTION, SU) ;
		end

		% Revenue Column
		if yr <= LENGTH_CONSTRUCTION_TABLE % ?? 
			cash_flow_matrix(row, REVENUE) = npv.mainProductRevenue * construction_matrix(row, VCOP);
		else
			cash_flow_matrix(row, REVENUE) = npv.mainProductRevenue * construction_matrix(LAST_ROW_CONSTRUCTION, VCOP);
		end

		% COM Column 
		if yr <= LENGTH_CONSTRUCTION_TABLE
			cash_flow_matrix(row, COM) = npv.VCOP * construction_matrix(row, VCOP) + ...
											npv.FCOP * construction_matrix(row, FCOP);
		else
			cash_flow_matrix(row, COM) = npv.VCOP * construction_matrix(LAST_ROW_CONSTRUCTION, VCOP) + ...
									npv.FCOP * construction_matrix(LAST_ROW_CONSTRUCTION, FCOP);
		end

		% Gross Profit
		cash_flow_matrix(row, GROSS_PROFIT) = cash_flow_matrix(row,REVENUE) - cash_flow_matrix(row, COM);
		
		% Depreciation
		if yr >= YEARS_OF_CONSTUCTION
			cash_flow_matrix(row, DEPRECIATION) = 0.1*(npv.totalFixedCapitalCost + npv.startupCost - 0.05*npv.totalFixedCapitalCost);
		end

		% Taxable Inc
		if yr >= YEARS_OF_CONSTUCTION
			cash_flow_matrix(row, TAXABLE_INC) = cash_flow_matrix(row, GROSS_PROFIT) - cash_flow_matrix(row,DEPRECIATION);
		end

		% Taxes Paid 
		if yr >= YEARS_OF_CONSTUCTION
			cash_flow_matrix(row, TAXES_PAID) = cash_flow_matrix(row, TAXABLE_INC) * npv.taxRate;
		end

		% Cash Flow
		cash_flow_matrix(row, CASH_FLOW) = -cash_flow_matrix(row, CAPITAL_EXPENSE) + ...
				( cash_flow_matrix(row,REVENUE) ...
					- cash_flow_matrix(row, COM) ...
					- cash_flow_matrix(row, DEPRECIATION) ... 
				) * ( 1 - npv.taxRate) + cash_flow_matrix(row, DEPRECIATION);
		
		% Cummulative Cash Flow
		cash_flow_matrix(row, CUM_CASH_FLOW) = sum( cash_flow_matrix( 1 : row, CASH_FLOW) );
		
		% PV of CF 
		cash_flow_matrix(row, PV_OF_CF) = cash_flow_matrix(row, CASH_FLOW) / ( 1 + npv.discountRate)^yr;

		% Cummulative PV of CF
		cash_flow_matrix(row , CUM_PV_OF_CF) = sum( cash_flow_matrix(1:row, PV_OF_CF) );

		% NPV
		if row > 1
			cash_flow_matrix(row , NPV) = cash_flow_matrix(row - 1, NPV) + cash_flow_matrix(row, PV_OF_CF);
		else
			cash_flow_matrix(row, NPV) = cash_flow_matrix(row, PV_OF_CF);
		end
	end

	% RETURN 
	cf.matrix = cash_flow_matrix;
	cf.lifetime_npv = cash_flow_matrix(LAST_ROW_CASHFLOW, NPV);
	
	lifetime_npv = cf.lifetime_npv;
end
