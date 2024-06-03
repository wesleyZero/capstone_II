
function fxns = get_economic_functions()
    fxns.get_npv = @get_npv;
	fxns.get_work_min_npv = @get_work_min_npv;

end

function lifetime_npv = get_work_min_npv(F, T)
	sep_fxns = separation_fxns();
	w_min = sep_fxns.get_work_min(F, T)

	lifetime_npv = -2;

end

function value = get_main_product_revenue()
	% Products are the value of DMC
	value = 0;
end

function value = get_biproduct_revenue()
	% byproducts are EG
	value = 0;

end

function cost = get_raw_material_cost()
	% raw material costs are EO, MeOH, C02, water
	% ask TJ how to cost the water for the sep unit

	cost = 0;
end

function cost = get_utilities_cost() 
	% utilities are electricity, fuel??

	cost = 0;
end

function value = get_CO2_sustainability_value()
	% 45 / MT, confirm with TJ

	value = 0;

end


function lifetime_npv = get_npv(npv)
	% USER_INPUTS | All inputs are in units of $MM
		% npv.mainProductRevenue = value_ethylene(P_ethylene);
		% npv.byProductRevenue = value_h2_chem(P_hydrogen - combusted_hydrogen); 
		% npv.rawMaterialsCost = value_ethane(F_fresh_ethane);
		% npv.utilitiesCost = cost_steam(F_steam, COST_RATES_STEAM(STEAM_CHOICE, STEAM_COST_COL)); 
		% npv.CO2sustainabilityCharge = tax_C02(combusted_fuel_flow_rates, F_natural_gas); 
		% npv.conversion = conversion(i);
		% npv.isbl = cost_rxt_vec + cost_separation_system(P_flowrates, F_steam, R_ethane);

	const = get_user_inputs();

	YEARS_IN_OPERATION = const.user.npv.period_plant_operation;

	% % Economic Assumptions 
	% npv.discountRate = 0.15;		% [ % in decimal ]
	% npv.taxRate = 0.27;				% [ % in decimal ]
	% npv.salvageValue = 0.05;		% [ % in decimal ]	

	npv.discountRate = const.user.npv.enterprise_rate;	% [ % in decimal ]
	npv.taxRate = const.user.npv.total_tax_rate;		% [ % in decimal ]
	npv.salvageValue = const.user.salvage_value;		% [ % in decimal ]	

	WORKING_CAP_PERCENT_OF_FCI = 0.15; 		% [ % in decimal ]
	STARTUP_COST_PERCENT_OF_FCI = 0.10;		% [ % in decimal ]
	LENGTH_CONSTRUCTION_TABLE = 6;
	LAST_ROW_CONSTRUCTION = LENGTH_CONSTRUCTION_TABLE; 
	YEARS_OF_CONSTUCTION = const.user.npv.construction_period;

	% Revenues & Production Costs	
	npv.consummablesCost = 0;
	npv.VCOP = npv.rawMaterialsCost + npv.utilitiesCost + ...
				npv.consummablesCost + npv.CO2sustainabilityCharge - ...
														npv.byProductRevenue;
	npv.salaryAndOverhead = 0;
	npv.maintenenace = 0;
	npv.interest = 15;
	npv.AGS = (npv.mainProductRevenue + npv.byProductRevenue)*0.05;		% ~5% revenue
	npv.FCOP = npv.salaryAndOverhead + npv.maintenenace +...
						 npv.AGS + npv.interest;

	% Capital Costs 
	npv.OSBLcapitalCost = npv.ISBLcapitalCost * 0.40;
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
	
	lifetime_npv = -1;
end
