
function fxns = separation_fxns()
	fxns.work_min = @get_work_min

end

function w_min = get_work_min(F)
	% Assumes that all separation products streams are pure
	% w_min has units of watts (J / s)
	const = get_constants();
	R = const.thermo.R;
	T = const.temp.c_to_k(F.T);

	w_min = 0;	
	fieldNames = fieldnames(F);
	for i = 1:length(fieldNames)
		species = fieldNames{i};
		x = 1; % assumes that each separated flowstream is pure
		z = F.(species).x;
		F_i = F.(species).mol;
		w_min = w_min + (F_i * R * T * x * log(x/z));
		
	end

end



function [sep_top1, sep_btm1] = flash_v100(sep)

	sep.heat = 0;

	K.ethane = 3.760 * 10^9;
	K.ethylene = 7.266 * 10^8;
	K.hydrogen = 3.193 * 10^6;
	K.methane = 8.488 * 10^7;
	K.propane = 5.252 * 10^11;
	K.butane = 3.978 * 10^14;
	K.water = 1.561 * 10^-2;

	[sep_top1, sep_btm1]= rachford_rice(sep, K);
	
end

function phi = underwood(z, r, s, alpha, y, x, q)
	% y are the distillate compositions 
	% alpha has the relative volatilities 
	% r is reflux ratio
	% s is boilup ratio 
	r_min_factor = 1.2;

	% Doherty & Malone eq 4.21  
	eqn_421_top = @(phi, r) -r - 1 + (alpha.a * y.a / (alpha.a - phi)) + (alpha.b * y.b / (alpha.b - phi));
	eqn_421_bot = @(phi, s) s + (alpha.a * x.a / (alpha.a - phi)) + (alpha.b * x.b / (alpha.b - phi));
	
	init_phi = alpha.b + (alpha.a / 2);
		% alpha_a > phi > alpha_b 
	phi_1_top =  fzero( @(phi) eqn_421_top(phi, r), init_phi);

	init_phi = alpha.b / 2;
		% alpha_b > phi > 0
	phi_2_top = fzero( @(phi) eqn_421_top(phi,r), init_phi);

	% Doherty and Malone eq 4.25
	term1 = alpha.a * z.a / (alpha.a - phi_2_top);
	term2 = alpha.b * z.b / (alpha.b - phi_2_top);
	term3 = alpha.a * z.a / (alpha.a - phi_1_top);
	term4 = alpha.b * z.b / (alpha.b - phi_1_top);
	trays_above_feed = log((term1 + term2) / (term3 + term4)) / log(phi_1_top / phi_2_top);


	init_phi = alpha.a * 2 ;
	init_phi = alpha.b * (0.5 * (alpha.a - alpha.b));
		% alpha_a > phi > alpha_b
	phi_2_bar =  fzero( @(phi) eqn_421_bot(phi, s), init_phi);

	init_phi = alpha.a * 1.5;
		% inf > phi > alpha_a
	phi_1_bar = fzero( @(phi) eqn_421_bot(phi, s), init_phi);

	term1 = alpha.a * z.a / (alpha.a - phi_1_bar);
	term2 = alpha.b * z.b / (alpha.b - phi_1_bar);
	term3 = alpha.a * z.a / (alpha.a - phi_2_bar);
	term4 = alpha.b * z.b / (alpha.b - phi_2_bar);
	trays_below_feed = log((term1 + term2) / (term3 + term4)) / log(phi_1_bar / phi_2_bar);

	% Doherty and Malone eq4.29
	find_theta = @(theta) q - 1 + (alpha.a * z.a / (alpha.a - theta)) + ...
								(alpha.b * z.b / (alpha.b - theta));
	init_theta = phi_1_top + (phi_2_bar - phi_1_top)*0.5; % ??? 
	theta = fzero( @(theta) find_theta(theta), init_theta);
	
	r_min = -1 + (alpha.a * y.a / (alpha.a - theta)) + (alpha.b * y.b / (alpha.b - theta));
	r = r_min_factor * r_min; % WHAT MULTIPLE OF R_MIN SHOULD WE USE? 
	phi = 0;
end


function W = compressor_work_TJ(sep, P_f)
	R = 8.314;		% [ J / mol K]
	n = total_molar_flowrate(sep.F);
	T = sep.T;
	P_0 = sep.P;

	W = - n * R * T * log10(P_f / P_0);	% ?? I think that it's base 10

	
	% ?? THIS ALWAYS RETURNS 0 OR NULL, NOT IMPLEMENTED YET
end



function T_f = adiabatic_temp(T_0, P_0, P_f)

	T_f = T_0 * ( P_0 / P_f);
end



function [sep_top, sep_btm]= rachford_rice(sep, K)
	global KT_PER_G MOLMASS_BUTANE MOLMASS_ETHANE MOLMASS_ETHYLENE MOLMASS_HYDROGEN MOLMASS_METHANE MOLMASS_PROPANE MOLMASS_WATER

	sep.x = all_mol_fractions(sep.F);	
	
	f_phi = @(phi, sep, K) ((sep.x.methane * (K.methane - 1)) / (1 + phi*(K.methane - 1))) + ...
		((sep.x.ethane * (K.ethane - 1)) / (1 + phi*(K.ethane - 1))) + ...
		((sep.x.ethylene * (K.ethylene - 1)) / (1 + phi*(K.ethylene - 1))) + ...
		((sep.x.hydrogen * (K.hydrogen - 1)) / (1 + phi*(K.hydrogen - 1))) + ...
		((sep.x.propane * (K.propane - 1)) / (1 + phi*(K.propane - 1))) + ...
		((sep.x.butane * (K.butane - 1)) / (1 + phi*(K.butane - 1))) + ...
		((sep.x.water * (K.water - 1)) / (1 + phi*(K.water - 1)));

	init_cond = 0.5;

	phi = fzero(@(phi) f_phi(phi, sep, K), init_cond);

	if phi > 1
		phi = 1;
	elseif phi < 0
		phi = 0;
	end
	

	% Liquid compositions 
	x.hydrogen = sep.x.hydrogen / (1 + phi*(K.hydrogen - 1));
	x.methane = sep.x.methane / (1 + phi*(K.methane - 1));
	x.ethane = sep.x.ethane / (1 + phi*(K.ethane - 1));
	x.ethylene = sep.x.ethylene / (1 + phi*(K.ethylene - 1));
	x.propane = sep.x.propane / (1 + phi*(K.propane - 1));
	x.butane = sep.x.butane / (1 + phi*(K.butane - 1));
	x.water = sep.x.water / (1 + phi*(K.water - 1));

	% Vapor compositions 
	y.hydrogen = K.hydrogen * x.hydrogen;
	y.methane = K.methane * x.methane;
	y.ethane = K.ethane * x.ethane;
	y.ethylene = K.ethylene * x.ethylene;
	y.propane = K.propane * x.propane;
	y.butane = K.butane * x.butane;
	y.water = K.water * x.water;

	% Splitting
	F_tot = total_molar_flowrate(sep.F);
	V = phi * F_tot; 
	L = (1 - phi) * F_tot; 

	% Tops 
	sep_top = sep;
	sep_top.y = y;
	sep_top.x = NaN;

	
	% kta             = (mol/yr) * (mol / mol) * (g / mol)   * (kt / g)
	sep_top.F.hydrogen = V * y.hydrogen * (MOLMASS_HYDROGEN) * KT_PER_G;
	sep_top.F.methane = V * y.methane * (MOLMASS_METHANE) * KT_PER_G;
	sep_top.F.ethane = V * y.ethane * (MOLMASS_ETHANE) * KT_PER_G;
	sep_top.F.ethylene = V * y.ethylene * (MOLMASS_ETHYLENE) * KT_PER_G;
	sep_top.F.propane = V * y.propane * (MOLMASS_PROPANE) * KT_PER_G;
	sep_top.F.butane = V * y.butane * (MOLMASS_BUTANE) * KT_PER_G;
	sep_top.F.water = V * y.water * (MOLMASS_WATER) * KT_PER_G;
	
	% Bottoms 
	sep_btm = sep; 
	sep_btm.x = x;
	sep_btm.y = NaN;

	% kta     = (mol/yr) * (mol / mol) * (g / mol)   * (kt / g)
	sep_btm.F.hydrogen = L * x.hydrogen * (MOLMASS_HYDROGEN) * KT_PER_G;
	sep_btm.F.methane = L * x.methane * (MOLMASS_METHANE) * KT_PER_G;
	sep_btm.F.ethane = L * x.ethane * (MOLMASS_ETHANE) * KT_PER_G;
	sep_btm.F.ethylene = L * x.ethylene * (MOLMASS_ETHYLENE) * KT_PER_G;
	sep_btm.F.propane = L * x.propane * (MOLMASS_PROPANE) * KT_PER_G;
	sep_btm.F.butane = L * x.butane * (MOLMASS_BUTANE) * KT_PER_G;
	sep_btm.F.water = L * x.water * (MOLMASS_WATER) * KT_PER_G;
	
	
end