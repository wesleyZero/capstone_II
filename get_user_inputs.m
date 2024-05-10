
function user = get_user_inputs()
 
    % Level 2
    user.level2.precision = 100;
    user.level2.selectivity_range = ...
        linspace(1/user.level2.precision, 1, user.level2.precision);
    F = get_feed_flowrates();
    user.level2.feed_stream = flowrate_fxns().set_F_kta(F);    

    % Level 3
    user.level3.precision = 100;
    user.level3.conversion_range = ...
        linspace(1/user.level3.precision, 1, user.level3.precision);

    

end 

function F = get_feed_flowrates()
    
    F.carbon_dioxide.kta = 1;
    F.ethylene_oxide.kta = 1; 
    F.methanol.kta = 1;
    F.ethylene_carbonate.kta = 0;
    F.ethylene_glycol.kta = 0;
    F.methoxy_ethanol.kta = 0;
    F.dimethyl_carbonate.kta = 0;

    F = flowrate_fxns().set_F_kta(F);

end 
