
function user = get_user_inputs()
 

    user.level2.selectivity_range = linspace(0.1, 1,100);
    F = get_feed_flowrates();
    user.level2.feed_stream = flowrate_fxns().set_F_kta(F);    
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
