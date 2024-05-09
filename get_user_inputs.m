



function user = get_user_inputs()
 

    user.level2.selectivity_range = linspace(0, 100,100);
    
    F.carbon_dioxide.kta = 1;
    F.ethylene_oxide.kta = 1; 
    F.methanol.kta = 1;
    F.ethylene_carbonate.kta = 0;
    F.ethylene_glycol.kta = 0;
    F.methoxy_ethanol.kta = 0;

    F_fxns = flowrate_fxns();
    F = F_fxns.set_F_kta(F);    
    user.level2.feed_stream = F; 

end 