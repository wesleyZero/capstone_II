clc; clear; close all; 



% SCRIPT________________________________________________________________________

level2();

% FUNCTIONS_____________________________________________________________________


function void = level2()

    console = get_console();
    const = get_constants(); 
    user = get_user_inputs();
    F_fxns = flowrate_fxns();
    F = user.level2.feed_stream;

    console.section("Starting Level 2 calculations")

    for s = user.level2.selectivity_range()

        F = level2_flowrates(F, s);

    end
    
    console.section("Level 2 calculations are complete")

end


function F = level2_flowrates(F, s)

    % Specify Constraints in kta
    F.dimethyl_carbonate.kta = 100;

    % Update the flowstream  
    F = flowrate_fxns().set_F_kta(F);
    
    % Calculate the molar flowrates as a result of reaction
    F.ethylene_glycol.mol = F.dimethyl_carbonate.mol;
    F.carbon_dioxide.mol = F.dimethyl_carbonate.mol;
    F.ethylene_oxide.mol = F.dimethyl_carbonate.mol / s;
    F.methoxy_ethanol.mol = F.dimethyl_carbonate.mol * ((1/s) - 1);
    F.methanol.mol = F.dimethyl_carbonate.mol * ( 1 + (1/s));

    % Update the flowstream
    F = flowrate_fxns().set_F_mol(F);

end
