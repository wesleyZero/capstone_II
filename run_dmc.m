clc; clear; close all; 

% SCRIPT________________________________________________________________________

level2();
level3();

% FUNCTIONS_____________________________________________________________________

function void = level3()

    console = get_console();
    const = get_constants(); 
    user = get_user_inputs();
    F_fxns = flowrate_fxns();
    F = user.level2.feed_stream;

    console.section("Starting Level 3 calculations")

    for s = user.level2.selectivity_range
        for chi = user.level3.conversion_range

            F = level3_flowrates(F, s, chi);
        end

    end
    
    console.section("Level 3 calculations are complete")

end

function [F, F_rxtr, R] = level3_flowrates(F, s, chi)
    % I need to really check the F_fresh and P flows because I think I am messing them up
    user = get_user_inputs();

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
    F.ethylene_carbonate.mol = 0 ;
        % Ethylene carbonate is recycled to completion
        % ?? Level 3 tho ??? 
    F = flowrate_fxns().set_F_mol(F);

    R = get_recycle_flowrates(F, s, chi);

    F_rxtr = get_reactor_flowrates(F, R);
end

function F_rxtr = get_reactor_flowrates(F, R)
    F_rxtr = flowrate_fxns().get_blank_flowstream();
    F_rxtr.ethylene_carbonate.mol = F.ethylene_carbonate.mol + R.ethylene_carbonate.mol;
    F_rxtr.ethylene_oxide.mol = F.ethylene_carbonate.mol + R.ethylene_carbonate.mol;
    F_rxtr.methanol.mol = F.methanol.mol + R.methanol.mol;
    F_rxtr.carbon_dioxide.mol = F.carbon_dioxide.mol + R.carbon_dioxide.mol;
    F_rxtr = flowrate_fxns().set_F_mol(F_rxtr);
end

function F = get_recycle_flowrates(F, s, chi)
    user = get_user_inputs();

    R = flowrate_fxns().get_blank_flowstream(); 
    MR = user.level3.molar_ratio_methanol;
    R.methanol.mol = F.dimethyl_carbonate.mol * ((MR/s) - 1 - (1/s));
    MR = user.level3.molar_ratio_carbon_dioxide;
    R.carbon_dioxide.mol = F.dimethyl_carbonate.mol * (((MR + 1)/s) - 1);
    R.ethylene_carbonate.mol = (F.dimethyl_carbonate.mol / s) * ( (1-chi) / chi);
    R = flowrate_fxns().set_F_mol(R);
end

function void = level2()
    console = get_console();
    const = get_constants(); 
    user = get_user_inputs();
    F_fxns = flowrate_fxns();
    F = user.level2.feed_stream;

    console.section("Starting Level 2 calculations")

    for s = user.level2.selectivity_range
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
    F.ethylene_carbonate.mol = 0 ;
        % Ethylene carbonate is recycled to completion

    % Update the flowstream
    F = flowrate_fxns().set_F_mol(F);
end
