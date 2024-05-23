clc; clear; close all; 

% SCRIPT________________________________________________________________________

level3_isobaric();

% FUNCTIONS_____________________________________________________________________

function void = level3_isobaric()

    console = get_console();
    const = get_constants(); 
    user = get_user_inputs();
    F_fxns = flowrate_fxns();
    P = user.level3.isobaric_press.bar;
    opt = 'isobaric';
    console.section("Starting Level 3 calculations")

    for T = user.level3.temp_range
        console.subsection(sprintf("T = %3.2f", T), 1);
        for tau = user.level3.tau_range
            [F_in, F_out, R] = level3_flowrates(tau, T, P, opt); 
        end
    end
    
    console.section("Level 3 calculations are complete")
end


function [F_in, F_out, R] = level3_flowrates(tau, temp, P, opt)
    user = get_user_inputs(); 
    flow_fxns = flowrate_fxns();
    rxtr_fxns = reactor_fxns();

    % Basis calculations 
    F_basis = flow_fxns.get_basis_feed_flowrates();
    [F_in, F_out, R] = rxtr_fxns.get_reactor_flows(F_basis, temp, P, opt, tau);
    

end



% function [F, F_rxtr, R] = level3_flowrates(F, s, chi)
%     % I need to really check the F_fresh and P flows because I think I am messing them up
%     user = get_user_inputs();

%     % Specify Constraints in kta
%     F.dimethyl_carbonate.kta = 100;

%     % Update the flowstream  
%     F = flowrate_fxns().set_F_kta(F);
    
%     % Calculate the molar flowrates as a result of reaction
%     F.ethylene_glycol.mol = F.dimethyl_carbonate.mol;
%     F.carbon_dioxide.mol = F.dimethyl_carbonate.mol;
%     F.ethylene_oxide.mol = F.dimethyl_carbonate.mol / s;
%     F.methoxy_ethanol.mol = F.dimethyl_carbonate.mol * ((1/s) - 1);
%     F.methanol.mol = F.dimethyl_carbonate.mol * ( 1 + (1/s));
%     F.ethylene_carbonate.mol = 0 ;
%         % Ethylene carbonate is recycled to completion
%         % ?? Level 3 tho ??? 
%     F = flowrate_fxns().set_F_mol(F);
%     % I SHOULD CALL THIS P

%     R = get_recycle_flowrates(F, s, chi);

%     F_rxtr = get_reactor_flowrates(F, R);
% end

% function F_rxtr = get_reactor_flowrates(F, R)
%     % Depreciated I think 
%     F_rxtr = flowrate_fxns().get_blank_flowstream();
%     F_rxtr.ethylene_carbonate.mol = F.ethylene_carbonate.mol + R.ethylene_carbonate.mol;
%     F_rxtr.ethylene_oxide.mol = F.ethylene_carbonate.mol + R.ethylene_carbonate.mol;
%     F_rxtr.methanol.mol = F.methanol.mol + R.methanol.mol;
%     F_rxtr.carbon_dioxide.mol = F.carbon_dioxide.mol + R.carbon_dioxide.mol;
%     F_rxtr = flowrate_fxns().set_F_mol(F_rxtr);
% end

% function F = get_recycle_flowrates(F, s, chi)
%     % Depreciated I think 
%     user = get_user_inputs();

%     R = flowrate_fxns().get_blank_flowstream(); 
%     MR = user.level3.molar_ratio_methanol;
%     R.methanol.mol = F.dimethyl_carbonate.mol * ((MR/s) - 1 - (1/s));
%     MR = user.level3.molar_ratio_carbon_dioxide;
%     R.carbon_dioxide.mol = F.dimethyl_carbonate.mol * (((MR + 1)/s) - 1);
%     R.ethylene_carbonate.mol = (F.dimethyl_carbonate.mol / s) * ( (1-chi) / chi);
%     R = flowrate_fxns().set_F_mol(R);
% end
