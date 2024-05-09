clc; clear; close all; 



% SCRIPT________________________________________________________________________

level_2();

% FUNCTIONS_____________________________________________________________________


function void = level_2()

    % F = get_feed_stream();

    const = get_constants(); 

    user = get_user_inputs();

    F_fxns = flowrate_fxns();

    F = user.level2.feed_stream;
    

    for s = user.level2.selectivity_range();
        


    end
    
    disp("Level 2 calculations are complete")

end