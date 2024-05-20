
function user = get_user_inputs()
 
    % Level 2
    user.level2.precision = 100;
    user.level2.selectivity_range = ...
        linspace(1/user.level2.precision, 1, user.level2.precision);

    % Level 3 | hard coded
    user.level3.precision = 100;
    user.level3.molar_ratio_methanol_EO = 15;
    user.level3.molar_ratio_carbon_dioxide_EO = 13; 

    % Level 3 | Soft coded
    user.level3.conversion_range = ...
        linspace(1/user.level3.precision, 1, user.level3.precision);
end 
