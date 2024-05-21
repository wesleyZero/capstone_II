
function user = get_user_inputs()
 
    % Level 3 | hard coded
    user.level3.precision = 100;
    user.level3.molar_ratio_methanol_EO = 15;
    user.level3.molar_ratio_carbon_dioxide_EO = 12; 
        % one mol equiv is consume in the virtual reactor
    user.level3.tau_precision = 100;
    user.level3.tau_max = 100;
    user.level3.temp_precision = 10;
    user.level3.temp_min.C = 80;
    user.level3.temp_max.C = 140;
    user.level3.pressure_precision = 10;
    user.level3.press_min.bar = 80;
    user.level3.press_max.bar = 150;
    user.level3.isothermal_temp.C = 140;
    user.level3.isobaric_press.bar = 150;

    % Level 3 | Soft coded
    user.level3.tau_range = ...
                linspace(1, user.level3.tau_max, user.level3.tau_precision);
    user.level3.conversion_range = ...
                linspace(1/user.level3.precision, 1, user.level3.precision);
    user.level3.temp_range = ...
        linspace(user.level3.temp_min.C, user.level3.temp_max.C, ...
                    user.level3.temp_precision);
    user.level3.press_range = ...
        linspace(user.level3.press_min.bar, user.level3.press_max.bar, ...
                    user.level3.pressure_precision);
end 
