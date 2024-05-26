
function user = get_user_inputs()
 
    % Level 3 | hard coded
    user.level3.molar_ratio_methanol_EO = 15;
    user.level3.molar_ratio_carbon_dioxide_EO = 12; 
        % one mol equiv is consumed in the virtual reactor
    user.level3.tau_precision = 1000;
    user.level3.tau_max = 10^4;
    user.level3.temp_precision = 10;
    user.level3.temp_min.C = 80;
    user.level3.temp_max.C = 140;
    user.level3.pressure_precision = 10;
    user.level3.press_min.bar = 80;
    user.level3.press_max.bar = 150;
    user.level3.isothermal_temp.C = 140;
    user.level3.isobaric_press.bar = 150;
    % user.level3.fsolveOpt = optimoptions('fsolve', 'Display', 'off'); 
    user.level3.fsolveOpt = optimoptions('fsolve', ...
        'Algorithm', 'trust-region-dogleg', ...       % Choose the algorithm
        'Display', 'off', ...                         % Control display level
        'FunctionTolerance', 1e-8, ...                % Termination tolerance on the function value
        'StepTolerance', 1e-8, ...                    % Termination tolerance on the step size
        'MaxFunctionEvaluations', 1000, ...           % Maximum number of function evaluations allowed
        'MaxIterations', 500, ...                     % Maximum number of iterations allowed
        'OptimalityTolerance', 1e-8, ...              % Termination tolerance on the first-order optimality
        'UseParallel', false, ...                     % Use parallel computing
        'CheckGradients', true, ...                  % Check gradients numerically
        'FiniteDifferenceStepSize', 1e-3, ...         % Step size for finite differences
        'FiniteDifferenceType', 'central', ...        % Finite difference type ('forward' or 'central')
        'FunctionTolerance', 1e-6, ...                % Termination tolerance on the function value
        'PlotFcn', [], ...                            % Plot functions
        'OutputFcn', []);                             % Output functions
    user.dmc_production_rate = 100; % [kta]

    % Level 3 | Soft coded
    user.level3.tau_range = ...
                linspace(1, user.level3.tau_max, user.level3.tau_precision);
    % user.level3.conversion_range = ...
    %             linspace(1/user.level3.precision, 1, user.level3.precision);
    user.level3.temp_range = ...
        linspace(user.level3.temp_min.C, user.level3.temp_max.C, ...
                    user.level3.temp_precision);
    user.level3.press_range = ...
        linspace(user.level3.press_min.bar, user.level3.press_max.bar, ...
                    user.level3.pressure_precision);
end 
