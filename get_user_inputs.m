
function user = get_user_inputs()
 
    % Level 3 | hard coded
    user.level3.molar_ratio_methanol_EO = 15;
    user.level3.molar_ratio_methanol_EC = 15;
    user.level3.molar_ratio_carbon_dioxide_EO = 13; 
    user.level3.molar_ratio_carbon_dioxide_EC = 12; 
        % one mol equiv is consumed in the virtual reactor
    user.level3.tau_precision = 500;
    % user.level3.tau_min = 50;
    % user.level3.tau_max = 500; 

    % isothermal tau ranges
    user.level3.isothermal.tau_min = 1;
    user.level3.isothermal.tau_max = 500;
    user.level3.isothermal.P_specify.tau_min = 50;
    user.level3.isothermal.P_specify.tau_max = 500;
    user.level3.isobaric.tau_min = 1;
    user.level3.isobaric.tau_max = 2000;


    user.level3.temp_precision = 10;
    user.level3.temp_min.C = 80;
    user.level3.temp_max.C = 140;

    user.level3.pressure_precision = 10;
    user.level3.press_min.bar = 50;
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

    % Plotting 
    user.plot.image_dpi = '-r600'; % [600 dpi]
    user.plot.isothermal.x_point = 0.65;
    user.plot.isothermal.y_point = 2 * 10^4;
    user.plot.supercritical_c02_pressure = 73.8;    %  Bar

    % Level 3 | Soft coded
    % user.level3.tau_range = ...
                % linspace(user.level3.isobaric.tau_min, user.level3.isobaric.tau_max, user.level3.tau_precision);

    user.level3.tau_range.isobaric = ...
                linspace(user.level3.isobaric.tau_min, user.level3.isobaric.tau_max, user.level3.tau_precision);

    user.level3.tau_range.isothermal = ...
                linspace(user.level3.isothermal.tau_min, user.level3.isothermal.tau_max, user.level3.tau_precision);

    user.level3.tau_range.P_specify.isothermal = ...
        linspace(user.level3.isothermal.P_specify.tau_min, user.level3.isothermal.P_specify.tau_max, user.level3.tau_precision);







    user.level3.temp_range = ...
        linspace(user.level3.temp_min.C, user.level3.temp_max.C, ...
                    user.level3.temp_precision);
    user.level3.press_range = ...
        linspace(user.level3.press_min.bar, user.level3.press_max.bar, ...
                    user.level3.pressure_precision);

    % NPV 
    user.npv.enterprise_rate = 0.15;                % per year
    user.npv.total_tax_rate = 0.27;                 % per year
    user.npv.admin_and_general_services = 0.05;     % 5% of annual revenues 
    user.npv.construction_period = 3;               % years
    user.npv.period_plant_operation = 12;           % years
    user.npv.project_life = ...
        user.npv.construction_period + user.npv.period_plant_operation;
                                                    % years (constr + op yrs)
    user.depreciation_schedule = 10;                % years, linear deprec.
    user.marshall_and_swift_index = 1800;           % M&S 
    user.salvage_value = 0.05;                      % 5% of fixed capital investment

    % user.npv.discount_rate = 
end 
