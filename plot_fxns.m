


function fxns = plot_fxns()

    fxns.get_empty_vector = @get_empty_vector;
    fxns.get_plot_struct = @get_plot_structure;
    fxns.set_plot_row = @set_plot_row;
    fxns.plot_reactor_volume_conversion = @plot_reactor_volume_conversion;
    fxns.plot_reactor_volume_conversion_allT = @plot_reactor_volume_conversion_allT;
end


function void = plot_reactor_volume_conversion_allT(all_temp_data)
    user = get_user_inputs();
    figure 
    hold on
    for i = 1:length(all_temp_data)
        plot_struct = all_temp_data(i);
        x = plot_struct.data.conversion(:);
        y = plot_struct.data.V_rxtr(:);

        % figure
        plot(x, y);
        % plt_title = sprintf('V_{rxtr} [L] %3.0f [Bar]', plot_struct.P);
        title(sprintf('V_{rxtr} [L] %3.0f [Bar]', plot_struct.P), 'Interpreter', 'tex');
        xlabel('\chi', 'Interpreter', 'tex');
        ylabel('V_{rxtr} [L]', 'Interpreter', 'tex')
        % Create the legend entry for this plot
        legendEntries{i} = sprintf('%3.0f°C', plot_struct.T);
    end
    legend(legendEntries, 'Interpreter', 'tex', 'location', 'northwest');
    hold off

    dir = [pwd  '/plots/' ];
    print(fullfile(dir, "isobaric_V_reactor"), '-dpng', user.plot.image_dpi) ;  % Save as PNG with 300 DPI

end 



function void = plot_reactor_volume_conversion(plot_struct)

    x = plot_struct.data.conversion(:);
    y = plot_struct.data.V_rxtr(:);

    hold on
    figure
    plot(x, y);
    title(sprintf('V_{rxtr} [L] vs \\chi %3.0f [C] %3.0f [Bar]', plot_struct.T, plot_struct.P), 'Interpreter', 'tex');

    xlabel('\chi', 'Interpreter', 'tex');
    ylabel('V_{rxtr} [L]', 'Interpreter', 'tex')
    hold off
end 

function plot_struct = set_plot_row(plot_struct, plot_row)

    if isnan(plot_row.conversion)
        return
    end

    % F_fresh.mol
    plot_struct.data.F_fresh.carbon_dioxide.mol(plot_row.row_number, 1) = plot_row.F_fresh.carbon_dioxide.mol;
    plot_struct.data.F_fresh.ethylene_oxide.mol(plot_row.row_number, 1) = plot_row.F_fresh.ethylene_oxide.mol;
    plot_struct.data.F_fresh.methanol.mol(plot_row.row_number, 1) = plot_row.F_fresh.methanol.mol;
    plot_struct.data.F_fresh.ethylene_carbonate.mol(plot_row.row_number, 1) = plot_row.F_fresh.ethylene_carbonate.mol;
    plot_struct.data.F_fresh.ethylene_glycol.mol(plot_row.row_number, 1) = plot_row.F_fresh.ethylene_glycol.mol;
    plot_struct.data.F_fresh.methoxy_ethanol.mol(plot_row.row_number, 1) = plot_row.F_fresh.methoxy_ethanol.mol;
    plot_struct.data.F_fresh.dimethyl_carbonate.mol(plot_row.row_number, 1) = plot_row.F_fresh.dimethyl_carbonate.mol;

    % F_fresh.kta
    plot_struct.data.F_fresh.carbon_dioxide.kta(plot_row.row_number, 1) = plot_row.F_fresh.carbon_dioxide.kta;
    plot_struct.data.F_fresh.ethylene_oxide.kta(plot_row.row_number, 1) = plot_row.F_fresh.ethylene_oxide.kta;
    plot_struct.data.F_fresh.methanol.kta(plot_row.row_number, 1) = plot_row.F_fresh.methanol.kta;
    plot_struct.data.F_fresh.ethylene_carbonate.kta(plot_row.row_number, 1) = plot_row.F_fresh.ethylene_carbonate.kta;
    plot_struct.data.F_fresh.ethylene_glycol.kta(plot_row.row_number, 1) = plot_row.F_fresh.ethylene_glycol.kta;
    plot_struct.data.F_fresh.methoxy_ethanol.kta(plot_row.row_number, 1) = plot_row.F_fresh.methoxy_ethanol.kta;
    plot_struct.data.F_fresh.dimethyl_carbonate.kta(plot_row.row_number, 1) = plot_row.F_fresh.dimethyl_carbonate.kta;

    % F_rxtr.mol
    plot_struct.data.F_rxtr.carbon_dioxide.mol(plot_row.row_number, 1) = plot_row.F_rxtr.carbon_dioxide.mol;
    plot_struct.data.F_rxtr.ethylene_oxide.mol(plot_row.row_number, 1) = plot_row.F_rxtr.ethylene_oxide.mol;
    plot_struct.data.F_rxtr.methanol.mol(plot_row.row_number, 1) = plot_row.F_rxtr.methanol.mol;
    plot_struct.data.F_rxtr.ethylene_carbonate.mol(plot_row.row_number, 1) = plot_row.F_rxtr.ethylene_carbonate.mol;
    plot_struct.data.F_rxtr.ethylene_glycol.mol(plot_row.row_number, 1) = plot_row.F_rxtr.ethylene_glycol.mol;
    plot_struct.data.F_rxtr.methoxy_ethanol.mol(plot_row.row_number, 1) = plot_row.F_rxtr.methoxy_ethanol.mol;
    plot_struct.data.F_rxtr.dimethyl_carbonate.mol(plot_row.row_number, 1) = plot_row.F_rxtr.dimethyl_carbonate.mol;

    % F_rxtr.kta
    plot_struct.data.F_rxtr.carbon_dioxide.kta(plot_row.row_number, 1) = plot_row.F_rxtr.carbon_dioxide.kta;
    plot_struct.data.F_rxtr.ethylene_oxide.kta(plot_row.row_number, 1) = plot_row.F_rxtr.ethylene_oxide.kta;
    plot_struct.data.F_rxtr.methanol.kta(plot_row.row_number, 1) = plot_row.F_rxtr.methanol.kta;
    plot_struct.data.F_rxtr.ethylene_carbonate.kta(plot_row.row_number, 1) = plot_row.F_rxtr.ethylene_carbonate.kta;
    plot_struct.data.F_rxtr.ethylene_glycol.kta(plot_row.row_number, 1) = plot_row.F_rxtr.ethylene_glycol.kta;
    plot_struct.data.F_rxtr.methoxy_ethanol.kta(plot_row.row_number, 1) = plot_row.F_rxtr.methoxy_ethanol.kta;
    plot_struct.data.F_rxtr.dimethyl_carbonate.kta(plot_row.row_number, 1) = plot_row.F_rxtr.dimethyl_carbonate.kta;

    % F_out.mol
    plot_struct.data.F_out.carbon_dioxide.mol(plot_row.row_number, 1) = plot_row.F_out.carbon_dioxide.mol;
    plot_struct.data.F_out.ethylene_oxide.mol(plot_row.row_number, 1) = plot_row.F_out.ethylene_oxide.mol;
    plot_struct.data.F_out.methanol.mol(plot_row.row_number, 1) = plot_row.F_out.methanol.mol;
    plot_struct.data.F_out.ethylene_carbonate.mol(plot_row.row_number, 1) = plot_row.F_out.ethylene_carbonate.mol;
    plot_struct.data.F_out.ethylene_glycol.mol(plot_row.row_number, 1) = plot_row.F_out.ethylene_glycol.mol;
    plot_struct.data.F_out.methoxy_ethanol.mol(plot_row.row_number, 1) = plot_row.F_out.methoxy_ethanol.mol;
    plot_struct.data.F_out.dimethyl_carbonate.mol(plot_row.row_number, 1) = plot_row.F_out.dimethyl_carbonate.mol;

    % F_out.kta
    plot_struct.data.F_out.carbon_dioxide.kta(plot_row.row_number, 1) = plot_row.F_out.carbon_dioxide.kta;
    plot_struct.data.F_out.ethylene_oxide.kta(plot_row.row_number, 1) = plot_row.F_out.ethylene_oxide.kta;
    plot_struct.data.F_out.methanol.kta(plot_row.row_number, 1) = plot_row.F_out.methanol.kta;
    plot_struct.data.F_out.ethylene_carbonate.kta(plot_row.row_number, 1) = plot_row.F_out.ethylene_carbonate.kta;
    plot_struct.data.F_out.ethylene_glycol.kta(plot_row.row_number, 1) = plot_row.F_out.ethylene_glycol.kta;
    plot_struct.data.F_out.methoxy_ethanol.kta(plot_row.row_number, 1) = plot_row.F_out.methoxy_ethanol.kta;
    plot_struct.data.F_out.dimethyl_carbonate.kta(plot_row.row_number, 1) = plot_row.F_out.dimethyl_carbonate.kta;

    % R.mol
    plot_struct.data.R.carbon_dioxide.mol(plot_row.row_number, 1) = plot_row.R.carbon_dioxide.mol;
    plot_struct.data.R.ethylene_oxide.mol(plot_row.row_number, 1) = plot_row.R.ethylene_oxide.mol;
    plot_struct.data.R.methanol.mol(plot_row.row_number, 1) = plot_row.R.methanol.mol;
    plot_struct.data.R.ethylene_carbonate.mol(plot_row.row_number, 1) = plot_row.R.ethylene_carbonate.mol;
    plot_struct.data.R.ethylene_glycol.mol(plot_row.row_number, 1) = plot_row.R.ethylene_glycol.mol;
    plot_struct.data.R.methoxy_ethanol.mol(plot_row.row_number, 1) = plot_row.R.methoxy_ethanol.mol;
    plot_struct.data.R.dimethyl_carbonate.mol(plot_row.row_number, 1) = plot_row.R.dimethyl_carbonate.mol;

    % R.kta
    plot_struct.data.R.carbon_dioxide.kta(plot_row.row_number, 1) = plot_row.R.carbon_dioxide.kta;
    plot_struct.data.R.ethylene_oxide.kta(plot_row.row_number, 1) = plot_row.R.ethylene_oxide.kta;
    plot_struct.data.R.methanol.kta(plot_row.row_number, 1) = plot_row.R.methanol.kta;
    plot_struct.data.R.ethylene_carbonate.kta(plot_row.row_number, 1) = plot_row.R.ethylene_carbonate.kta;
    plot_struct.data.R.ethylene_glycol.kta(plot_row.row_number, 1) = plot_row.R.ethylene_glycol.kta;
    plot_struct.data.R.methoxy_ethanol.kta(plot_row.row_number, 1) = plot_row.R.methoxy_ethanol.kta;
    plot_struct.data.R.dimethyl_carbonate.kta(plot_row.row_number, 1) = plot_row.R.dimethyl_carbonate.kta;

    % reactor data
    plot_struct.data.conversion(plot_row.row_number) = plot_row.conversion;
    plot_struct.data.tau(plot_row.row_number) = plot_row.tau;
    plot_struct.data.V_rxtr(plot_row.row_number) = plot_row.V_rxtr;
     
end

function plot_struct = get_plot_structure(T, P, opt)
    user = get_user_inputs();
    % Thermodynamic state
    plot_struct.T = T;
    plot_struct.P = P;
    plot_struct.opt = opt;

    % Flowstreams
    % F_fresh.mol
    plot_struct.data.F_fresh.carbon_dioxide.mol = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_fresh.carbon_dioxide.mol = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_fresh.ethylene_oxide.mol = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_fresh.methanol.mol = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_fresh.ethylene_carbonate.mol = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_fresh.ethylene_glycol.mol = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_fresh.methoxy_ethanol.mol = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_fresh.dimethyl_carbonate.mol = get_empty_vector(length(user.level3.tau_range));
   
    % F_fresh.kta
    plot_struct.data.F_fresh.carbon_dioxide.kta = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_fresh.ethylene_oxide.kta = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_fresh.methanol.kta = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_fresh.ethylene_carbonate.kta = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_fresh.ethylene_glycol.kta = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_fresh.methoxy_ethanol.kta = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_fresh.dimethyl_carbonate.kta = get_empty_vector(length(user.level3.tau_range));

    % F_rxtr.mol
    plot_struct.data.F_rxtr.carbon_dioxide.mol = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_rxtr.ethylene_oxide.mol = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_rxtr.methanol.mol = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_rxtr.ethylene_carbonate.mol = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_rxtr.ethylene_glycol.mol = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_rxtr.methoxy_ethanol.mol = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_rxtr.dimethyl_carbonate.mol = get_empty_vector(length(user.level3.tau_range));

    % F_rxtr.kta
    plot_struct.data.F_rxtr.carbon_dioxide.kta = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_rxtr.ethylene_oxide.kta = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_rxtr.methanol.kta = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_rxtr.ethylene_carbonate.kta = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_rxtr.ethylene_glycol.kta = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_rxtr.methoxy_ethanol.kta = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_rxtr.dimethyl_carbonate.kta = get_empty_vector(length(user.level3.tau_range));

    % F_out.mol
    plot_struct.data.F_out.carbon_dioxide.mol = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_out.ethylene_oxide.mol = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_out.methanol.mol = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_out.ethylene_carbonate.mol = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_out.ethylene_glycol.mol = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_out.methoxy_ethanol.mol = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_out.dimethyl_carbonate.mol = get_empty_vector(length(user.level3.tau_range));

    % F_out.kta
    plot_struct.data.F_out.carbon_dioxide.kta = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_out.ethylene_oxide.kta = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_out.methanol.kta = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_out.ethylene_carbonate.kta = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_out.ethylene_glycol.kta = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_out.methoxy_ethanol.kta = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.F_out.dimethyl_carbonate.kta = get_empty_vector(length(user.level3.tau_range));

    % R.mol
    plot_struct.data.R.carbon_dioxide.mol = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.R.ethylene_oxide.mol = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.R.methanol.mol = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.R.ethylene_carbonate.mol = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.R.ethylene_glycol.mol = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.R.methoxy_ethanol.mol = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.R.dimethyl_carbonate.mol = get_empty_vector(length(user.level3.tau_range));

    % R.kta
    plot_struct.data.R.carbon_dioxide.kta = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.R.ethylene_oxide.kta = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.R.methanol.kta = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.R.ethylene_carbonate.kta = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.R.ethylene_glycol.kta = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.R.methoxy_ethanol.kta = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.R.dimethyl_carbonate.kta = get_empty_vector(length(user.level3.tau_range));

    % Reactor Data 
    plot_struct.data.conversion = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.tau = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.V_rxtr = get_empty_vector(length(user.level3.tau_range));

    % Flowstreams
    % F_fresh.mol
    plot_struct.data.F_fresh.carbon_dioxide.mol(:) = NaN;
    plot_struct.data.F_fresh.ethylene_oxide.mol(:) = NaN;
    plot_struct.data.F_fresh.methanol.mol(:) = NaN;
    plot_struct.data.F_fresh.ethylene_carbonate.mol(:) = NaN;
    plot_struct.data.F_fresh.ethylene_glycol.mol(:) = NaN;
    plot_struct.data.F_fresh.methoxy_ethanol.mol(:) = NaN;
    plot_struct.data.F_fresh.dimethyl_carbonate.mol(:) = NaN;

    % F_fresh.kta
    plot_struct.data.F_fresh.carbon_dioxide.kta(:) = NaN;
    plot_struct.data.F_fresh.ethylene_oxide.kta(:) = NaN;
    plot_struct.data.F_fresh.methanol.kta(:) = NaN;
    plot_struct.data.F_fresh.ethylene_carbonate.kta(:) = NaN;
    plot_struct.data.F_fresh.ethylene_glycol.kta(:) = NaN;
    plot_struct.data.F_fresh.methoxy_ethanol.kta(:) = NaN;
    plot_struct.data.F_fresh.dimethyl_carbonate.kta(:) = NaN;

    % F_rxtr.mol
    plot_struct.data.F_rxtr.carbon_dioxide.mol(:) = NaN;
    plot_struct.data.F_rxtr.ethylene_oxide.mol(:) = NaN;
    plot_struct.data.F_rxtr.methanol.mol(:) = NaN;
    plot_struct.data.F_rxtr.ethylene_carbonate.mol(:) = NaN;
    plot_struct.data.F_rxtr.ethylene_glycol.mol(:) = NaN;
    plot_struct.data.F_rxtr.methoxy_ethanol.mol(:) = NaN;
    plot_struct.data.F_rxtr.dimethyl_carbonate.mol(:) = NaN;

    % F_rxtr.kta
    plot_struct.data.F_rxtr.carbon_dioxide.kta(:) = NaN;
    plot_struct.data.F_rxtr.ethylene_oxide.kta(:) = NaN;
    plot_struct.data.F_rxtr.methanol.kta(:) = NaN;
    plot_struct.data.F_rxtr.ethylene_carbonate.kta(:) = NaN;
    plot_struct.data.F_rxtr.ethylene_glycol.kta(:) = NaN;
    plot_struct.data.F_rxtr.methoxy_ethanol.kta(:) = NaN;
    plot_struct.data.F_rxtr.dimethyl_carbonate.kta(:) = NaN;

    % F_out.mol
    plot_struct.data.F_out.carbon_dioxide.mol(:) = NaN;
    plot_struct.data.F_out.ethylene_oxide.mol(:) = NaN;
    plot_struct.data.F_out.methanol.mol(:) = NaN;
    plot_struct.data.F_out.ethylene_carbonate.mol(:) = NaN;
    plot_struct.data.F_out.ethylene_glycol.mol(:) = NaN;
    plot_struct.data.F_out.methoxy_ethanol.mol(:) = NaN;
    plot_struct.data.F_out.dimethyl_carbonate.mol(:) = NaN;

    % F_out.kta
    plot_struct.data.F_out.carbon_dioxide.kta(:) = NaN;
    plot_struct.data.F_out.ethylene_oxide.kta(:) = NaN;
    plot_struct.data.F_out.methanol.kta(:) = NaN;
    plot_struct.data.F_out.ethylene_carbonate.kta(:) = NaN;
    plot_struct.data.F_out.ethylene_glycol.kta(:) = NaN;
    plot_struct.data.F_out.methoxy_ethanol.kta(:) = NaN;
    plot_struct.data.F_out.dimethyl_carbonate.kta(:) = NaN;

    % R.mol
    plot_struct.data.R.carbon_dioxide.mol(:) = NaN;
    plot_struct.data.R.ethylene_oxide.mol(:) = NaN;
    plot_struct.data.R.methanol.mol(:) = NaN;
    plot_struct.data.R.ethylene_carbonate.mol(:) = NaN;
    plot_struct.data.R.ethylene_glycol.mol(:) = NaN;
    plot_struct.data.R.methoxy_ethanol.mol(:) = NaN;
    plot_struct.data.R.dimethyl_carbonate.mol(:) = NaN;

    % R.kta
    plot_struct.data.R.carbon_dioxide.kta(:) = NaN;
    plot_struct.data.R.ethylene_oxide.kta(:) = NaN;
    plot_struct.data.R.methanol.kta(:) = NaN;
    plot_struct.data.R.ethylene_carbonate.kta(:) = NaN;
    plot_struct.data.R.ethylene_glycol.kta(:) = NaN;
    plot_struct.data.R.methoxy_ethanol.kta(:) = NaN;
    plot_struct.data.R.dimethyl_carbonate.kta(:) = NaN;

    % Reactor Data 
    plot_struct.data.conversion(:) = NaN;
    plot_struct.data.tau(:) = NaN;
    plot_struct.data.V_rxtr(:) = NaN;

    % Plot info (maybe delete ??)
    plot_struct.name.title = '';
    plot_struct.name.y = '';
    plot_struct.name.x = '';
    % plot_struct.length = NaN;
    plot_struct.color = '';
end

% function plot_struct = get_plot_structure()
%     plot_struct.name.title = '';
%     plot_struct.name.y = '';
%     plot_struct.name.x = '';
%     plot_struct.length = NaN;
%     plot_struct.color = '';
% end

function vector = get_empty_vector(length)
    vector = zeros(length, 1);

end
