


function fxns = plot_fxns()

    fxns.get_empty_vector = @get_empty_vector;
    fxns.get_plot_struct = @get_plot_structure;
    fxns.set_plot_row = @set_plot_row;
    fxns.plot_reactor_volume_conversion = @plot_reactor_volume_conversion;
    fxns.plot_reactor_volume_conversion_allT = @plot_reactor_volume_conversion_allT;
    fxns.plot_reactor_volume_conversion_allP = @plot_reactor_volume_conversion_allP;
    fxns.delete_old_plots = @delete_old_plots;
    fxns.plot_molar_flowrates_conversion = @plot_molar_flowrates_conversion;
    fxns.plot_fresh_feed_conversion = @plot_fresh_feed_conversion;
    fxns.plot_recycle_flowrates = @plot_recycle_flowrates;
    fxns.plot_effluent_composition = @plot_effluent_composition;
    fxns.plot_total_reactor_feed = @plot_total_reactor_feed;
    fxns.plot_total_separation_feed = @plot_total_separation_feed;
    fxns.plot_selectivity = @plot_selectivity;
    fxns.plot_npv = @plot_npv

end

function void = plot_npv(plot_struct)
    void = NaN;
    user = get_user_inputs();
    const = get_constants();
    flow_fxns = flowrate_fxns();
    figure 
    hold on
    % F_total = flow_fxns.get_total_flowrate(plot_struct.data.F_out, 'mol');
    % fieldNames = fieldnames(plot_struct.data.F_rxtr);
    % for i = 1:length(fieldNames)
        x = plot_struct.data.conversion(:);
        % y = plot_struct.data.F_out.(fieldNames{i}).x(:);
        y = plot_struct.data.npv(:) / 15;

        % figure
        plot(x, y);
        title(sprintf('Nominal Present Value at [ %3.0f Bar ] [ %3.0f °C ]',plot_struct.P, plot_struct.T), 'Interpreter', 'tex');
        xlabel('\chi', 'Interpreter', 'tex');
        ylabel('NPV [ $MM ]', 'Interpreter', 'tex');
        % Create the legend entry for this plot
        % legendEntries{i} = sprintf('%s', strrep(fieldNames{i}, '_', " "));
    % end
    % The design variable point
    % legendEntries{i + 1} = sprintf('\\chi = %0.2f, %3.1f m^3' , user.plot.isothermal.x_point, ...
                                % (user.plot.isothermal.y_point * const.units.volume.m3_per_l));
    xline(user.plot.isothermal.x_point, '--k', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    % yline(user.plot.isothermal.y_point * const.units.volume.m3_per_l, '--k', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    % plot(user.plot.isothermal.x_point, user.plot.isothermal.y_point * const.units.volume.m3_per_l, ...
    %             'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize', 3);

    % Add Legend
    % legend(legendEntries, 'Interpreter', 'tex', 'location', 'west');
    hold off

    fileName =  "isothermal_npv_"+ string(plot_struct.P) + "Bar";
    save_plot(fileName);
end


function void = plot_selectivity(plot_struct)
    void = NaN;
    user = get_user_inputs();
    const = get_constants();
    flow_fxns = flowrate_fxns();
    figure 
    hold on
    % F_total = flow_fxns.get_total_flowrate(plot_struct.data.F_out, 'mol');
    S = plot_struct.data.F_out.dimethyl_carbonate.mol(:) ./ plot_struct.data.F_fresh.ethylene_oxide.mol(:);
        x = plot_struct.data.conversion(:);
        % y = plot_struct.data.F_out.(fieldNames{i}).x(:);
        y = S(:);

        % figure
        plot(x, y);
        title(sprintf('Selectivity at [ %3.0f Bar ] [ %3.0f °C ]',plot_struct.P, plot_struct.T), 'Interpreter', 'tex');
        xlabel('\chi', 'Interpreter', 'tex');
        ylabel('Selectivity', 'Interpreter', 'tex');
        % Create the legend entry for this plot
        % legendEntries{i} = sprintf('%s', strrep(fieldNames{i}, '_', " "));
    % end
    % The design variable point
    % legendEntries{i + 1} = sprintf('\\chi = %0.2f, %3.1f m^3' , user.plot.isothermal.x_point, ...
                                % (user.plot.isothermal.y_point * const.units.volume.m3_per_l));
    xline(user.plot.isothermal.x_point, '--k', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    % yline(user.plot.isothermal.y_point * const.units.volume.m3_per_l, '--k', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    % plot(user.plot.isothermal.x_point, user.plot.isothermal.y_point * const.units.volume.m3_per_l, ...
    %             'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize', 3);

    % Add Legend
    % legend(legendEntries, 'Interpreter', 'tex', 'location', 'west');
    hold off

    fileName =  "isothermal_selectivity" + string(plot_struct.P) + "Bar";
    save_plot(fileName);


end

function void = plot_total_separation_feed(plot_struct)
    void = NaN;
    user = get_user_inputs();
    const = get_constants();
    flow_fxns = flowrate_fxns();
    figure 
    hold on
    F_total = flow_fxns.get_total_flowrate(plot_struct.data.F_out, 'mol');
    % fieldNames = fieldnames(plot_struct.data.F_rxtr);
    % for i = 1:length(fieldNames)
        x = plot_struct.data.conversion(:);
        % y = plot_struct.data.F_out.(fieldNames{i}).x(:);
        y = F_total(:);

        % figure
        plot(x, y);
        title(sprintf('F_{Total Separator Feed} [ mol / s ] at [ %3.0f Bar ] [ %3.0f °C ]',plot_struct.P, plot_struct.T), 'Interpreter', 'tex');
        xlabel('\chi', 'Interpreter', 'tex');
        ylabel('F_{Total Separator Feed} [ mol / s ] ', 'Interpreter', 'tex');
        % Create the legend entry for this plot
        % legendEntries{i} = sprintf('%s', strrep(fieldNames{i}, '_', " "));
    % end
    % The design variable point
    % legendEntries{i + 1} = sprintf('\\chi = %0.2f, %3.1f m^3' , user.plot.isothermal.x_point, ...
                                % (user.plot.isothermal.y_point * const.units.volume.m3_per_l));
    xline(user.plot.isothermal.x_point, '--k', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    % yline(user.plot.isothermal.y_point * const.units.volume.m3_per_l, '--k', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    % plot(user.plot.isothermal.x_point, user.plot.isothermal.y_point * const.units.volume.m3_per_l, ...
    %             'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize', 3);

    % Add Legend
    % legend(legendEntries, 'Interpreter', 'tex', 'location', 'west');
    hold off

    fileName =  "isothermal_F_separator_total" + string(plot_struct.P) + "Bar";
    save_plot(fileName);
end

function void = plot_total_reactor_feed(plot_struct)
    void = NaN;
    user = get_user_inputs();
    const = get_constants();
    flow_fxns = flowrate_fxns();
    figure 
    hold on
    F_total = flow_fxns.get_total_flowrate(plot_struct.data.F_rxtr, 'mol');
    % fieldNames = fieldnames(plot_struct.data.F_rxtr);
    % for i = 1:length(fieldNames)
        x = plot_struct.data.conversion(:);
        % y = plot_struct.data.F_out.(fieldNames{i}).x(:);
        y = F_total(:);

        % figure
        plot(x, y);
        title(sprintf('F_{Total Reactor Feed} [ mol / s ] at [ %3.0f Bar ] [ %3.0f °C ]',plot_struct.P, plot_struct.T), 'Interpreter', 'tex');
        xlabel('\chi', 'Interpreter', 'tex');
        ylabel('F_{Total Reactor Feed} [ mol / s ] ', 'Interpreter', 'tex');
        % Create the legend entry for this plot
        % legendEntries{i} = sprintf('%s', strrep(fieldNames{i}, '_', " "));
    % end
    % The design variable point
    % legendEntries{i + 1} = sprintf('\\chi = %0.2f, %3.1f m^3' , user.plot.isothermal.x_point, ...
                                % (user.plot.isothermal.y_point * const.units.volume.m3_per_l));
    xline(user.plot.isothermal.x_point, '--k', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    % yline(user.plot.isothermal.y_point * const.units.volume.m3_per_l, '--k', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    % plot(user.plot.isothermal.x_point, user.plot.isothermal.y_point * const.units.volume.m3_per_l, ...
    %             'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize', 3);

    % Add Legend
    % legend(legendEntries, 'Interpreter', 'tex', 'location', 'west');
    hold off

    fileName =  "isothermal_F_rxtr_total" + string(plot_struct.P) + "Bar";
    save_plot(fileName);
end


function void = plot_effluent_composition(plot_struct)
    void = NaN;
    user = get_user_inputs();
    const = get_constants();
    flow_fxns = flowrate_fxns();
    figure 
    hold on
    fieldNames = fieldnames(plot_struct.data.F_out);
    % plot_struct.data.F_out = flow_fxns.set_mol_fractions(plot_struct.data.F_out);
    for i = 1:length(fieldNames)
        x = plot_struct.data.conversion(:);
        y = plot_struct.data.F_out.(fieldNames{i}).x(:);

        % figure
        plot(x, y);
        title(sprintf('X_i [ mol / mol ] at [ %3.0f Bar ] [ %3.0f °C ]',plot_struct.P, plot_struct.T), 'Interpreter', 'tex');
        xlabel('\chi', 'Interpreter', 'tex');
        ylabel('X_i [ mol / mol ]', 'Interpreter', 'tex')
        % Create the legend entry for this plot
        legendEntries{i} = sprintf('%s', strrep(fieldNames{i}, '_', " "));
    end
    % The design variable point
    legendEntries{i + 1} = sprintf('\\chi = %0.2f, %3.1f m^3' , user.plot.isothermal.x_point, ...
                                (user.plot.isothermal.y_point * const.units.volume.m3_per_l));
    xline(user.plot.isothermal.x_point, '--k', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    % yline(user.plot.isothermal.y_point * const.units.volume.m3_per_l, '--k', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    % plot(user.plot.isothermal.x_point, user.plot.isothermal.y_point * const.units.volume.m3_per_l, ...
    %             'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize', 3);

    % Add Legend
    legend(legendEntries, 'Interpreter', 'tex', 'location', 'west');
    hold off

    fileName =  "isothermal_effluent_composition_" + string(plot_struct.P) + "Bar";
    save_plot(fileName);
end

function void = plot_recycle_flowrates(plot_struct)
    void = NaN;
    user = get_user_inputs();
    const = get_constants();
    figure 
    hold on
    fieldNames = fieldnames(plot_struct.data.F_out);
    for i = 1:length(fieldNames)
        x = plot_struct.data.conversion(:);
        y = plot_struct.data.R.(fieldNames{i}).mol(:);

        % figure
        plot(x, y);
        title(sprintf('R [ mol / s ] at [ %3.0f Bar ] [ %3.0f °C ]',plot_struct.P, plot_struct.T), 'Interpreter', 'tex');
        xlabel('\chi', 'Interpreter', 'tex');
        ylabel('R [ mol / s ]', 'Interpreter', 'tex')
        % Create the legend entry for this plot
        legendEntries{i} = sprintf('%s', strrep(fieldNames{i}, '_', " "));
    end
    % The design variable point
    legendEntries{i + 1} = sprintf('\\chi = %0.2f, %3.1f m^3' , user.plot.isothermal.x_point, ...
                                (user.plot.isothermal.y_point * const.units.volume.m3_per_l));
    xline(user.plot.isothermal.x_point, '--k', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    % yline(user.plot.isothermal.y_point * const.units.volume.m3_per_l, '--k', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    % plot(user.plot.isothermal.x_point, user.plot.isothermal.y_point * const.units.volume.m3_per_l, ...
    %             'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize', 3);

    % Add Legend
    legend(legendEntries, 'Interpreter', 'tex', 'location', 'west');
    hold off

    fileName =  "isothermal_recycle_" + string(plot_struct.P) + "Bar";
    save_plot(fileName);
end

function void = plot_fresh_feed_conversion(plot_struct)
    void = NaN;
    user = get_user_inputs();
    const = get_constants();
    figure 
    hold on
    fieldNames = fieldnames(plot_struct.data.F_out);
    for i = 1:length(fieldNames)
        x = plot_struct.data.conversion(:);
        y = plot_struct.data.F_fresh.(fieldNames{i}).mol(:);

        % figure
        plot(x, y);
        title(sprintf('F_{fresh feed} [ mol / s ] at [ %3.0f Bar ] [ %3.0f °C ]',plot_struct.P, plot_struct.T), 'Interpreter', 'tex');
        xlabel('\chi', 'Interpreter', 'tex');
        ylabel('F_{fresh feed} [ mol / s ]', 'Interpreter', 'tex')
        % Create the legend entry for this plot
        legendEntries{i} = sprintf('%s', strrep(fieldNames{i}, '_', " "));
    end
    % The design variable point
    legendEntries{i + 1} = sprintf('\\chi = %0.2f, %3.1f m^3' , user.plot.isothermal.x_point, ...
                                (user.plot.isothermal.y_point * const.units.volume.m3_per_l));
    xline(user.plot.isothermal.x_point, '--k', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    % yline(user.plot.isothermal.y_point * const.units.volume.m3_per_l, '--k', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    % plot(user.plot.isothermal.x_point, user.plot.isothermal.y_point * const.units.volume.m3_per_l, ...
    %             'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize', 3);

    % Add Legend
    legend(legendEntries, 'Interpreter', 'tex', 'location', 'west');
    hold off

    fileName =  "isothermal_F_fresh_feed_" + string(plot_struct.P) + "Bar";
    save_plot(fileName);
end

function void = plot_molar_flowrates_conversion(plot_struct)
    void = NaN;
    user = get_user_inputs();
    const = get_constants();
    figure 
    hold on
    fieldNames = fieldnames(plot_struct.data.F_out);
    for i = 1:length(fieldNames)
        x = plot_struct.data.conversion(:);
        y = plot_struct.data.F_out.(fieldNames{i}).mol(:);

        % figure
        plot(x, y);
        title(sprintf('F_{effluent} [ mol / s ] at [ %3.0f Bar ] [ %3.0f °C ]',plot_struct.P, plot_struct.T), 'Interpreter', 'tex');
        xlabel('\chi', 'Interpreter', 'tex');
        ylabel('F_{effluent} [ mol / s ]', 'Interpreter', 'tex')
        % Create the legend entry for this plot
        legendEntries{i} = sprintf('%s', strrep(fieldNames{i}, '_', " "));
    end
    % The design variable point
    legendEntries{i + 1} = sprintf('\\chi = %0.2f, %3.1f m^3' , user.plot.isothermal.x_point, ...
                                (user.plot.isothermal.y_point * const.units.volume.m3_per_l));
    xline(user.plot.isothermal.x_point, '--k', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    % yline(user.plot.isothermal.y_point * const.units.volume.m3_per_l, '--k', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    % plot(user.plot.isothermal.x_point, user.plot.isothermal.y_point * const.units.volume.m3_per_l, ...
    %             'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize', 3);

    % Add Legend
    legend(legendEntries, 'Interpreter', 'tex', 'location', 'west');
    hold off

    fileName =  "isothermal_F_effluent_" + string(plot_struct.P) + "Bar";
    save_plot(fileName);
end

function void = save_plot(fileName)
    user = get_user_inputs();
    if ispc
        dir = [pwd  '\plots\' ];
    elseif ismac
        dir = [pwd  '/plots/' ];
    else
        dir = [pwd  '/plots/' ];
    end

    if ~exist(dir, 'dir')
        mkdir(dir);
    end

    delete_png_if_exists(fileName);
    print(fullfile(dir, fileName), '-dpng', user.plot.image_dpi) ;  % Save as PNG with 300 DPI
end

function void = delete_png_if_exists(fileName)
    % Check platform and set directory
    if ispc
        dirName = [pwd  '\plots\' ];
    elseif ismac
        dirName = [pwd  '/plots/' ];
    else
        dirName = [pwd  '/plots/' ];
    end

    % Construct the full file path
    filePath = fullfile(dirName, fileName);

    % Check if the file exists and delete if it does
    if exist(filePath, 'file') == 2 % 2 indicates it is a file
        delete(filePath);
    end
end

function void = delete_old_plots()
    % disp("test ")
    if ispc
        dirName = [pwd  '\plots\' ];
    elseif ismac
        dirName = [pwd  '/plots/' ];
    else
        dirName = [pwd  '/plots/' ];
    end

    pngFiles = dir(fullfile(dirName, '*.png'));

    for k = 1:length(pngFiles)
        filePath = fullfile(pngFiles(k).folder, pngFiles(k).name);
        delete(filePath);
    end

end

function void = plot_reactor_volume_conversion_allP(all_pressure_data)
    user = get_user_inputs();
    const = get_constants();
    figure 
    hold on
    for i = 1:length(all_pressure_data)
        plot_struct = all_pressure_data(i);
        x = plot_struct.data.conversion(:);
        y = plot_struct.data.V_rxtr(:) .* const.units.volume.m3_per_l;

        % figure
        plot(x, y);
        title(sprintf('V_{reactor} [ m^3 ] at [ %3.0f °C ]', plot_struct.T), 'Interpreter', 'tex');
        xlabel('\chi', 'Interpreter', 'tex');
        ylabel('V_{reactor} [ m^3 ]', 'Interpreter', 'tex')
        % Create the legend entry for this plot
        legendEntries{i} = sprintf('%3.0f Bar', plot_struct.P);
    end
    % The design variable point
    legendEntries{i + 1} = sprintf('\\chi = %0.2f, %3.1f m^3' , user.plot.isothermal.x_point, ...
                                (user.plot.isothermal.y_point * const.units.volume.m3_per_l));
    xline(user.plot.isothermal.x_point, '--k', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    yline(user.plot.isothermal.y_point * const.units.volume.m3_per_l, '--k', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    plot(user.plot.isothermal.x_point, user.plot.isothermal.y_point * const.units.volume.m3_per_l, ...
                'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize', 3);

    % Add Legend
    legend(legendEntries, 'Interpreter', 'tex', 'location', 'northwest');
    hold off

    fileName =  "isothermal_V_reactor";
    save_plot(fileName);
end 


function void = plot_reactor_volume_conversion_allT(all_temp_data)
    user = get_user_inputs();
    const = get_constants();
    figure 
    hold on
    for i = 1:length(all_temp_data)
        plot_struct = all_temp_data(i);
        x = plot_struct.data.conversion(:) ;
        y = plot_struct.data.V_rxtr(:) * const.units.volume.m3_per_l;

        % figure
        plot(x, y);
        title(sprintf('V_{reactor} [ m^3 ] [ %3.0f Bar ]', plot_struct.P), 'Interpreter', 'tex');
        xlabel('\chi', 'Interpreter', 'tex');
        ylabel('V_{reactor} [ m^3 ]', 'Interpreter', 'tex')
        % Create the legend entry for this plot
        legendEntries{i} = sprintf('%3.0f°C', plot_struct.T);
    end
    legend(legendEntries, 'Interpreter', 'tex', 'location', 'northwest');
    hold off

    % if ispc
    %     dir = [pwd  '\plots\' ];
    % elseif ismac
    %     dir = [pwd  '/plots/' ];
    % else
    %     dir = [pwd  '/plots/' ];
    % end

    % if ~exist(dir, 'dir')
    %     mkdir(dir);
    % end

    % print(fullfile(dir, "isobaric_V_reactor"), '-dpng', user.plot.image_dpi) ;  % Save as PNG with 300 DPI

    fileName =   "isobaric_V_reactor";
    save_plot(fileName);
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

    fields = {'carbon_dioxide', 'ethylene_oxide', 'methanol', 'ethylene_carbonate', 'ethylene_glycol', 'methoxy_ethanol', 'dimethyl_carbonate'};
    subfields = {'mol', 'kta', 'x'};

    % Iterate over F_fresh
    for i = 1:length(fields)
        for j = 1:length(subfields)
            plot_struct.data.F_fresh.(fields{i}).(subfields{j})(plot_row.row_number, 1) = plot_row.F_fresh.(fields{i}).(subfields{j});
        end
    end

    % Iterate over F_rxtr
    for i = 1:length(fields)
        for j = 1:length(subfields)
            plot_struct.data.F_rxtr.(fields{i}).(subfields{j})(plot_row.row_number, 1) = plot_row.F_rxtr.(fields{i}).(subfields{j});
        end
    end

    % Iterate over F_out
    for i = 1:length(fields)
        for j = 1:length(subfields)
            plot_struct.data.F_out.(fields{i}).(subfields{j})(plot_row.row_number, 1) = plot_row.F_out.(fields{i}).(subfields{j});
        end
    end

    % Iterate over R
    for i = 1:length(fields)
        for j = 1:length(subfields)
            plot_struct.data.R.(fields{i}).(subfields{j})(plot_row.row_number, 1) = plot_row.R.(fields{i}).(subfields{j});
        end
    end

    % reactor data
    plot_struct.data.conversion(plot_row.row_number) = plot_row.conversion;
    plot_struct.data.tau(plot_row.row_number) = plot_row.tau;
    plot_struct.data.V_rxtr(plot_row.row_number) = plot_row.V_rxtr;
    plot_struct.data.npv(plot_row.row_number) = plot_row.npv;
     
end

function plot_struct = get_plot_structure(T, P, opt)
    user = get_user_inputs();
    % Thermodynamic state
    plot_struct.T = T;
    plot_struct.P = P;
    plot_struct.opt = opt;

    fields = {'carbon_dioxide', 'ethylene_oxide', 'methanol', 'ethylene_carbonate', 'ethylene_glycol', 'methoxy_ethanol', 'dimethyl_carbonate'};
    subfields = {'mol', 'kta'};
    categories = {'F_fresh', 'F_rxtr', 'F_out', 'R'};

    % Initialize the fields with get_empty_vector
    for c = 1:length(categories)
        for i = 1:length(fields)
            for j = 1:length(subfields)
                plot_struct.data.(categories{c}).(fields{i}).(subfields{j}) = get_empty_vector(length(user.level3.tau_range));
            end
        end
    end

    % Initialize NaN values for the fields
    for c = 1:length(categories)
        for i = 1:length(fields)
            for j = 1:length(subfields)
                plot_struct.data.(categories{c}).(fields{i}).(subfields{j})(:) = NaN;
            end
        end
    end

    % Initialize reactor data
    plot_struct.data.conversion = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.tau = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.V_rxtr = get_empty_vector(length(user.level3.tau_range));
    plot_struct.data.npv = get_empty_vector(length(user.level3.tau_range));

    % Reactor Data 
    plot_struct.data.conversion(:) = NaN;
    plot_struct.data.tau(:) = NaN;
    plot_struct.data.V_rxtr(:) = NaN;
    plot_struct.data.npv(:) = NaN;

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
