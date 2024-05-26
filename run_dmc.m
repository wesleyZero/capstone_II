clc; clear; close all; 

% SCRIPT________________________________________________________________________

level3();


% FUNCTIONS_____________________________________________________________________

function void = level3()
    level3_isobaric();
    level3_isothermal();
    void = NaN;

end

function void = level3_isobaric()
    void = NaN;

    console = get_console();
    const = get_constants(); 
    user = get_user_inputs();
    F_fxns = flowrate_fxns();
    rxtr_fxns = reactor_fxns();
    plt_fxns = plot_fxns();
    P = user.level3.isobaric_press.bar;
    opt = 'isobaric';
    % test1 = opt + "test q"
    console.section("Starting Level 3 " + opt + " calculations")
    
    i = 1;
    for T = user.level3.temp_range
        console.subsection(sprintf("T = %3.2f", T), 1);
        row = 1;
        isoBar_plt = plt_fxns.get_plot_struct(T, P, opt);
        

        for tau = user.level3.tau_range
            console.subsection(sprintf("tau = %3.2f", tau), 2)
            [F_fresh, F_rxtr, F_out, R, V_rxtr] = level3_flowrates(tau, T, P, opt); 
            conversion = rxtr_fxns.get_conversion(F_rxtr, F_out);

            if isnan(conversion)
                disp("ERROR : COMPLEX CONC. BREAKING TO NEXT TEMP")
                break;
            end
            % Store row data
            plot_row.F_fresh = F_fresh;
            plot_row.F_rxtr = F_rxtr;
            plot_row.F_out = F_out;
            plot_row.R = R;
            plot_row.conversion = conversion;
            plot_row.row_number = row;
            plot_row.tau = tau;
            plot_row.V_rxtr = V_rxtr;
            isoBar_plt = plt_fxns.set_plot_row(isoBar_plt, plot_row);
            % increment
            row = row + 1; 
        end

        all_temp_data(i) = isoBar_plt;
        i = i + 1;
        % plt_fxns.plot_reactor_volume_conversion(isoBar_plt);
    end

    plt_fxns.plot_reactor_volume_conversion_allT(all_temp_data);
    
    console.section("Level 3 " + opt + " calculations are complete")
end

function void = level3_isothermal()
    void = NaN;

    console = get_console();
    const = get_constants(); 
    user = get_user_inputs();
    F_fxns = flowrate_fxns();
    rxtr_fxns = reactor_fxns();
    plt_fxns = plot_fxns();
    T = user.level3.isothermal_temp.C;
    opt = 'isothermal';
    % test1 = opt + "test q"
    console.section("Starting Level 3 " + opt + " calculations")
    
    for P = user.level3.press_range
        console.subsection(sprintf("P = %3.2f", P), 1);
        row = 1;
        isoBar_plt = plt_fxns.get_plot_struct(T, P, opt);
        

        for tau = user.level3.tau_range
            console.subsection(sprintf("tau = %3.2f", tau), 2)
            [F_fresh, F_rxtr, F_out, R, V_rxtr] = level3_flowrates(tau, T, P, opt); 
            conversion = rxtr_fxns.get_conversion(F_rxtr, F_out);

            if isnan(conversion)
                disp("ERROR : COMPLEX CONC. BREAKING TO NEXT TEMP")
                break;
            end
            % Store row data
            plot_row.F_fresh = F_fresh;
            plot_row.F_rxtr = F_rxtr;
            plot_row.F_out = F_out;
            plot_row.R = R;
            plot_row.conversion = conversion;
            plot_row.row_number = row;
            plot_row.tau = tau;
            plot_row.V_rxtr = V_rxtr;
            isoBar_plt = plt_fxns.set_plot_row(isoBar_plt, plot_row);
            % increment
            row = row + 1; 
        end

    end
    
    console.section("Level 3 " + opt + " calculations are complete")
end

function [F_fresh, F_rxtr, F_out, R, V_rxtr] = level3_flowrates(tau, temp, P, opt)
    F_fresh = NaN; F_rxtr = NaN; F_out = NaN; R = NaN;
    user = get_user_inputs(); 
    flow_fxns = flowrate_fxns();
    rxtr_fxns = reactor_fxns();

    F_basis = flow_fxns.get_basis_feed_flowrates();
    [F_fresh, F_rxtr, F_out, R, V_rxtr] = rxtr_fxns.get_reactor_flows(F_basis, temp, P, opt, tau);
    

end


