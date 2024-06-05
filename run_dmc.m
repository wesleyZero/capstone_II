clc; clear; close all; 

% SCRIPT________________________________________________________________________
level3() 

% FUNCTIONS_____________________________________________________________________

function void = level3()
    void = NaN;
    P = 74; % bar
    plt_fxns = plot_fxns();
    level3_isothermal_aspen_compare(P);
    level3_isothermal(P);
    level3_isobaric();
    level3_isothermal(NaN);

end



function void = level3_isothermal_aspen_compare(P_specify)
    void = NaN;

    console = get_console();
    const = get_constants(); 
    user = get_user_inputs();
    F_fxns = flowrate_fxns();
    rxtr_fxns = reactor_fxns();
    plt_fxns = plot_fxns();
    econ_fxns = economic_fxns();
    T = user.level3.isothermal_temp.C;
    opt = 'isothermal';
    console.section("Starting Level 3 " + opt + " calculations")
    
    if ~isnan(P_specify)
        % i = 1;
        P = P_specify;
        % for P = user.level3.press_range
        console.subsection(sprintf("P = %3.2f", P), 1);
        row = 1;
        isoTherm_plt = plt_fxns.get_plot_struct(T, P, opt);

        for tau = user.level3.tau_range.P_specify.isothermal
            console.subsection(sprintf("tau = %3.2f", tau), 2)
            [F_fresh, F_rxtr, F_out, R, V_rxtr] = level3_flowrates(tau, T, P, opt); 
            conversion = rxtr_fxns.get_conversion(F_rxtr, F_out);

            if isnan(conversion)
                disp("ERROR : COMPLEX CONC. BREAKING TO NEXT TEMP")
                break;
            end

            npv_opt = 'matlab';
            % npv_opt = 'matlab';
            npv = econ_fxns.get_work_min_npv(F_out, T, P, V_rxtr, conversion, npv_opt);

            % Store row data
            plot_row.F_fresh = F_fresh;
            plot_row.F_rxtr = F_rxtr;
            plot_row.F_out = F_out;
            plot_row.R = R;
            plot_row.conversion = conversion;
            plot_row.row_number = row;
            plot_row.tau = tau;
            plot_row.V_rxtr = V_rxtr ;
            plot_row.npv = npv;
            isoTherm_plt = plt_fxns.set_plot_row(isoTherm_plt, plot_row);
            % increment
            row = row + 1; 
        end

        % all_pressure_data(i) = isoTherm_plt;
        % i = i + 1;
        % plt_fxns.plot_selectivity(isoTherm_plt);
        % plt_fxns.plot_effluent_composition(isoTherm_plt);
        % plt_fxns.plot_total_separation_feed(isoTherm_plt);
        % plt_fxns.plot_total_reactor_feed(isoTherm_plt);
        % plt_fxns.plot_recycle_flowrates(isoTherm_plt);
        % plt_fxns.plot_fresh_feed_conversion(isoTherm_plt);
        % plt_fxns.plot_molar_flowrates_conversion(isoTherm_plt);
        % plt_fxns.plot_npv(isoTherm_plt);
        plt_fxns.plot_npv_with_aspen_data(isoTherm_plt);
    % else
    %     i = 1;
    %     for P = user.level3.press_range
    %         console.subsection(sprintf("P = %3.2f", P), 1);
    %         row = 1;
    %         isoTherm_plt = plt_fxns.get_plot_struct(T, P, opt);
            
    %             console.subsection(sprintf("tau = %3.2f", tau), 2)
    %             [F_fresh, F_rxtr, F_out, R, V_rxtr] = level3_flowrates(tau, T, P, opt); 
    %             conversion = rxtr_fxns.get_conversion(F_rxtr, F_out);

    %             if isnan(conversion)
    %                 disp("ERROR : COMPLEX CONC. BREAKING TO NEXT TEMP")
    %                 break;
    %             end
    %             npv_opt = 'aspen';
    %             npv = econ_fxns.get_work_min_npv(F_out, T, P, V_rxtr, conversion, npv_opt);
    
    %             % Store row data
    %             plot_row.F_fresh = F_fresh;
    %             plot_row.F_rxtr = F_rxtr;
    %             plot_row.F_out = F_out;
    %             plot_row.R = R;
    %             plot_row.conversion = conversion;
    %             plot_row.row_number = row;
    %             plot_row.tau = tau;
    %             plot_row.V_rxtr = V_rxtr ;
    %             plot_row.npv = npv;
    %             isoTherm_plt = plt_fxns.set_plot_row(isoTherm_plt, plot_row);
    %             % increment
    %             row = row + 1; 
    %         end

    %         all_pressure_data(i) = isoTherm_plt;
    %         i = i + 1;
    %     end
    %     plt_fxns.plot_reactor_volume_conversion_allP(all_pressure_data);
    end
    console.section("Level 3 " + opt + " calculations are complete")
end


function void = level3_isothermal(P_specify)
    void = NaN;

    console = get_console();
    const = get_constants(); 
    user = get_user_inputs();
    F_fxns = flowrate_fxns();
    rxtr_fxns = reactor_fxns();
    plt_fxns = plot_fxns();
    econ_fxns = economic_fxns();
    T = user.level3.isothermal_temp.C;
    opt = 'isothermal';
    console.section("Starting Level 3 " + opt + " calculations")
    
    if ~isnan(P_specify)
        i = 1;
        P = P_specify;
        % for P = user.level3.press_range
        console.subsection(sprintf("P = %3.2f", P), 1);
        row = 1;
        isoTherm_plt = plt_fxns.get_plot_struct(T, P, opt);

        for tau = user.level3.tau_range.P_specify.isothermal
            console.subsection(sprintf("tau = %3.2f", tau), 2)
            [F_fresh, F_rxtr, F_out, R, V_rxtr] = level3_flowrates(tau, T, P, opt); 
            conversion = rxtr_fxns.get_conversion(F_rxtr, F_out);

            if isnan(conversion)
                disp("ERROR : COMPLEX CONC. BREAKING TO NEXT TEMP")
                break;
            end

            % npv_opt = 'aspen';
            npv_opt = 'matlab';
            npv = econ_fxns.get_work_min_npv(F_out, T, P, V_rxtr, conversion, npv_opt);

            % Store row data
            plot_row.F_fresh = F_fresh;
            plot_row.F_rxtr = F_rxtr;
            plot_row.F_out = F_out;
            plot_row.R = R;
            plot_row.conversion = conversion;
            plot_row.row_number = row;
            plot_row.tau = tau;
            plot_row.V_rxtr = V_rxtr ;
            plot_row.npv = npv;
            isoTherm_plt = plt_fxns.set_plot_row(isoTherm_plt, plot_row);
            % increment
            row = row + 1; 
        end

        all_pressure_data(i) = isoTherm_plt;
        i = i + 1;
        % end
        plt_fxns.plot_selectivity(isoTherm_plt);
        plt_fxns.plot_effluent_composition(isoTherm_plt);
        plt_fxns.plot_total_separation_feed(isoTherm_plt);
        plt_fxns.plot_total_reactor_feed(isoTherm_plt);
        plt_fxns.plot_recycle_flowrates(isoTherm_plt);
        plt_fxns.plot_fresh_feed_conversion(isoTherm_plt);
        plt_fxns.plot_molar_flowrates_conversion(isoTherm_plt);
        plt_fxns.plot_npv(isoTherm_plt);
    else
        i = 1;
        for P = user.level3.press_range
            console.subsection(sprintf("P = %3.2f", P), 1);
            row = 1;
            isoTherm_plt = plt_fxns.get_plot_struct(T, P, opt);
            

            for tau = user.level3.tau_range.isothermal
                console.subsection(sprintf("tau = %3.2f", tau), 2)
                [F_fresh, F_rxtr, F_out, R, V_rxtr] = level3_flowrates(tau, T, P, opt); 
                conversion = rxtr_fxns.get_conversion(F_rxtr, F_out);

                if isnan(conversion)
                    disp("ERROR : COMPLEX CONC. BREAKING TO NEXT TEMP")
                    break;
                end
                npv_opt = 'matlab';
                npv = econ_fxns.get_work_min_npv(F_out, T, P, V_rxtr, conversion, npv_opt);
    
                % Store row data
                plot_row.F_fresh = F_fresh;
                plot_row.F_rxtr = F_rxtr;
                plot_row.F_out = F_out;
                plot_row.R = R;
                plot_row.conversion = conversion;
                plot_row.row_number = row;
                plot_row.tau = tau;
                plot_row.V_rxtr = V_rxtr ;
                plot_row.npv = npv;
                isoTherm_plt = plt_fxns.set_plot_row(isoTherm_plt, plot_row);
                % increment
                row = row + 1; 
            end

            all_pressure_data(i) = isoTherm_plt;
            i = i + 1;
        end
        plt_fxns.plot_reactor_volume_conversion_allP(all_pressure_data);
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



function void = level3_isobaric()
    void = NaN;

    console = get_console();
    const = get_constants(); 
    user = get_user_inputs();
    F_fxns = flowrate_fxns();
    rxtr_fxns = reactor_fxns();
    plt_fxns = plot_fxns();
    econ_fxns = economic_fxns();

    P = user.level3.isobaric_press.bar;
    opt = 'isobaric';
    console.section("Starting Level 3 " + opt + " calculations")
    
    i = 1;
    for T = user.level3.temp_range
        console.subsection(sprintf("T = %3.2f", T), 1);
        row = 1;
        isoBar_plt = plt_fxns.get_plot_struct(T, P, opt);
        

        for tau =  user.level3.tau_range.isobaric 
            console.subsection(sprintf("tau = %3.2f", tau), 2)
            [F_fresh, F_rxtr, F_out, R, V_rxtr] = level3_flowrates(tau, T, P, opt); 
            conversion = rxtr_fxns.get_conversion(F_rxtr, F_out);
			
            if isnan(conversion)
                disp("ERROR : COMPLEX CONC. BREAKING TO NEXT TEMP")
                break;
            end

            npv_opt = 'aspen';
            npv = econ_fxns.get_work_min_npv(F_out, T, P, V_rxtr, conversion, npv_opt);

            % Store row data
            plot_row.F_fresh = F_fresh;
            plot_row.F_rxtr = F_rxtr;
            plot_row.F_out = F_out;
            plot_row.R = R;
            plot_row.conversion = conversion;
            plot_row.row_number = row;
            plot_row.tau = tau;
            plot_row.V_rxtr = V_rxtr;
            plot_row.npv = npv;
            isoBar_plt = plt_fxns.set_plot_row(isoBar_plt, plot_row);
            % increment
            row = row + 1; 
        end

        all_temp_data(i) = isoBar_plt;
        i = i + 1;
    end

    plt_fxns.plot_reactor_volume_conversion_allT(all_temp_data);
    
    console.section("Level 3 " + opt + " calculations are complete")
end


