% HYSYS Costing Function
% Order of Components: CO2, EO, MeOH, DMC, ME, EG, EC

function [ISBL, OPEX, CO2e, feedCosts, productRev] = HYSYS_Cost(filename, econparam)

    % Extract values from econparam
    flowprices = econparam.flowprices;
    fuelprice = econparam.fuelprice;
    elec_price = econparam.elec_price;
    
    % Connect to HYSYS
    HYSYS = actxserver('Hysys.Application');
    simCase = HYSYS.SimulationCase.Open([cd,strcat('\',filename,'.hsc')]);
    simCase.Visible = true;
    
    % Import all streams and unit operations from flow sheet
    streams = simCase.Flowsheet.Streams;
    unitOps = simCase.Flowsheet.Operations;
    
    %% FEEDS AND PRODUCTS
    
    % Get access to fresh flow streams
    F_EO = get(streams, 'item', 'Fresh_EO');
    F_CO2 = get(streams, 'item', 'Fresh_CO2');
    F_MeOH = get(streams, 'item', 'Fresh_MeOH');
    
    % Determine fresh feed flow rates
    F = [F_CO2.MolarFlowValue, F_EO.MolarFlowValue,...
        F_MeOH.MolarFlowValue, 0, 0, 0, 0] * 3600; % [kmol/hr]
    
    % Get access to product streams
    P_DMC = get(streams, 'item', 'P_DMC');
    P_EG = get(streams, 'item', 'P_EG');
    P_ME = get(streams, 'item', 'D1_S');
    
    % Determine product flow rates
    P = [0, 0, 0, P_DMC.MolarFlowValue, P_ME.MolarFlowValue,...
        P_EG.MolarFlowValue, 0] * 3600; % [kmol/hr]
    
    % Get feed costs and product revenue
    feedCosts = flowCosts(F, flowprices);
    productRev = flowCosts(P, flowprices);
    
    %% REACTOR
    
    % Get access to reactor and inlet flow
    ToCSTR = get(streams, 'item', 'ToCSTR');
    CSTR = get(unitOps,'item','CSTR-100');
    
    % Get total heat duty of reactor
    Xheat = get(streams, 'item', 'HeatRxn');
    Rheat = get(streams, 'item', 'Reactor_E');
    heatduty = (Rheat.HeatFlowValue + Xheat.HeatFlowValue) * 3600; % [kJ/hr]
    
    % Cost vessel as jacketed CSTR
    rISBL = IC_pvessel(ToCSTR.PressureValue, CSTR.VolumeValue, "SS", false);
    rISBL = rISBL + IC_pvessel(100, CSTR.VolumeValue + 1, "CS", false);
    
    % Calculate reactor OPEX
    rOPEX = cost_cooler(-heatduty, econparam, ToCSTR.TemperatureValue);
    
    %% HEATERS AND COOLERS
    
    % Get all heaters and coolers
    heater = get(unitOps, 'item', 'E-100'); % Heater pre-reactor
    cooler = get(unitOps, 'item', 'E-103'); % Cooler post-reactor
    
    % Define vectors with all heater/cooler costs
    duty = [heater.DutyValue, -cooler.DutyValue]*3600; % [kJ/hr]
    T = [heater.ProductTemperatureValue, cooler.ProductTemperatureValue]; % [C]
    P = [heater.ProductPressureValue, cooler.ProductPressureValue]; % [kPa]
    material = ["CS/SS", "CS/SS"];
    
    % Initialize utilities ISBL, OPEX, and CO2e
    utilISBL = 0;
    utilOPEX = 0;
    utilCO2e = 0;
    
    % Calculate total utilities ISBL, OPEX, and CO2e
    for i = 1:length(duty)
    
        % Calculate individual ISBL, OPEX, and CO2e
        [uISBL, uOPEX, uCO2e] =...
            trimCosts(duty(i), T(i), P(i), material(i), "U", econparam);
    
        % Add to totals
        utilISBL = utilISBL + uISBL;
        utilOPEX = utilOPEX + uOPEX;
        utilCO2e = utilCO2e + uCO2e;
    
    end
    clear('uISBL', 'uOPEX', 'uCO2e')
    
    %% GAS COMPRESSOR
    
    K1 = get(unitOps, 'item', 'K-101');
    P_in = K1.FeedPressureValue;
    P_out = K1.ProductPressureValue;
    T_in = K1.FeedTemperatureValue;
    
    K1feed = get(streams, 'item', 'D1_DV');
    Q_in = K1feed.ActualVolumeFlowValue*3600; % [m^3/hr]
    cp = K1feed.MolarHeatCapacityValue;
    
    cISBL = sizeCompressor(P_in, P_out, T_in, Q_in, cp, "CM");
    
    % OPEX of compressor
    cOPEX = elec_price*K1.EnergyValue*3600*8400/power(10,6);
    
    %% DISTILLATION COLUMNS
    
    % Get access to all distillation columns
    Col1 = get(unitOps, 'item', 'DCol1');
    ECEGCol = get(unitOps, 'item', 'ECEG-Col');
    area_sheet = get(unitOps, 'item', 'MATLAB-Help');
    
    % Get access to units within distillation columns
    Col1T = get(Col1.ColumnFlowsheet.Operations, 'item', 'Main Tower');
    ECEGColT = get(ECEGCol.ColumnFlowsheet.Operations, 'item', 'Main Tower');
    
    Col1C = get(Col1.ColumnFlowsheet.Operations, 'item', 'Condenser');
    ECEGColC = get(ECEGCol.ColumnFlowsheet.Operations, 'item', 'Condenser');
    
    Col1R = get(Col1.ColumnFlowsheet.Operations, 'item', 'Reboiler');
    ECEGColR = get(ECEGCol.ColumnFlowsheet.Operations, 'item', 'Reboiler');
    
    % Get specified column values
    distP = [Col1C.VesselPressureValue, ECEGColC.VesselPressureValue]; % [kPa]
    distA = [area_sheet.Cell('B10').CellValue, area_sheet.Cell('C10').CellValue]; % [m^2]
    stages = [Col1T.NumberOfStages, ECEGColT.NumberOfStages];
    QC = [Col1C.HeatFlowValue, ECEGColC.HeatFlowValue]*3600; % [kJ/hr]
    TC = [Col1C.VesselTemperatureValue, ECEGColC.VesselTemperatureValue]; % [C]
    QRB = [Col1R.HeatFlowValue, ECEGColR.HeatFlowValue]*3600; % [kJ/hr]
    
    % Initialize distillation ISBL, OPEX, and CO2e
    distISBL = 0;
    distOPEX = 0;
    distCO2e = 0;
    
    % Calculate total distillation ISBL, OPEX, and CO2e
    for i = 1:length(distP)
    
        % Calculate individual ISBL, OPEX, and CO2e
        [dISBL, dOPEX, dCO2e] = distillation_costs(distP(i), distA(i), stages(i),...
            QC(i), TC(i), QRB(i), econparam);
    
        % Add to totals
        distISBL = distISBL + dISBL;
        distOPEX = distOPEX + dOPEX;
        distCO2e = distCO2e + dCO2e;
    
    end
    clear('dISBL', 'dOPEX', 'dCO2e')
    
    ISBL = rISBL + utilISBL + distISBL + cISBL;
    OPEX = rOPEX + utilOPEX + distOPEX + cOPEX;
    CO2e = utilCO2e + distCO2e;
    
    disp("ISBL: " + ISBL)
    disp("OPEX: " + OPEX)
    disp("CO2e: " + CO2e)

end % function






