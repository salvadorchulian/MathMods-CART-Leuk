function interactive_ODE_model
    clc; close all;
    
    %% Parameters
    params = struct( ...
        'Tend', 460, 'Yend', 1.5e10, 'eps', 0, 'delta', 0, ...
        'L0', 1e7, 'mrdPct', 0.1, 'scenario', 1, ... 
        'BStem', 1e8, 'B0', 1e7, 'CAR0', 1e7, 'CM0', 0, ...
        'alpha', 3e-10, 'k', 1e11, 'h', 1e10, 'rhoCA', 0.9, ...
        'Bmax', 5e9, 'rhoL', 1/30, 'rhoB', log(2)/8, ...
        'tauB', 60, 'tauCA', 6.5, 'tauCM', 300, ...
        'gammaAM', 0.001, 'gammaMA', 0.33, 'Lmax', 1e12);
    defaultParams = params; 

    %% Figure & Layout
    f = figure('Name', 'Interactive Bivariate ODE CAR-T Model', 'Position', [50 50 1400 900], 'Color', 'w');
    
    axTumor = axes('Parent', f, 'Position', [0.06 0.55 0.60 0.38]); 
    axCAR   = axes('Parent', f, 'Position', [0.06 0.08 0.60 0.38]);  
    
    % Control Panel setup
    panel = uipanel('Parent', f, 'Title', 'Controls', 'Units', 'normalized', 'Position', [0.72 0.08 0.26 0.85]);
    
    ypos = 0.94;
    
    % --- TOP BUTTONS (Reset & Export) ---
    uicontrol('Parent', panel, 'Style', 'pushbutton', 'String', 'Reset Controls', 'Units', 'normalized', ...
        'Position', [0.15 0.95 0.3 0.04], 'FontSize', 9, 'Callback', @resetSliders);
        
    uicontrol('Parent', panel, 'Style', 'pushbutton', 'String', 'Export (Data & Images)', 'Units', 'normalized', ...
        'Position', [0.55 0.95 0.35 0.04], 'FontSize', 9, 'Callback', @exportAll);
        
    ypos = ypos - 0.05;
    
    % --- SCENARIO DROPDOWN ---
    uicontrol('Parent', panel, 'Style', 'text', 'String', 'Tumor Initial State', 'Units', 'normalized', ...
        'Position', [0.02 ypos 0.25 0.04], 'FontSize', 9, 'HorizontalAlignment', 'right');
    ui_scenario = uicontrol('Parent', panel, 'Style', 'popupmenu', ...
        'String', {'100% Antigen-Positive', 'Mixed Clones', '100% Antigen-Negative'}, ...
        'Units', 'normalized', 'Position', [0.30 ypos 0.67 0.04], 'Value', params.scenario, 'Callback', @updatePlot);
    ypos = ypos - 0.055;
    
    % Separator Line
    annotation('line', [0.73 0.97], [ypos+0.04 ypos+0.04], 'Color', [0.8 0.8 0.8]);
    
    % --- SLIDERS ---
    sliders.mrdPct  = createSlider(panel, 'MRD Thresh (%)', 1e-4, 1, params.mrdPct, @updatePlot, ypos); ypos = ypos - 0.055;
    sliders.L0      = createSlider(panel, 'L0 (Total)', 0, 1e11, params.L0, @updatePlot, ypos); ypos = ypos - 0.055;
    sliders.BStem   = createSlider(panel, 'BStem', 0, 1e9, params.BStem, @updatePlot, ypos); ypos = ypos - 0.055;
    sliders.B0      = createSlider(panel, 'B0 (Initial)', 1e7, 1e9, params.B0, @updatePlot, ypos); ypos = ypos - 0.055;
    sliders.CAR0    = createSlider(panel, 'CAR0', 1e7, 1.5e7, params.CAR0, @updatePlot, ypos); ypos = ypos - 0.055;
    sliders.eps     = createSlider(panel, 'eps', 0, 1e-2, params.eps, @updatePlot, ypos); ypos = ypos - 0.055;
    sliders.delta   = createSlider(panel, 'delta', 0, 1e-20, params.delta, @updatePlot, ypos); ypos = ypos - 0.055;
    sliders.alpha   = createSlider(panel, 'alpha', 1e-11, 1e-9, params.alpha, @updatePlot, ypos); ypos = ypos - 0.055;
    sliders.Tend    = createSlider(panel, 'Tend', 30, 12000, params.Tend, @updatePlot, ypos); ypos = ypos - 0.055;
    sliders.Yend    = createSlider(panel, 'Yend Limit L', 1e6, params.Lmax, params.Yend, @updatePlot, ypos); ypos = ypos - 0.055;
    sliders.YendCAR    = createSlider(panel, 'Yend Limit CAR', 1e6, params.Lmax, params.Yend, @updatePlot, ypos); ypos = ypos - 0.055;
    sliders.k       = createSlider(panel, 'k', 1e8, 1e12, params.k, @updatePlot, ypos); ypos = ypos - 0.055;
    sliders.h       = createSlider(panel, 'h', 1e8, 1e12, params.h, @updatePlot, ypos); ypos = ypos - 0.055;
    sliders.rhoCA   = createSlider(panel, 'rhoCA', 0.01, 1, params.rhoCA, @updatePlot, ypos); ypos = ypos - 0.055;
    sliders.Bmax    = createSlider(panel, 'Bmax', 1e9, 2e10, params.Bmax, @updatePlot, ypos); 
    
    % Initial Plot
    updatePlot();

    %% Callbacks & Solvers
    function updatePlot(~, ~)
        flds = fieldnames(sliders);
        for idx = 1:numel(flds)
            params.(flds{idx}) = get(sliders.(flds{idx}), 'Value'); 
        end
        params.scenario = get(ui_scenario, 'Value');
        
        [t, y] = solveSystem(params);
       
        LN = y(:,1); LP = y(:,2); B = y(:,3); CA = y(:,4); CM = y(:,5);
        TotalTumor = LN + LP;
        
        mrd_threshold = (params.mrdPct / 100) * params.Lmax;
        
        cla(axTumor); hold(axTumor, 'on');
        
        yline(axTumor, mrd_threshold, 'r--', 'MRD Detection Limit', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
        
        h_tot = plot(axTumor, t, TotalTumor, '-', 'Color', [0 0 0], 'LineWidth', 1.5); % Black Solid
        h_ln  = plot(axTumor, t, LN, '--', 'Color', [77 175 74]/255, 'LineWidth', 2);  % Green Dashed
        h_lp  = plot(axTumor, t, LP, '-', 'Color', [77 175 74]/255, 'LineWidth', 2);   % Green Solid
        h_b   = plot(axTumor, t, B, ':', 'Color', [55 126 184]/255, 'LineWidth', 2);   % Blue Dotted
        
        xlabel(axTumor, 'Time (Days)'); 
        ylabel(axTumor, 'Cell Count'); 
        title(axTumor, 'Tumor and B-Cell Dynamics');
        ylim(axTumor, [0, params.Yend]);
        
        legend(axTumor, [h_tot, h_ln, h_lp, h_b], ...
            {'Total Tumor', '$L_N$ (Antigen-Negative)', '$L_P$ (Antigen-Positive)', '$B$ (Healthy)'}, ...
            'Location', 'best','Interpreter','latex');
        grid(axTumor, 'on');
        
        cla(axCAR); hold(axCAR, 'on');
        h_ca = plot(axCAR, t, CA, '-', 'Color', [255 127 0]/255, 'LineWidth', 2);     % Orange Solid
        h_cm = plot(axCAR, t, CM, '--', 'Color', [255 165 0]/255, 'LineWidth', 2);    % Orange Dashed
        
        xlabel(axCAR, 'Time (Days)'); 
        ylabel(axCAR, 'CAR-T Cell Count'); 
        title(axCAR, 'CAR-T Cell Dynamics');
        legend(axCAR, [h_ca, h_cm], {'$C_A$ (Active)', '$C_M$ (Memory)'}, 'Location', 'best', 'Interpreter','latex');
        ylim(axCAR, [0, params.YendCAR]);
        grid(axCAR, 'on');
    end

function exportAll(~, ~)
        % Prompt user for a base filename
        [filename, pathname] = uiputfile('*.csv', 'Save Base Name For Params & Images', 'CAR_T_Sim.csv');
        
        if isequal(filename, 0) || isequal(pathname, 0)
            return; 
        end
        
        % Extract just the name (without extension) to base all 3 files on
        [~, name, ~] = fileparts(filename);
        
        %% 1. EXPORT CSV DATA (PARAMETERS ONLY)
        flds = fieldnames(params);
        vals = zeros(length(flds), 1);
        for i = 1:length(flds)
            vals(i) = params.(flds{i});
        end
        paramTable = table(flds, vals, 'VariableNames', {'Parameter', 'Value'});
        writetable(paramTable, fullfile(pathname, [name, '.csv']));
        
        %% 2. EXPORT CLEAN IMAGES (Using copyobj instead of re-plotting)
        img_tumor_path = fullfile(pathname, [name, '_Tumor_Plot.png']);
        img_car_path   = fullfile(pathname, [name, '_CART_Plot.png']);
        
        try
            % Create a temporary, hidden figure with perfect proportions
            tempFig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 800, 600]);
            
            % --- Export Plot 1: Tumor Dynamics ---
            % Copy the exact existing axis to the hidden figure
            axTemp1 = copyobj(axTumor, tempFig);
            % Reset its position so it fills the standard 800x600 image perfectly
            set(axTemp1, 'Units', 'normalized', 'Position', [0.15 0.15 0.75 0.75]);
            exportgraphics(axTemp1, img_tumor_path, 'Resolution', 300);
            
            % Delete the copied axis to make room for the next one
            delete(axTemp1); 
            
            % --- Export Plot 2: CAR-T Dynamics ---
            % Copy the exact existing axis to the hidden figure
            axTemp2 = copyobj(axCAR, tempFig);
            set(axTemp2, 'Units', 'normalized', 'Position', [0.15 0.15 0.75 0.75]);
            exportgraphics(axTemp2, img_car_path, 'Resolution', 300);
            
            % Close the hidden figure
            close(tempFig);
            
            disp(['Success! Saved Parameters and 2 perfectly cropped Images to: ', pathname]);
            
        catch ME
            disp(['Error saving images: ', ME.message]);
            % Ensure the hidden figure closes even if there is an error
            if exist('tempFig', 'var') && isvalid(tempFig)
                close(tempFig);
            end
        end
    end

    function resetSliders(~, ~)
        flds = fieldnames(sliders);
        for idx = 1:numel(flds)
            set(sliders.(flds{idx}), 'Value', defaultParams.(flds{idx}));
        end
        set(ui_scenario, 'Value', defaultParams.scenario);
        
        updatePlot();
        
        for idx = 1:numel(flds)
            notify(sliders.(flds{idx}), 'Action'); 
        end
    end

    function [t, y] = solveSystem(p)
        switch p.scenario
            case 1 % 100% Antigen-Positive
                LN0 = 0;
                LP0 = p.L0;
            case 2 % Mixed Clones
                LN0 = 1;
                LP0 = p.L0 ;
            case 3 % 100% Antigen-Negative
                LN0 = p.L0;
                LP0 = 0;
        end
        
        opts = odeset('NonNegative', 1:5, 'RelTol', 1e-4, 'AbsTol', 1e-6);
        y0 = [LN0, LP0, p.B0, p.CAR0, p.CM0];
        
        [t, y] = ode15s(@(t, y) odes(t, y, p), [0, p.Tend], y0, opts);
        
        y(y(:,1) < 1, 1) = 0; % LN
        y(y(:,2) < 1, 2) = 0; % LP
    end

    function dydt = odes(~, y, p)
        LN = y(1); LP = y(2); B = y(3); CA = y(4); CM = y(5);
        
        fExpr = (B + LP) / (p.k + B + LP);
        efExpr = CA * (B + LP) / (p.h + B + LP);
        
        dLN = p.rhoL * LN * (1 - (LN + LP)/p.Lmax) + efExpr * p.delta * LP - (1 - efExpr) * p.eps * LN;
        dLP = p.rhoL * LP * (1 - (LP + LN)/p.Lmax) - efExpr * (p.alpha + p.delta) * LP + (1 - efExpr) * p.eps * LN;
        dB = p.BStem + p.rhoB * B * (1 - (B + LP + LN)/(p.Bmax + LP + LN)) - B/p.tauB - p.alpha * efExpr * B;
        dCA = fExpr * (p.rhoCA * CA + p.gammaMA * CM) - CA/p.tauCA - (1 - fExpr) * p.gammaAM * CA;
        dCM = (1 - fExpr) * p.gammaAM * CA - CM/p.tauCM - fExpr * p.gammaMA * CM;
        
        dydt = [dLN; dLP; dB; dCA; dCM];
    end

    %% UI Helper Functions
    function sl = createSlider(p, l, minv, maxv, initv, cb, yp)
        uicontrol('Parent', p, 'Style', 'text', 'String', l, 'Units', 'normalized', ...
            'Position', [0.02 yp 0.25 0.04], 'FontSize', 9, 'HorizontalAlignment', 'right');
        
        sl = uicontrol('Parent', p, 'Style', 'slider', 'Min', minv, 'Max', maxv, 'Value', initv, ...
            'Units', 'normalized', 'Position', [0.30 yp 0.40 0.04], 'Callback', cb);
            
        et = uicontrol('Parent', p, 'Style', 'edit', 'String', num2str(initv, '%g'), ...
            'Units', 'normalized', 'Position', [0.72 yp 0.25 0.04], 'FontSize', 9, ...
            'Callback', @(src, event) editBoxCallback(src, sl, minv, maxv, cb));
            
        addlistener(sl, 'Value', 'PostSet', @(src, event) set(et, 'String', num2str(get(sl, 'Value'), '%g')));
    end

    function editBoxCallback(src, sl, minv, maxv, cb)
        val = str2double(get(src, 'String'));
        if isnan(val) || val < minv || val > maxv
            set(src, 'String', num2str(get(sl, 'Value'), '%g')); 
        else
            set(sl, 'Value', val);
            cb(); 
        end
    end
end