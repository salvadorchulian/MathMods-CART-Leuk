function CAR_T_ODE_Antigen_Proportion
    clc; close all;

    %% 1. Interactive Choice
    scenarios = {'100percent Antigen-Positive', '50_50 Mixed Clones', '100percent Antigen-Negative'};
    choice = questdlg('Select Initial Leukaemia Condition:', ...
        'Initial Condition', ...
        scenarios{1}, scenarios{2}, scenarios{3}, ...
        scenarios{2});
    
    if isempty(choice), return; end 
    scenario_idx = find(strcmp(scenarios, choice));

    %% 2. Configuration (10-Year Timeline)
    years = [0, 1/3, 1/2, 1, 2, 3, 5, 7, 10]; 
    timePoints = years * 365;               
    numT = length(timePoints);
    res = 20; % Resolution of the heatmap (20x20 grid)
    
    % Sweeping h vs k (using logspace for both due to large magnitudes)
    h_range = logspace(13, 15, res);      
    k_range = logspace(8, 12, res);    
    [H, K] = meshgrid(h_range, k_range);
    
    Total_L_Store = zeros(res, res, numT);
    Composition_Store = zeros(res, res, numT);
    
    % Core Model Parameters from Interactive Bivariate ODE
    p = struct('Tend', max(timePoints), 'eps', 0.00348, 'delta', 1e-20, ...
        'L0', 2e7, 'BStem', 1e8, 'B0', 1e8, 'CAR0', 1e7, 'CM0', 0, ...
        'alpha', 0.5, 'k', 1e11, 'h', 5.23e14, 'rhoCA', 0.95, ...
        'Bmax', 1e10, 'rhoL', 1/30, 'rhoB', log(2)/8, ...
        'tauB', 60, 'tauCA', 6.5, 'tauCM', 300, ...
        'gammaAM', 0.001, 'gammaMA', 0.33, 'Lmax', 1e12, ...
        'scenario', scenario_idx);

    %% 3. Parameter Sweep
    fprintf('Simulating 10-year horizon for: %s\n', choice);
    tic
    for i = 1:res
        for j = 1:res
            cp = p;
            cp.h = H(i,j); % Now sweeping h
            cp.k = K(i,j);
            
            [t, y] = solveSystem(cp);
            
            for t_idx = 1:numT
                % Interpolate to find state at specific time point
                if timePoints(t_idx) == 0
                    y_at_t = y(1,:);
                else
                    y_at_t = interp1(t, y, timePoints(t_idx), 'linear', 'extrap');
                end
                
                % Prevent negative populations (artifact of ODE solver)
                y_at_t(y_at_t < 1) = 0; 
                
                LN = y_at_t(1);
                LP = y_at_t(2);
                sumL = LN + LP;
                
                Total_L_Store(i,j,t_idx) = sumL;
                
                % Grey Area Logic: If tumor is essentially cleared (<= 100 cells)
                if sumL > 100 
                    Composition_Store(i,j,t_idx) = LP / sumL;
                else
                    Composition_Store(i,j,t_idx) = NaN; % Becomes Grey
                end
            end
        end
        if mod(i,4)==0, fprintf('Progress: %d%%\n', round(i/res*100)); end
    end
    toc
%% 4. Visualization 
    % Adjust figure size: narrow and tall [Left Bottom Width Height]
    fig = figure('Color', 'w', 'Position', [100 50 800 1800]); 
    %sgtitle(['Long-term Evolution: ', choice], 'FontWeight', 'bold', 'FontSize', 14);
    
    % Layout Math
    top_margin = 0.06;
    bottom_margin = 0.05;
    left_margin = 0.2;
    right_margin = 0.12;
    spacing_v = 0.01; % Vertical spacing between time points
    spacing_h = 0.01; % Horizontal spacing between the two columns
    
    available_height = 1 - top_margin - bottom_margin;
    % Height per row (each row is one time point)
    height = (available_height - (numT-1)*spacing_v) / numT;
    % Width per column (2 columns total)
    width = (1 - left_margin - right_margin - spacing_h) / 2;
    
    max_burden = max(log10(Total_L_Store(:) + 1));
    if max_burden < 1, max_burden = 12; end 
    
    for t_idx = 1:numT
        % Calculate vertical position (starting from top)
        y_pos = 1 - top_margin - t_idx*height - (t_idx-1)*spacing_v;
        
        % --- Smart Label Logic ---
        if years(t_idx) == 0, timeLabel = 'Day 0';
        elseif years(t_idx) < 1, timeLabel = sprintf('%d Mo', round(years(t_idx) * 12));
        else, timeLabel = sprintf('Year %d', years(t_idx)); 
        end
        
        % --- Column 1: Tumor Burden ---
        ax1 = axes('Position', [left_margin, y_pos, width, height]);
        imagesc(h_range, log10(k_range), log10(Total_L_Store(:,:,t_idx) + 1));
        set(ax1, 'YDir', 'normal'); clim(ax1, [0, max_burden]);
        colormap(ax1, parula);
        
        % Add Time Label to the left of the first column
        ylabel(['\bf', timeLabel, '\rm', newline, 'log_{10}(k)'], 'FontSize', 9);
        if t_idx == numT, xlabel('h'); else, set(ax1, 'XTickLabel', []); end
        if t_idx == 1, title('Tumor Burden', 'FontSize', 11); end

        % --- Column 2: Mean Antigen ---
        ax2 = axes('Position', [left_margin + width + spacing_h, y_pos, width, height]);
        h_img = imagesc(h_range, log10(k_range), Composition_Store(:,:,t_idx));
        set(h_img, 'AlphaData', ~isnan(Composition_Store(:,:,t_idx)));
        set(ax2, 'YDir', 'normal', 'Color', [0.8 0.8 0.8], 'YTickLabel', []); 
        colormap(ax2, custom_rb_map());
        clim(ax2, [0 1]);
        
        if t_idx == numT, xlabel('h'); else, set(ax2, 'XTickLabel', []); end
        if t_idx == 1, title('Mean Antigen', 'FontSize', 11); end
        
        % Colorbars at the bottom or top? Let's put them next to the first/last plots
        if t_idx == floor(numT/2)
             cb1 = colorbar(ax1, 'westoutside', 'Position', [0.1, y_pos, 0.02, height*2]);
             ylabel(cb1, 'log_{10}(Cells)');
             cb2 = colorbar(ax2, 'eastoutside', 'Position', [1-0.1, y_pos, 0.02, height*2]);
             ylabel(cb2, 'Antigen Expr.');
        end
    end
    %% Exporting figures
    root="/Users/salvador/Library/CloudStorage/Dropbox/MATLAB/CART";
    exportgraphics(fig, root+"/ODE_Ant_"+choice+".png", 'Resolution', 300);
end

%% --- Helper Functions ---
function [t, y] = solveSystem(p)
    switch p.scenario
        case 1 % 100% Antigen-Positive
            LN0 = 0;
            LP0 = p.L0;
        case 2 % Mixed Clones
            LN0 = p.L0 * 0.5;
            LP0 = p.L0 * 0.5;
        case 3 % 100% Antigen-Negative
            LN0 = p.L0;
            LP0 = 0;
    end
    
    y0 = [LN0, LP0, p.B0, p.CAR0, p.CM0];
    tspan = [0 p.Tend];
    
    % Using 'AbsTol', 1e-15 to prevent sub-cell fraction integration crashes
    opts = odeset('NonNegative', 1:5, 'RelTol', 1e-3, 'AbsTol', 1e-13);
    [t, y] = ode15s(@(t, y) odes(t, y, p), tspan, y0, opts);
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

function cmap = custom_rb_map()
    % Restored Original Blue-White-Red Colormap
    r = [linspace(0,1,32)'; linspace(1,1,32)'];
    g = [linspace(0,1,32)'; linspace(1,0,32)'];
    b = [linspace(1,1,32)'; linspace(1,0,32)'];
    cmap = [r, g, b];
end