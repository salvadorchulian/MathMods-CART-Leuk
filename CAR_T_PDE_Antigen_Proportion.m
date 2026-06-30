function CAR_T_PDE_Antigen_Proportion
    clc; close all;

    %% 1. Interactive Choice
    choice = questdlg('Select Initial Leukaemia Condition:', ...
        'Initial Condition', ...
        'Single Clone (0.6)', 'Single Clone (0.4)', 'Two Clones (0.4 & 0.6)', ...
        'Two Clones (0.4 & 0.6)');
    
    if isempty(choice), return; end 

    %% 2. Configuration (10-Year Timeline)
    % Using fractions to ensure 4 and 6 month marks are precise
    years = [0, 1/3, 1/2, 1, 2, 3, 5, 7, 10]; 
    timePoints = years * 365;               
    numT = length(timePoints);
    res = 20; 
    
    h_range = linspace(0, 1, res);      
    k_range = logspace(8, 12, res);    
    [H, K] = meshgrid(h_range, k_range);
    
    Total_L_Store = zeros(res, res, numT);
    Mean_Antigen_Store = zeros(res, res, numT);
    
    % Core Model Parameters
    p = struct('tend', max(timePoints), 'Nx', 100, 'sigma', 1/20, 'alpha', 3e-10, ...
        'rhoCA', 0.9, 'rhoB', log(2)/8, 'rhoL', 1/30, ...
        'tauCA', 6.5, 'tauCM', 300, 'tauB', 60, ...
        'gammaAM', 0.001, 'gammaMA', 0.33, 'h', 0.75, ...
        'k', 1e10, 'LMax', 1e12, 'CA0', 1e7, 'CM0', 0, ...
        'L0', 1e7, 'B0', 1e7, 'BMax', 5e9, 'm', 1, ...
        'x0', 50, 's', 1e-3, 'BStem', 1e8, 'scenario', choice);

    %% 3. Parameter Sweep
    fprintf('Simulating 10-year horizon for: %s\n', choice);
    tic
    for i = 1:res
        for j = 1:res
            cp = p;
            cp.h = H(i,j); cp.k = K(i,j);
            [t, ~, ~, L_dist, ~] = solveSystem(cp);
            x_vec = linspace(0, 1, p.Nx);
            
            for t_idx = 1:numT
                % Interpolate to find state at specific time point
                if timePoints(t_idx) == 0
                    L_at_t = L_dist(1,:);
                else
                    L_at_t = interp1(t, L_dist, timePoints(t_idx), 'linear', 'extrap');
                end
                
                L_at_t(L_at_t < 1) = 0; 
                sumL = sum(L_at_t);
                Total_L_Store(i,j,t_idx) = sumL;
                
                % Grey Area Logic: sumL <= 100 cells
                if sumL > 100 
                    Mean_Antigen_Store(i,j,t_idx) = sum(x_vec .* L_at_t) / sumL;
                else
                    Mean_Antigen_Store(i,j,t_idx) = NaN; % Becomes Grey
                end
            end
        end
        if mod(i,5)==0, fprintf('Progress: %d%%\n', round(i/res*100)); end
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
        h_img = imagesc(h_range, log10(k_range), Mean_Antigen_Store(:,:,t_idx));
        set(h_img, 'AlphaData', ~isnan(Mean_Antigen_Store(:,:,t_idx)));
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
    exportgraphics(fig, root+"/PDE_Ant_"+choice+".png", 'Resolution', 300);
end

%% --- Helper Functions ---
function [t, CA, CM, L, B] = solveSystem(p)
    tend_days = p.tend; % Already in days
    tspan = [0 tend_days];
    xspan = linspace(0, 1, p.Nx);
    dx = xspan(2) - xspan(1); 
    
    % Initial Tumor (L) using normalized Gaussian
    switch p.scenario
        case 'Single Clone (0.6)', ICL = gaussian(xspan, 0.6, p.sigma)*p.L0;
        case 'Single Clone (0.4)', ICL = gaussian(xspan, 0.4, p.sigma)*p.L0;
        case 'Two Clones (0.4 & 0.6)', ICL = (gaussian(xspan,0.4,p.sigma)*0.5 + gaussian(xspan,0.6,p.sigma)*0.5)*p.L0;
    end
  
    % Initial B-Cells and the B-Cell Niche (hB)
    ICB = gaussian(xspan, 0.2, p.sigma)*p.B0;
    hB_vec = (gaussian(xspan,0.25,p.sigma)*0.1 + gaussian(xspan,0.5,p.sigma)*0.65 + gaussian(xspan,0.75,p.sigma)*0.25)*p.BMax;
    
    IC = [p.CA0; p.CM0; ICL(:); ICB(:)];
    opts = odeset('RelTol', 1e-5, 'AbsTol', 1e-8);
    [t, res] = ode45(@(t,e) odefun_core(t, e, dx, p, hB_vec(:)), tspan, IC, opts);
    
    CA = res(:,1); 
    CM = res(:,2); 
    L  = res(:, 3 : p.Nx+2); 
    B  = res(:, p.Nx+3 : end);
end

function ne = odefun_core(~, e, dx, p, hB_vec)
    CA = e(1); CM = e(2); 
    L = e(3:p.Nx+2); B = e(p.Nx+3:end); 
    
    % Thresholding to prevent negative populations
    L(L<1) = 0; B(B<1) = 0; 
    x = (dx*(1:p.Nx))';
    
    % Global interactions
    IL = sum(L)*dx; 
    IB = sum(B)*dx; 
    IN = IL + IB; 
    F = IN / (p.k + IN); % CAR-T Activation Factor
    
    % Spatial Stem Cell Source
    Sx = p.s ./ (1 + exp(p.m * (x - p.x0)));
    
    % Dynamics
    dL = p.rhoL.*L.*(1 - IL/p.LMax) - p.alpha.*(x./(x + p.h + 1e-5)).*L.*CA;
    dB = p.BStem*Sx + p.rhoB.*B.*(1 - (B+IL)./(hB_vec+IL)) - p.alpha.*(x./(x+p.h+1e-5)).*B.*CA - B/p.tauB;
    
    dCA = F*(p.rhoCA*CA + p.gammaMA*CM) - (1/p.tauCA)*CA - (1-F)*p.gammaAM*CA;
    dCM = (1-F)*p.gammaAM*CA - (1/p.tauCM)*CM - F*p.gammaMA*CM;
    
    ne = [dCA; dCM; dL; dB]; 
end

function g = gaussian(x, m, s)
    g = exp(-((x-m)/s).^2/2) ./ (s*sqrt(2*pi)); 
end

function cmap = custom_rb_map()
    r = [linspace(0,1,32)'; linspace(1,1,32)'];
    g = [linspace(0,1,32)'; linspace(1,0,32)'];
    b = [linspace(1,1,32)'; linspace(1,0,32)'];
    cmap = [r, g, b];
end