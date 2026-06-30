%% Sobol Sensitivity Analysis - Stabilized ODE Bivariate Version
clear; clc; close all;

%% 1. Parameter Definitions
paramNames = {'$B_\text{Stem}$', '$B_\max$', '$k$', '$\rho_{CA}$', '$\gamma_{AM}$', '$\tau_{CA}$', '$\alpha$', '$h$', '$B_0$', '$L_0$', '$\epsilon$','$\delta$'};
nparam = numel(paramNames);
L = 1000; % Base number of samples

% Define parameter ranges (Log-sampled parameters use the exponent: e.g. 8 = 10^8)
paramRanges = [
    1e6,   1e9;    % 1: BStem
    1e9,   1e11;   % 2: Bmax
    8,     12;     % 3: k (Log-sampled: 10^8 to 10^12)
    0.01,  1.0;    % 4: rhoCA
    1e-4,  1e-2;   % 5: gammaAM
    1,     15;     % 6: tauCA
    3e-11, 3e-9;    % 7: alpha 
    9,     12;     % 8: h (Log-sampled: 10^13 to 10^15)
    1e6,   1e9;    % 9: B0
    3,     9;      % 10: L0 (Log-sampled: 10^6 to 10^9)
    0,     1e-2    % 11: eps
    0,     1e-15    % 11: delta

];

%% 2. Saltelli Sampling Strategy (Matrices A, B, and Ab)
A = zeros(L, nparam); B = zeros(L, nparam);
for i = 1:nparam
    A(:, i) = paramRanges(i, 1) + rand(L, 1) .* (paramRanges(i, 2) - paramRanges(i, 1));
    B(:, i) = paramRanges(i, 1) + rand(L, 1) .* (paramRanges(i, 2) - paramRanges(i, 1));
end

Ab = cell(1, nparam);
for i = 1:nparam
    Ab{i} = A; 
    Ab{i}(:, i) = B(:, i); 
end

%% 3. Evaluate Model (Parallel Computing)
Tdays = linspace(0, 365*5, 10); % Matching the PDE timeline
outputNames = {'Antigen-Negative (LN)', 'Antigen-Positive (LP)', 'Total Tumor', ...
               'Healthy B-Cells', 'Active CAR-T', 'Memory CAR-T'};
numOutputs = length(outputNames);
p = length(Tdays);

YA = zeros(L, numOutputs, p);
YB = zeros(L, numOutputs, p);
YAb = zeros(L, nparam, numOutputs, p);

fprintf('Starting %d parallel ODE evaluations...\n', L * (2 + nparam));
tic;
parfor j = 1:L
    YA(j, :, :) = runODE(A(j, :), Tdays);
    YB(j, :, :) = runODE(B(j, :), Tdays);
    % Temporary storage for Ab iterations to keep parfor happy
    tempYAb = zeros(nparam, numOutputs, p);
    for i = 1:nparam
        tempYAb(i, :, :) = runODE(Ab{i}(j, :), Tdays);
    end
    YAb(j, :, :, :) = tempYAb;
end
toc;

%% 4. Variance-Based Index Calculation 
Si = zeros(p, nparam, numOutputs); 
STi = zeros(p, nparam, numOutputs);

for c = 1:numOutputs
    for t_idx = 1:p
        ya = YA(:, c, t_idx);
        yb = YB(:, c, t_idx);
        
        % Clean base vectors
        ya(isnan(ya) | isinf(ya)) = 0;
        yb(isnan(yb) | isinf(yb)) = 0;
        
        all_Y = [ya; yb];
        VY = var(all_Y); 
        
        if VY < 1e-12 || isnan(VY), continue; end 
        
        for j = 1:nparam
            yab = squeeze(YAb(:, j, c, t_idx));
            yab(isnan(yab) | isinf(yab)) = 0;
            
            
            % STi: Uses (A - Ab) because they ONLY differ in parameter j
            V_Total_j = mean((ya - yab).^2) / 2;
            STi(t_idx, j, c) = V_Total_j / VY;
            
            % Si: Uses 1 - (B - Ab) because they are the SAME only in parameter j
            V_minus_j = mean((yb - yab).^2) / 2;
            Si(t_idx, j, c) = 1 - (V_minus_j / VY);
        end
    end
end

%% 5. Rainbow Visualization 
rainbowMap = turbo(nparam); 
for c = 1:numOutputs
    figure('Color','w','Name',['Sensitivity Analysis: ' outputNames{c}],'Position', [100 100 1100 500]);
    
    cleanSi  = min(1, max(0, Si(:,:,c))); cleanSi(isnan(cleanSi)) = 0;
    cleanSTi = min(1, max(0, STi(:,:,c))); cleanSTi(isnan(cleanSTi)) = 0;


     % --- Dynamic Y-Limit Calculation for Si ---
    stackSi = sum(cleanSi, 2);
    limitSi = max(stackSi);
    if limitSi > 0, limitSi = limitSi * 1.1; else, limitSi = 1.2; end
    
    % --- Dynamic Y-Limit Calculation for STi ---
    stackSTi = sum(cleanSTi, 2);
    limitSTi = max(stackSTi);
    if limitSTi > 0, limitSTi = limitSTi * 1.1; else, limitSTi = 1.2; end

    subplot(1,2,1);
    h1 = bar(Tdays/365, cleanSi, 'stacked', 'EdgeColor', 'none');
    for k=1:nparam, h1(k).FaceColor = rainbowMap(k,:); end
    ylim([0 limitSi]); grid on;
    title(['Direct Impact (S_i): ', outputNames{c}]);
    ylabel('Sensitivity Index'); xlabel('Years');
    
    subplot(1,2,2);
    h2 = bar(Tdays/365, cleanSTi, 'stacked', 'EdgeColor', 'none');
    for k=1:nparam, h2(k).FaceColor = rainbowMap(k,:); end
    ylim([0 limitSTi]); grid on;
    title(['Total Impact (S_{Ti}): ', outputNames{c}]);
    ylabel('Total Sensitivity Index'); xlabel('Years');
    
    legend(paramNames, 'Location', 'eastoutside', 'Interpreter','latex');
end

%% 6. Final Interpretation Table & TXT Export
% 
folderName = 'Sensitivity_Results_ODE';
if ~exist(folderName, 'dir'), mkdir(folderName); end
txtFile = fullfile(folderName, 'Reporte_Parametros.txt');
fid = fopen(txtFile, 'w');

% 
header = sprintf('\n========================================================\n');
header = [header, sprintf('   TOP 3 PARAMETERS BY CELL TYPE AND YEAR (S_Ti Index)   \n')];
header = [header, sprintf('========================================================\n')];

fprintf(header); 
fprintf(fid, '%s', header);

checkYears = [1, 3, 5]; 
% 
if exist('cellTypes','var'), names = cellTypes; else, names = outputNames; end

for c = 1:length(names)
    targetHeader = sprintf('\n--- TARGET: %s ---\n', upper(names{c}));
    tableCols = sprintf('%-10s | %-15s | %-15s | %-15s\n', 'Timeline', 'Rank 1', 'Rank 2', 'Rank 3');
    divider = sprintf('--------------------------------------------------------\n');
    
    % 
    fprintf('%s%s%s', targetHeader, tableCols, divider);
    fprintf(fid, '%s%s%s', targetHeader, tableCols, divider);
    
    for y = checkYears
        [~, tIdx] = min(abs((Tdays/365) - y));
        currentST = min(1, max(0, STi(tIdx, :, c)));
        [sortedVals, sortedIdx] = sort(currentST, 'descend');
        
        row = sprintf('Year %d     | %s (%.2f)  | %s (%.2f)  | %s (%.2f)\n', ...
            y, ...
            paramNames{sortedIdx(1)}, sortedVals(1), ...
            paramNames{sortedIdx(2)}, sortedVals(2), ...
            paramNames{sortedIdx(3)}, sortedVals(3));
            
        fprintf(row);
        fprintf(fid, '%s', row);
    end
end

%% 7. Heatmap Visualization (WITH NUMBERS)
for c = 1:numOutputs
    figure('Color','w','Name',['Heatmap: ' outputNames{c}]);
    
    heatmapData = min(1, max(0, STi(:,:,c)')); 
    
    imagesc(Tdays/365, 1:nparam, heatmapData);
    colormap(turbo); colorbar; clim([0 1]); 
    
    set(gca, 'YTick', 1:nparam, 'YTickLabel', paramNames , 'TickLabelInterpreter','latex');
    xlabel('Time (Years)'); ylabel('Parameters');
    title(['S_{Ti} Intensity Heatmap: ', outputNames{c}]);
    
    % Add the numeric labels back to the heatmap
    [rows, cols] = size(heatmapData);
    for i = 1:rows
        for j = 1:cols
            if heatmapData(i,j) > 0.05 
                text(Tdays(j)/365, i, num2str(heatmapData(i,j), '%0.2f'), ...
                    'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', 'w');
            end
        end
    end
end


%% 8. Export All Figures
% Create a folder to save figures if it doesn't exist
folderName = 'Sensitivity_Results_ODE';
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

% Get handles to all open figures
figHandles = findobj('Type', 'figure');

% Loop through figures and save
for i = 1:length(figHandles)
    % Select the figure
    fig = figHandles(i);
    
    % Clean up the name for the filename (remove LaTeX or special characters)
    rawName = get(fig, 'Name');
    safeName = regexprep(rawName, '[^a-zA-Z0-9]', '_'); 
    
    % Use the figure number if no name is assigned
    if isempty(safeName)
        safeName = sprintf('Figure_%d', fig.Number);
    end
    
    % Define the full path
    fileName = fullfile(folderName, [safeName, '.png']);
    
    % Save the figure (Resolution set to 300 DPI for clarity)
    exportgraphics(fig, fileName, 'Resolution', 300);
    
    fprintf('Exported: %s\n', fileName);
end

fprintf('Done! All %d figures are saved in the "%s" folder.\n', length(figHandles), folderName);

%% ---------- CORE SOLVER FUNCTIONS ----------
function OutputMatrix = runODE(pv, Tdays)
    p.delta = 1e-20; p.CAR0 = 1e7; p.CM0 = 0; p.rhoL = 1/30; 
    p.rhoB = log(2)/8; p.tauCM = 300; p.gammaMA = 0.33; p.Lmax = 1e12;
    p.tauB = 60;
    
    p.BStem   = pv(1);
    p.Bmax    = pv(2);
    p.k       = 10^pv(3);   
    p.rhoCA   = pv(4);
    p.gammaAM = pv(5);
    p.tauCA   = pv(6);
    p.alpha   = pv(7);
    p.h       = 10^pv(8);   
    p.B0      = pv(9);
    L0_val    = 10^pv(10);  
    p.eps     = pv(11);
    p.delta   = pv(12);
    
    LN0 = L0_val * 0.5;
    LP0 = L0_val * 0.5;
    y0 = [LN0, LP0, p.B0, p.CAR0, p.CM0];
    
    opts = odeset('RelTol', 1e-3, 'AbsTol', 1e-5);
    try
        [t, y] = ode15s(@(t, y) odes(t, y, p), [0, max(Tdays)], y0, opts);
        
        LN = interp1(t, y(:,1), Tdays, 'linear', 'extrap');
        LP = interp1(t, y(:,2), Tdays, 'linear', 'extrap');
        B  = interp1(t, y(:,3), Tdays, 'linear', 'extrap');
        CA = interp1(t, y(:,4), Tdays, 'linear', 'extrap');
        CM = interp1(t, y(:,5), Tdays, 'linear', 'extrap');
        
        TotalTumor = LN + LP;
        raw_Y = [LN; LP; TotalTumor; B; CA; CM];
        raw_Y(raw_Y > 1e15) = 1e15; 
        
        % --- THE LOG10 FIX ---
        OutputMatrix = log10(max(0, raw_Y) + 1);
    catch
        OutputMatrix = NaN(6, length(Tdays)); 
    end
end

function dydt = odes(~, y, p)
    LN = max(0, y(1)); LP = max(0, y(2)); B  = max(0, y(3)); 
    CA = max(0, y(4)); CM = max(0, y(5));
    fExpr = (B + LP) / (p.k + B + LP);
    efExpr = CA * (B + LP) / (p.h + B + LP);
    
    dLN = p.rhoL * LN * (1 - (LN + LP)/p.Lmax) + efExpr * p.delta * LP - (1 - efExpr) * p.eps * LN;
    dLP = p.rhoL * LP * (1 - (LP + LN)/p.Lmax) - efExpr * (p.alpha + p.delta) * LP + (1 - efExpr) * p.eps * LN;
    dB  = p.BStem + p.rhoB * B * (1 - (B + LP + LN)/(p.Bmax + LP + LN)) - B/p.tauB - p.alpha * efExpr * B;
    dCA = fExpr * (p.rhoCA * CA + p.gammaMA * CM) - CA/p.tauCA - (1 - fExpr) * p.gammaAM * CA;
    dCM = (1 - fExpr) * p.gammaAM * CA - CM/p.tauCM - fExpr * p.gammaMA * CM;
    dydt = [dLN; dLP; dB; dCA; dCM];
end
