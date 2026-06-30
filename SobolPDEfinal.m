%% Sobol Sensitivity Analysis - Stabilized Global PDE Version (Log10 Output Fix)
clear; clc; close all;

%% 1. Parameter Definitions (Aligned with Program 1)
paramNames = {'$B_\text{Stem}$', '$B_\max$', '$s$', '$x_0$', '$m$', '$k$', '$\rho_{CA}$', '$\gamma_{AM}$', '$\tau_{CA}$', '$\alpha$', '$h$', '$B_0$', '$L_0$'};
nparam = numel(paramNames);
l = 1000; % Number of samples 

% Define ranges 
paramRanges = [
    1e6,  1e9;    % BStem
    1e9,  1e11;   % BMax
    1e-4, 1e-1;   % s
    0,    100;    % x0
    0,    10;     % m
    1e8,  1e12;   % k (Log-sampled)
    0.01,  1;    % rhoCA
    1e-4, 1e-2;   % gammaAM
    1,    15;     % tauCA
    3e-11, 3e-9; % alpha 
    0,    1;      % h
    1e6,  1e9;    % B0
    1e3,  1e9;    % L0 (Log-sampled)
];



%% 2. Saltelli Sampling Strategy (Matrices A, B, and Ab)
A = zeros(l, nparam); B = zeros(l, nparam);
for i = 1:nparam
    lb = paramRanges(i,1); ub = paramRanges(i,2);
    if ismember(i, [6, 10, 13]) 
        A(:,i) = 10.^(log10(lb) + (log10(ub)-log10(lb)).*rand(l,1));
        B(:,i) = 10.^(log10(lb) + (log10(ub)-log10(lb)).*rand(l,1));
    else
        A(:,i) = lb + (ub-lb).*rand(l,1);
        B(:,i) = lb + (ub-lb).*rand(l,1);
    end
end

Ab = cell(1, nparam);
for j = 1:nparam
    Ab{j} = A; Ab{j}(:,j) = B(:,j);
end

%% 3. Parallel Simulation Execution
Tdays = linspace(0, 365*5, 10); % Simulation time points
p = length(Tdays);
Y_A = zeros(l, p, 6); Y_B = zeros(l, p, 6); Y_Ab = zeros(l, p, nparam, 6);

fprintf('Starting %d parallel PDE simulations...\n', l * (nparam + 2));
tic
parfor i = 1:l
    Y_A(i,:,:) = runFastIteration(A(i,:), Tdays);
    Y_B(i,:,:) = runFastIteration(B(i,:), Tdays);
    for j = 1:nparam
        Y_Ab(i,:,j,:) = runFastIteration(Ab{j}(i,:), Tdays);
    end
end
toc


%% 4. Variance-Based Index Calculation 
Si = zeros(p, nparam, 6); STi = zeros(p, nparam, 6);

for c = 1:6
    for t = 1:p
        ya = Y_A(:,t,c);
        yb = Y_B(:,t,c);
        
        % Step 1: Clean base vectors (Remove any failed simulations)
        ya(isnan(ya) | isinf(ya)) = 0;
        yb(isnan(yb) | isinf(yb)) = 0;
        
        % Step 2: Calculate Total Variance of the system at this timepoint
        all_Y = [ya; yb];
        VY = var(all_Y); 
        
        if VY < 1e-12 || isnan(VY), continue; end 
        
        for j = 1:nparam
            yab = Y_Ab(:,t,j,c);
            yab(isnan(yab) | isinf(yab)) = 0;
            
            % Total-Order Index (STi): 
            % Measures direct effect + all interactions.
            % Uses (A - Ab) because they ONLY differ in parameter j.
            V_Total_j = mean((ya - yab).^2) / 2;
            STi(t,j,c) = V_Total_j / VY;
            
            % First-Order Index (Si): 
            % Measures the direct effect of parameter j only.
            % Uses 1 - (B - Ab) because they are the SAME only in parameter j.
            V_minus_j = mean((yb - yab).^2) / 2;
            Si(t,j,c) = 1 - (V_minus_j / VY);
        end
    end
end

%% 5. Rainbow Visualization 
cellTypes = {'Activated CAR-T','Memory CAR-T','Total Leukemia','Total B-Cells', ...
             'Leukemia (Antigen < 50%)', 'Leukemia (Antigen >= 50%)'};
rainbowMap = turbo(nparam); 

for c = 1:6
    figure('Color','w','Name',['Sensitivity Analysis: ' cellTypes{c}],'Position', [100 100 1100 500]);
    
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
    title(['Direct Impact (S_i): ', cellTypes{c}]);
    ylabel('Sensitivity Index'); xlabel('Years');
    
    subplot(1,2,2);
    h2 = bar(Tdays/365, cleanSTi, 'stacked', 'EdgeColor', 'none');
    for k=1:nparam, h2(k).FaceColor = rainbowMap(k,:); end
    ylim([0 limitSTi]); grid on;
    title(['Total Impact (S_{Ti}): ', cellTypes{c}]);
    ylabel('Total Sensitivity Index'); xlabel('Years');
    
    legend(paramNames, 'Location', 'eastoutside', 'Interpreter','latex');
end

%% 6. Final Interpretation Table & TXT Export
% 
folderName = 'Sensitivity_Results_PDE';
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

%% 7. Heatmap Visualization 
for c = 1:6
    figure('Color','w','Name',['Heatmap: ' cellTypes{c}]);
    
    heatmapData = min(1, max(0, STi(:,:,c)')); 
    
    imagesc(Tdays/365, 1:nparam, heatmapData);
    colormap(turbo);
    colorbar;
    clim([0 1]); 
    
    set(gca, 'YTick', 1:nparam, 'YTickLabel', paramNames, 'TickLabelInterpreter','latex');
    xlabel('Time (Years)');
    ylabel('Parameters');
    title(['S_{Ti} Intensity Heatmap: ', cellTypes{c}]);
    
    [rows, cols] = size(heatmapData);
    if cols < 25 
        for i = 1:rows
            for j = 1:cols
                if heatmapData(i,j) > 0.05 
                    text(Tdays(j)/365, i, num2str(heatmapData(i,j), '%0.2f'), ...
                        'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', 'w');
                end
            end
        end
    end
end

%% 8. Export All Figures
% Create a folder to save figures if it doesn't exist
folderName = 'Sensitivity_Results_PDE';
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

function Y_out = runFastIteration(pv, Tdays)
    Nx = 50; 
    dx = 1/Nx; xspan = linspace(0, 1, Nx);
    sigma = 0.05; LMax = 1e12; 
    
    ICL = (gaussian(xspan,0.4,sigma)*0.5 + gaussian(xspan,0.6,sigma)*0.5)*pv(13);
    ICB = gaussian(xspan,0.2,sigma)*pv(12);
    hB = (gaussian(xspan,0.25,sigma)*0.1 + gaussian(xspan,0.5,sigma)*0.65 + gaussian(xspan,0.75,sigma)*0.25)*pv(2);
    
    IC = [1e7; 0; ICL(:); ICB(:)]; 
    
    opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-6, 'NonNegative', 1:(2*Nx+2));
    [t, res] = ode45(@(t,e) fast_pde(t,e,dx,Nx,pv,hB(:),LMax), [0 max(Tdays)], IC, opts);


    %
    idx_mid = floor(Nx/2); % Mid point
    
    L_low  = sum(res(:, 3 : 3 + idx_mid - 1), 2) * dx;     % x < 0.5
    L_high = sum(res(:, 3 + idx_mid : Nx + 2), 2) * dx;    % x >= 0.5
    L_tot  = sum(res(:, 3:Nx+2), 2) * dx;  % 
    B_tot  = sum(res(:, Nx+3:end), 2) * dx; % 
    
    
    % Ahora raw_Y tiene 6 columnas
    raw_Y = [interp1(t, res(:,1), Tdays, 'linear', 'extrap')', ... % CAR-T A
             interp1(t, res(:,2), Tdays, 'linear', 'extrap')', ... % CAR-T M
             interp1(t, L_tot, Tdays, 'linear', 'extrap')', ...    % Total L
             interp1(t, B_tot, Tdays, 'linear', 'extrap')', ...    % Total B
             interp1(t, L_low, Tdays, 'linear', 'extrap')', ...    % L < 0.5 
             interp1(t, L_high, Tdays, 'linear', 'extrap')'];      % L >= 0.5 
             
    raw_Y(raw_Y > 1e15) = 1e15; 
    Y_out = log10(raw_Y + 1);
    
end

function ne = fast_pde(~,e,dx,Nx,pv,hB,LMax)
    CA=e(1); CM=e(2); L=e(3:Nx+2); B=e(Nx+3:end);
    x=(dx*(1:Nx))';
    IL=sum(L)*dx; IB=sum(B)*dx; IN=IL+IB;
    F = IN/(pv(6)+IN);
    Sx = pv(3)./(1+exp(pv(5)*(x-pv(4))));
    
    dL = 0.0333*L.*(1-IL/LMax) - pv(10).*(x./(x+pv(11)+1e-5)).*L.*CA;
    dB = pv(1)*Sx + 0.0866*B.*(1-(B+IL)./(hB+IL)) - pv(10).*(x./(x+pv(11)+1e-5)).*B.*CA - B/60;
    
    dCA = F*(pv(7)*CA + 0.33*CM) - (1/pv(9))*CA - (1-F)*pv(8)*CA;
    dCM = (1-F)*pv(8)*CA - (1/300)*CM - F*0.33*CM;
    ne = [dCA; dCM; dL; dB];
end

function g = gaussian(x, m, s), g = exp(-((x-m)/s).^2/2)./(s*sqrt(2.5066)); end