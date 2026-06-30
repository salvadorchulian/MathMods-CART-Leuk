function interactive_PDE_model
clc; close all;
%% Parameters
params = struct( ...
    'tend', 5, 'Nx', 100, 'sigma', 1/20, 'alpha', 3e-10, ...
    'rhoCA', 0.9, 'rhoB', log(2)/8, 'rhoL', 1/30, ...
    'tauCA', 6.5, 'tauCM', 300, 'tauB', 60, ...
    'gammaAM', 0.001, 'gammaMA', 0.33, 'h', 0.75, ...
    'k', 1e10, 'LMax', 1e12, 'CA0', 1e7, 'CM0', 0, ...
    'L0', 1e7, 'B0', 1e7, 'BMax', 5e9, 'm', 1, ...
    'x0', 50, 's', 1e-3, 'BStem', 1e8, ...
    'xmin', 0,  'xmax', 1, ...
    'scenario', 'Two Clones (0.4 & 0.6)', 'mrdPct', 0.1,'exh',1e9);
defaultParams = params; 

%% Figure & Layout
f = figure('Name','Bivariate CAR-T Model (Years)','Position', [50 50 1400 900], 'Color', 'w');
ax1 = axes('Parent',f,'Position',[0.04 0.55 0.55 0.38]); 
ax2 = axes('Parent',f,'Position',[0.04 0.10 0.55 0.38]);  
ax3 = axes('Parent',f,'Position',[0.72 0.08 0.25 0.35]); 
axLeg1 = axes('Parent',f,'Position',[0.62 0.65 0.06 0.18]); 
axLeg2 = axes('Parent',f,'Position',[0.62 0.20 0.06 0.18]);
panel = uipanel('Parent',f,'Title','Controls','Units','normalized','Position',[0.72 0.45 0.26 0.53]);

ypos = 0.9;
uicontrol('Parent',panel,'Style', 'pushbutton', 'String', 'Reset', 'Units', 'normalized', ...
    'Position', [0.55 ypos+0.05 0.2 0.04], 'FontSize', 8, 'Callback', @resetSliders);
uicontrol('Parent',panel,'Style', 'pushbutton', 'String', 'Export Data', 'Units', 'normalized', ...
    'Position', [0.77 ypos+0.05 0.2 0.04], 'FontSize', 8, 'Callback', @exportData);

uicontrol('Parent',panel,'Style','text','String','Condition:','Units','normalized',...
    'Position',[0.05 ypos+0.05 0.15 0.035],'FontSize',7,'HorizontalAlignment','left');
popScenario = uicontrol('Parent',panel,'Style','popupmenu',...
    'String',{'Single Clone (0.6)', 'Single Clone (0.4)', 'Two Clones (0.4 & 0.6)'},...
    'Value', 3, 'Units','normalized','Position',[0.22 ypos+0.05 0.32 0.035],'Callback',@updatePlot);

ypos = ypos - 0.05;
sliders.BStem   = createSlider(panel,'BStem',1e6,1e8,params.BStem,@updatePlot,ypos); 
sliders.BMax    = createSlider(panel,'BMax',1e9,1e11,params.BMax,@updatePlot,ypos); ypos = ypos - 0.04;
sliders.s       = createSlider(panel,'s',1e-4,1,params.s,@updatePlot,ypos); ypos = ypos - 0.04;
sliders.x0      = createSlider(panel,'x0',0,100,params.x0,@updatePlot,ypos); ypos = ypos - 0.04;
sliders.m       = createSlider(panel,'m',0,10,params.m,@updatePlot,ypos); ypos = ypos - 0.04;
sliders.k       = createSlider(panel,'k',1e8,1e13,params.k,@updatePlot,ypos); ypos = ypos - 0.04;
sliders.rhoCA   = createSlider(panel,'rhoCA',0,2,params.rhoCA,@updatePlot,ypos); ypos = ypos - 0.04;
sliders.gammaAM = createSlider(panel,'gammaAM',1e-3,1e-1,params.gammaAM,@updatePlot,ypos); ypos = ypos - 0.04;
sliders.tauCA   = createSlider(panel,'tauCA',1,15,params.tauCA,@updatePlot,ypos); ypos = ypos - 0.04;
sliders.alpha   = createSlider(panel,'alpha',3e-11,3e-9,params.alpha,@updatePlot,ypos); ypos = ypos - 0.04;
sliders.h       = createSlider(panel,'h',0,1,params.h,@updatePlot,ypos); ypos = ypos - 0.04;
sliders.B0      = createSlider(panel,'B0',1e6,1e8,params.B0,@updatePlot,ypos); ypos = ypos - 0.04;
sliders.L0      = createSlider(panel,'L0',0,1e9,params.L0,@updatePlot,ypos); ypos = ypos - 0.04;
sliders.CA0      = createSlider(panel,'CA0',0,1e7,params.CA0,@updatePlot,ypos); ypos = ypos - 0.04;
sliders.tend    = createSlider(panel,'Years end',0.1,10,params.tend,@updatePlot,ypos);ypos = ypos - 0.04;
sliders.xmin    = createSlider(panel,'xmin',0,1,params.xmin,@updatePlot,ypos);ypos = ypos - 0.04;
%sliders.exh    = createSlider(panel,'exh',1e6,1e10,params.exh,@updatePlot,ypos);ypos = ypos - 0.04;
sliders.mrdPct  = createSlider(panel,'MRD (%)', 0.001, 1, params.mrdPct, @updatePlot, ypos);

C_blue=[55,126,184]/255; C_green=[77,175,74]/255; C_brown=[166,86,40]/255; C_purple=[152,78,163]/255;
[H_mesh, D_mesh] = meshgrid(linspace(0,1,100), linspace(0,1,100));
ImgL = zeros(100,100,3); ImgB = zeros(100,100,3);
for i=1:100
    for j=1:100
        w = D_mesh(i,j)^0.4;
        ImgL(i,j,:) = ((1-H_mesh(i,j))*C_green + H_mesh(i,j)*C_brown)*w + [1 1 1]*(1-w);
        ImgB(i,j,:) = ((1-H_mesh(i,j))*C_purple + H_mesh(i,j)*C_blue)*w + [1 1 1]*(1-w);
    end
end

updatePlot();

function updatePlot(~,~)
    flds = fieldnames(sliders);
    for idx = 1:numel(flds), params.(flds{idx}) = get(sliders.(flds{idx}),'Value'); end
    contents = get(popScenario,'String');
    params.scenario = contents{get(popScenario,'Value')};
    
    [t, CA, CM, L, B] = solveSystem(params);
    t = t / 365; 
    xspan = linspace(0, 1, params.Nx); [X, T] = meshgrid(xspan, t);
    
    cla(ax3);
    plot(ax3, t, CA, 'Color', [1 0.5 0], 'LineWidth', 2); hold(ax3, 'on');
    plot(ax3, t, CM, '--', 'Color', [1 0.5 0], 'LineWidth', 2);
    xlabel(ax3, 'Years'); ylabel(ax3, 'CAR-T Cells'); grid(ax3, 'on');
    legend(ax3,["Activated CAR-T","Memory CAR-T"]);
    
    cla(ax1);
    RGB_L = zeros(size(L,1), size(L,2), 3); maxL = max(L(:)) + 1e-12;
    for jj = 1:params.Nx
        base = (1-xspan(jj))*C_green + xspan(jj)*C_brown;
        for ii = 1:length(t), ...%ww = (L(ii,jj)/maxL)^0.4
                ww=(min(L(ii,jj), 1e11) / 1e11)^0.4; RGB_L(ii,jj,:) = base*ww + [1 1 1]*(1-ww); end
    end
    surf(ax1, X, T, L, 'CData', RGB_L, 'FaceColor', 'interp', 'EdgeColor', 'none');
    set(ax1, 'YDir', 'reverse', 'YTick', 0:1:params.tend); % Force consistent ticks
    view(ax1, 30, 45); grid(ax1, 'on'); 
    xlabel(ax1, 'Antigen','Rotation',-5); ylabel(ax1, 'Years','Rotation',30); title(ax1, 'L Cells');
    zlim(ax1, [0 5e10]); % Fixed for L Cells scale

    cla(ax2);
    RGB_B = zeros(size(B,1), size(B,2), 3); maxB = max(B(:)) + 1e-12;
    for jj = 1:params.Nx
        base = (1-xspan(jj))*C_purple + xspan(jj)*C_blue;
        for ii = 1:length(t), ww = (B(ii,jj)/maxB)^0.4; RGB_B(ii,jj,:) = base*ww + [1 1 1]*(1-ww); end
    end
    surf(ax2, X, T, B, 'CData', RGB_B, 'FaceColor', 'interp', 'EdgeColor', 'none');
    set(ax2, 'YDir', 'reverse', 'YTick', 0:1:params.tend); % Force consistent ticks
    view(ax2, 30, 45); grid(ax2, 'on');
    xlabel(ax2, 'Antigen','Rotation',-5); ylabel(ax2, 'Years','Rotation',30); title(ax2, 'B Cells');
    zlim(ax2, [0 5e10]); % Fixed for B Cells scale
    
    image(axLeg1, [0 1], [1e9 maxL], ImgL); set(axLeg1, 'YDir', 'normal', 'FontSize', 7);
    xlabel(axLeg1, 'Antigen'); ylabel(axLeg1, 'Cells (L)');
    t1 = title(axLeg1, 'Legend (L cells)');
    p1 = get(t1, 'Position'); set(t1, 'Position', [p1(1), p1(2) + (maxL * 0.15), p1(3)]); 
    
    image(axLeg2, [0 1], [1e9 maxB], ImgB); set(axLeg2, 'YDir', 'normal', 'FontSize', 7);
    xlabel(axLeg2, 'Antigen'); ylabel(axLeg2, 'Cells (B)');
    t2 = title(axLeg2, 'Legend (B cells)');
    p2 = get(t2, 'Position'); set(t2, 'Position', [p2(1), p2(2) + (maxB * 0.15), p2(3)]);

    % MRD Threshold Value
mrdVal = (params.mrdPct / 100) * 1e12; 
colorred='#E23F44';
% Draw the line on the Back Walls of ax1
hold(ax1, 'on');
% Line along the Antigen wall (X)
plot3(ax1, [0 1], [0 0], [mrdVal mrdVal], '--', 'Color',colorred, 'LineWidth', 1); 
% Line along the Years wall (Y)
plot3(ax1, [1 1], [0 params.tend], [mrdVal mrdVal], '--','Color',colorred,  'LineWidth', 1);

% Repeat for ax2 if needed
hold(ax2, 'on');
plot3(ax2, [0 1], [0 0], [mrdVal mrdVal], '--', 'Color',colorred, 'LineWidth', 1);
plot3(ax2, [1 1], [0 params.tend], [mrdVal mrdVal], '--', 'Color',colorred, 'LineWidth', 1);

% This marks the specific tick as "MRD"

% --- For L Cells (ax1) ---
hold(ax1, 'on');
% This draws a single contour line at exactly the mrdVal height
[~, ~] = contour3(ax1, X, T, L, [mrdVal mrdVal], '--','Color',colorred,  'LineWidth', 1);

% --- For B Cells (ax2) ---
%hold(ax2, 'on');
%[~, ~] = contour3(ax2, X, T, B, [mrdVal mrdVal], '--','Color',colorred,  'LineWidth', 1);
end

function exportData(~,~)
    
    selPath = uigetdir(pwd, 'Select Folder to Save Exported Data');
    if selPath == 0, return; end 
    flds = fieldnames(sliders);
    for idx = 1:numel(flds), params.(flds{idx}) = get(sliders.(flds{idx}),'Value'); end
    contents = get(popScenario,'String');
    params.scenario = contents{get(popScenario,'Value')};
    
    paramTable = struct2table(params);
    writetable(paramTable, fullfile(selPath, 'Model_Parameters.csv'));
    
    exportAxisWithLegend(ax1, axLeg1, fullfile(selPath, 'L_Cells_Plot.pdf'), true);
    exportAxisWithLegend(ax2, axLeg2, fullfile(selPath, 'B_Cells_Plot.pdf'), true);
    exportAxisWithLegend(ax3, [], fullfile(selPath, 'CART_Kinetics_Plot.pdf'), false);
    msgbox('Data and PDFs exported successfully', 'Done');
end

    function exportAxisWithLegend(sourceAx, legendAx, fileName, is3D)
    tempFig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [0 0 1000 700], 'Color', 'w');
    newMain = copyobj(sourceAx, tempFig);
    set(newMain, 'Units', 'normalized', 'Position', [0.1 0.15 0.65 0.75]);
    
    if ~isempty(legendAx)
        newLeg = copyobj(legendAx, tempFig);
        set(newLeg, 'Units', 'normalized', 'Position', [0.82 0.35 0.05 0.35]);
    end
    
    % Use the 'tend' slider value directly for consistent axis limits and ticks
    currentTend = get(sliders.tend, 'Value');
    
    if is3D
        % Set identical limits and tick spacing for consistency
        set(newMain, 'YLim', [0 currentTend], 'YTick', 0:1:currentTend); 
        view(newMain, 30, 45); 
    else
        set(newMain, 'XLim', [0 currentTend]); % For ax3
        view(newMain, 2); 
    end
    
    grid(newMain, 'on');
    % Hybrid export: Vector text/axes with a High-Res raster surface to save space
    exportgraphics(tempFig, fileName, 'ContentType', 'auto', 'Resolution', 300);
    close(tempFig);
end

function [t,CA,CM,L,B] = solveSystem(p)
    tend_days = p.tend * 365;
    tspan = linspace(0, tend_days, min(round(tend_days)+1, 500)); 
    xspan = linspace(0,1,p.Nx);
    dx = xspan(2) - xspan(1);      % Proper spacing for integration
    switch p.scenario
        case 'Single Clone (0.6)', ICL = gaussian(xspan, 0.6, p.sigma)*p.L0;
        case 'Single Clone (0.4)', ICL = gaussian(xspan, 0.4, p.sigma)*p.L0;
        case 'Two Clones (0.4 & 0.6)', ICL = (gaussian(xspan,0.4,p.sigma)*0.5 + gaussian(xspan,0.6,p.sigma)*0.5)*p.L0;
    end
  
    ICB = gaussian(xspan, 0.2, p.sigma)*p.B0;
    IC = [p.CA0; p.CM0; ICL(:); ICB(:)];
    [t,res] = ode45(@(t,e) odefun_core(t,e,dx,p), tspan, IC, odeset('RelTol',1e-5));
    CA=res(:,1); CM=res(:,2); L=res(:,3:p.Nx+2); B=res(:,p.Nx+3:end);
end

function ne = odefun_core(~,e,dx,p)
    CA=e(1); CM=e(2); L=e(3:p.Nx+2); B=e(p.Nx+3:end); 
    L(L<1)=0; B(B<1)=0; x=linspace(0, 1, p.Nx)';
    mask = (x >= p.xmin);
    IL=sum(L(mask))*dx; IB=sum(B(mask))*dx; IN=IL+IB; F=IN/(p.k+IN);
    Sx = p.s./(1+exp(p.m*(x-p.x0)));
    hB_vec = (gaussian(x,0.25,p.sigma)*0.1 + gaussian(x,0.5,p.sigma)*0.65 + gaussian(x,0.75,p.sigma)*0.25)*p.BMax;
    dL = p.rhoL.*L.*(1-IL/p.LMax) - p.alpha.*(x./(x+p.h+1e-5)).*L.*CA;
    dB = p.BStem*Sx + p.rhoB.*B.*(1-(B+IL)./(hB_vec+IL)) - p.alpha.*(x./(x+p.h+1e-5)).*B.*CA - B/p.tauB;
    dCA = F*(p.rhoCA*CA + p.gammaMA*CM) - (1/p.tauCA)*CA - (1-F)*p.gammaAM*CA;...
    dCM = (1-F)*p.gammaAM*CA - (1/p.tauCM)*CM - F*p.gammaMA*CM;
    ne = [dCA; dCM; dL; dB]; 
end

function g = gaussian(x, m, s), g = exp(-((x-m)/s).^2/2)./(s*sqrt(2*pi)); end

function sl = createSlider(p,l,minv,maxv,initv,cb,yp)
    uicontrol('Parent',p,'Style','text','String',l,'Units','normalized',...
        'Position',[0.02 yp 0.25 0.035],'FontSize',7,'HorizontalAlignment','right');
    sl = uicontrol('Parent',p,'Style','slider','Min',minv,'Max',maxv,'Value',initv,...
        'Units','normalized','Position',[0.30 yp 0.40 0.03],'Callback',cb);
    et = uicontrol('Parent',p,'Style','edit','String',num2str(initv,'%g'),...
        'Units','normalized','Position',[0.72 yp 0.25 0.04],'FontSize',7,...
        'Callback', @(src, event) editBoxCallback(src, sl, minv, maxv, cb));
    addlistener(sl, 'Value', 'PostSet', @(s,e) set(et, 'String', num2str(get(sl,'Value'),'%g')));
end

function editBoxCallback(src, slider, minv, maxv, updatePlotCB)
    val = str2double(get(src, 'String'));
    if ~isnan(val)
        val = max(minv, min(val, maxv)); set(slider, 'Value', val);
        set(src, 'String', num2str(val, '%g')); updatePlotCB(); 
    else
        set(src, 'String', num2str(get(slider, 'Value'), '%g'));
    end
end

function resetSliders(~,~)
    fds = fieldnames(sliders); 
    for ii = 1:numel(fds), set(sliders.(fds{ii}), 'Value', defaultParams.(fds{ii})); end
    updatePlot(); 
end
end