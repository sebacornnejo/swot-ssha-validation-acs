%% =======================================================================
%  SWOT vs IN SITU SPECTRAL COHERENCE SYNTHESIS - 1 ROW x 3 COLUMNS
%  (a) gamma^2 heatmap   (b) Weighted phase lag (h)   (c) Variance ratio %
%  Stations ordered N->S. Panel positions set manually with axes().
%  Reads : ./Estadisticas_Spectra/Master_TEST1_Spectral_Coherence.xlsx
%  Saves : ./coherence_insituvsswot at 600 dpi.
%  =======================================================================
clear; close all; clc;

file_xlsx = './Estadisticas_Spectra/Master_TEST1_Spectral_Coherence.xlsx';

sheets    = {'G. San Matias','P. Deseado','Lander','Carina','Vega-Pleyade'};
sheet_lbl = {'G. San Matias','P. Deseado','Lander','Carina','Vega-Pleyade'};
corrs     = {'Trad','CalVal','MSAS','SIROCCO','ERA5','IB-MSLP'};

cmap_dac = [196 196 196;
            126 222 214;
            195 147 255;
            253 154  95;
            250 223  90;
            246  96  97] / 255;

FONT  = 'Times New Roman';
FS    = 20;
IB_ms = 12.35;   % 95% of 13

nS = numel(sheets); nC = numel(corrs);
[COH, PW, LW, RV] = deal(nan(nS, nC));

for i = 1:nS
    T      = readtable(file_xlsx, 'Sheet', sheets{i}, 'ReadVariableNames', true);
    cnames = strtrim(string(T.Correccion));
    for j = 1:nC
        k = find(strcmpi(cnames, corrs{j}), 1);
        if isempty(k), continue; end
        COH(i,j) = T.Coherencia_Media_2_10d(k);
        PW(i,j)  = T.Fase_Ponderada_deg(k);
        LW(i,j)  = T.Desfase_Ponderado_Horas(k);
        RV(i,j)  = T.Ratio_Varianza_SWOT_Insitu(k) * 100;
    end
end

off = linspace(-0.34, 0.34, nC);

%% --- Canvas and manual positions (no tiledlayout) ----------------------
scr = get(groot, 'ScreenSize');
fig = figure('Visible','off','Color','w','Units','pixels', ...
             'Position',[1 1 round(scr(3)*1.30) scr(4)]);

% Normalized layout: 3 columns. Each panel is shrunk to 90% (SHRINK) of its
% available width, keeping the height and centring within its cell.
margL = 0.085; margR = 0.045; gap = 0.075;   % margins and gap between panels
panB  = 0.135; panH = 0.80;                    % panel base and height
SHRINK = 0.90;                                 % panel width at 90%

totalW = 1 - margL - margR - 2*gap;            % usable width for 3 panels
cellW  = totalW / 3;                           % cell width per panel
panW   = cellW * SHRINK;                       % actual panel width (90%)
padW   = (cellW - panW) / 2;                   % padding on each side

xpos = margL + (0:2)*(cellW + gap) + padW;     % x of each panel
posA = [xpos(1) panB panW panH];
posB = [xpos(2) panB panW panH];
posC = [xpos(3) panB panW panH];

%% ---------- (a) gamma^2 heatmap ----------------------------------------
ax1 = axes('Parent',fig,'Position',posA); hold(ax1,'on');
imagesc(ax1, COH, [0.15 1]);
colormap(ax1, slanCM_fallback());
cb = colorbar(ax1);
cb.Label.String  = '\gamma^2 (2-10 d)';
cb.Label.FontName = FONT; cb.FontName = FONT; cb.FontSize = FS;
set(ax1,'YDir','reverse','XTick',1:nC,'XTickLabel',corrs, ...
        'YTick',1:nS,'YTickLabel',sheet_lbl,'TickLabelInterpreter','none', ...
        'FontName',FONT,'FontSize',FS,'Layer','top');
xlim(ax1,[0.5 nC+0.5]); ylim(ax1,[0.5 nS+0.5]);
for i = 1:nS
    for j = 1:nC
        if isnan(COH(i,j)), continue; end
        nrm    = (COH(i,j) - 0.15) / (1 - 0.15);
        txtcol = 'k'; if nrm < 0.5, txtcol = 'w'; end
        text(ax1, j, i, sprintf('%.2f', COH(i,j)), ...
             'HorizontalAlignment','center','FontName',FONT, ...
             'FontSize',FS,'FontWeight','bold','Color',txtcol);
    end
end
xline(ax1, nC-0.5, 'k-', 'LineWidth', 1.4);
title(ax1,'(a)','FontName',FONT,'FontSize',FS,'FontWeight','normal');

%% ---------- (b) Weighted phase lag (h) ---------------------------------
ax2 = axes('Parent',fig,'Position',posB); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
for i = 1:nS
    for j = 1:nC
        y = i + off(j);
        plot(ax2, [0 LW(i,j)], [y y], '-', 'Color', cmap_dac(j,:), 'LineWidth', 1.8);
        if j < nC
            plot(ax2, LW(i,j), y, 'o', 'MarkerFaceColor', cmap_dac(j,:), ...
                 'MarkerEdgeColor','k', 'MarkerSize', 10);
        else
            plot(ax2, LW(i,j), y, 'd', 'MarkerFaceColor', cmap_dac(j,:), ...
                 'MarkerEdgeColor','k', 'MarkerSize', IB_ms);
        end
    end
end
xtk2 = 0:5:15;
set(ax2,'YDir','reverse','YTick',1:nS,'YTickLabel',{}, ...
        'XTick',xtk2,'XTickLabel',arrayfun(@num2str,xtk2,'UniformOutput',false), ...
        'TickLabelInterpreter','none','FontName',FONT,'FontSize',FS);
xlim(ax2,[0 20]); ylim(ax2,[0.5 nS+0.5]);
xlabel(ax2,'Weighted phase lag (h)','FontName',FONT,'FontSize',FS);
title(ax2,'(b)','FontName',FONT,'FontSize',FS,'FontWeight','normal');

% Legend in (b)
hL = gobjects(1,nC);
for j = 1:nC
    if j < nC
        hL(j) = plot(ax2, nan, nan, 'o', 'MarkerFaceColor', cmap_dac(j,:), ...
                     'MarkerEdgeColor','k', 'MarkerSize', 10);
    else
        hL(j) = plot(ax2, nan, nan, 'd', 'MarkerFaceColor', cmap_dac(j,:), ...
                     'MarkerEdgeColor','k', 'MarkerSize', IB_ms);
    end
end
lg = legend(hL, corrs, 'Location','southeast', 'NumColumns', 2);
set(lg,'FontName',FONT,'FontSize',16,'Box','on');

%% ---------- (c) Variance ratio (%) -------------------------------------
ax3 = axes('Parent',fig,'Position',posC); hold(ax3,'on'); box(ax3,'on'); grid(ax3,'on');
patch(ax3,[90 110 110 90],[0.5 0.5 nS+0.5 nS+0.5],[0.92 0.92 0.92],'EdgeColor','none');
xline(ax3, 100, 'k-', 'LineWidth', 1);
for i = 1:nS
    for j = 1:nC
        x = RV(i,j); y = i + off(j);
        plot(ax3, [100 x], [y y], '-', 'Color', cmap_dac(j,:), 'LineWidth', 1.8);
        if j < nC
            plot(ax3, x, y, 'o', 'MarkerFaceColor', cmap_dac(j,:), ...
                 'MarkerEdgeColor','k', 'MarkerSize', 10);
        else
            plot(ax3, x, y, 'd', 'MarkerFaceColor', cmap_dac(j,:), ...
                 'MarkerEdgeColor','k', 'MarkerSize', IB_ms);
        end
    end
end
xtk3 = 0:25:125;
set(ax3,'YDir','reverse','YTick',1:nS,'YTickLabel',{}, ...
        'XTick',xtk3,'XTickLabel',arrayfun(@num2str,xtk3,'UniformOutput',false), ...
        'TickLabelInterpreter','none','FontName',FONT,'FontSize',FS);
xlim(ax3,[0 140]); ylim(ax3,[0.5 nS+0.5]);
xlabel(ax3,'var(SWOT)/var(in situ) [%]','FontName',FONT,'FontSize',FS);
title(ax3,'(c)','FontName',FONT,'FontSize',FS,'FontWeight','normal');

%% ------------------------------ EXPORT ---------------------------------
outDir = 'coherence_insituvsswot';
if ~exist(outDir,'dir'), mkdir(outDir); end
exportgraphics(fig, fullfile(outDir,'Fig04.png'), 'Resolution', 600);
fprintf('Figure saved in: %s\n', fullfile(outDir,'Fig04.png'));
close(fig);

%% -----------------------------------------------------------------------
function cm = slanCM_fallback()
anchors = [ 68  1 84;  59 82 139;  33 145 140;  94 201  98; 253 231 37] / 255;
x  = linspace(0, 1, size(anchors,1));
xi = linspace(0, 1, 256);
cm = [interp1(x,anchors(:,1),xi)' interp1(x,anchors(:,2),xi)' ...
      interp1(x,anchors(:,3),xi)'];
end
