%% ========================================================================
%  R^2 bubble maps, SWOT vs in situ, 3 panels:
%     (a) no tide/DAC correction   (b) tide + IB correction
%     (c) the 5 DAC series (Trad/CalVal/MSAS/SIROCCO/ERA5), 5 mini-bubbles
%         per station anchored in open ocean with a leader line to the real
%         location (avoids overlap when stacked on the geographic position).
%  Bubble colour and size scale with R^2 [0-100] (OrRd_09 palette).
% ========================================================================

clear all; close all; clc;

try
    addpath("../../m_map1.4f");
catch
    warning('m_map not found at ../../m_map1.4f');
end

outputDir = './r2_maps_3panels_mslp';
if ~exist(outputDir, 'dir')
   mkdir(outputDir);
end

%% ----------------------- Block 1: Data ---------------------------------
% Stations
stationNames = {'SBE37 - Lander', 'SBE26 - G. San Matias', 'TG - P. Deseado', ...
                'Radar - Carina', 'Radar - Vega-Pleyade'};
stationLon = [-66.67, -63.7803, -65.917, -67.21960, -67.74600];
stationLat = [-52.14, -41.1822, -47.75, -52.75718, -53.29900];
nSta = numel(stationLon);

% --- Panel (a): in situ / no tide or DAC correction ---
r2_noCorr = [95.38, 97.85, 92.12, 97.38, 96.39];

% --- Panel (b): IB ---
r2_IB     = [89.98035432,54.76748837,60.96758975,92.73367978,90.61093102];

% --- Panel (c): the 5 DAC series. Rows = stations, columns = DAC ---
%   columns: Trad | CalVal | MSAS | SIROCCO | ERA5
dacNames = {'Trad','CalVal','MSAS','SIROCCO','ERA5'};
nDAC = numel(dacNames);
% Each row corresponds to one station (same order as stationLon)
r2_DAC = [ ...
   51.77, 60.83, 65.06, 52.17, 50.17;   % Lander
   61.39, 62.92, 32.11, 34.01, 39.11;   % SBE26
   11.16, 34.92, 23.01,  4.67, 12.54;   % P. Deseado
   70.94, 75.95, 76.39, 71.09, 69.05;   % Carina
   70.37, 75.16, 73.19, 75.46, 68.16];  % Vega

%% ------------------- Block 2: Plot parameters --------------------------
map_lon_limits = [-69.25 -61.27];
map_lat_limits = [-53.99 -39.02];

% --- OrRd_09 colormap ---
try
    aux1 = colormap_cpt('OrRd_09');
catch
    % 9-class OrRd (ColorBrewer) embedded as fallback
    aux1 = [255 247 236; 254 232 200; 253 212 158; 253 187 132; ...
            252 141  89; 239 101  72; 215  48  31; 179   0   0; ...
            127   0   0] / 255;
end
aux2 = linspace(1,100,length(aux1(:,1))); aux3 = 1:100;
R = interp1(aux2,aux1(:,1),aux3,'linear');
G = interp1(aux2,aux1(:,2),aux3,'linear');
B = interp1(aux2,aux1(:,3),aux3,'linear');
cmapt = [R' G' B'];
cmapt(cmapt>1)=1; cmapt(cmapt<0)=0;
clear aux1 aux2 aux3 R G B;

% Map an R^2 value [0-100] to a palette colour
val2color = @(v) cmapt(min(100,max(1,round(v))), :);

% Bubble size scaling
min_bubble_size = 50;
max_bubble_size = 1500;
sz = @(v) min_bubble_size + (v/100)*(max_bubble_size - min_bubble_size);

% Pass 7 tracks
try
    Pass7 = kmz2struct('Pass 7.kmz');
    havePass = true;
catch
    warning('Pass 7.kmz not found. Tracks omitted.');
    havePass = false;
end

%% =======================================================================
%  PANEL (a): No Correction
%  =======================================================================
figure('Units','pixels','Position',[0 0 1920 1080],'Visible','off');
hA = gca;
m_proj('mercator','long',map_lon_limits,'lat',map_lat_limits);
hold on; set(gcf,'Renderer','painters');
if havePass
    m_plot(Pass7(2).Lon, Pass7(2).Lat, '-', 'linewidth',1,'Color',[.25 .25 .25]);
    m_plot(Pass7(1).Lon-0.135, Pass7(1).Lat, '-', 'linewidth',1,'Color',[.25 .25 .25]);
    m_plot(Pass7(1).Lon+0.135, Pass7(1).Lat, '-', 'linewidth',1,'Color',[.25 .25 .25]);
    m_plot(Pass7(1).Lon, Pass7(1).Lat, '-.', 'linewidth',1,'color',[248,161,90]/255);
end
m_gshhs_f('patch',[.83 .83 .83],'edgecolor','k');
m_grid('linestyle',':','FontName','Times New Roman','FontSize',20,'xtick',-68:2:-62);

hs = m_scatter(stationLon, stationLat, sz(r2_noCorr), r2_noCorr, ...
    'filled','MarkerEdgeColor',[.25 .25 .25],'LineWidth',1.5);
set(hs,'Clipping','off');
% labels: south (Lander/Carina/Vega) to the right, north (SBE26/PDeseado)
% above; latNudge lowers P. Deseado only
nudge_a = [0, 0, -0.18, 0, 0];
htxt = labelText(stationLon, stationLat, r2_noCorr, 0.70, 0.55, nudge_a);
colormap(hA, cmapt); caxis(hA,[0 100]);
reorderLayers(hA, hs, htxt);
title('(a)','FontName','Times New Roman','FontSize',20);
set(gca,'FontName','Times New Roman','FontSize',20);
print('-dpng', fullfile(outputDir,'Fig05a_NoCorrection.png'), '-r600');
hold off;

%% =======================================================================
%  PANEL (b): IB
%  =======================================================================
figure('Units','pixels','Position',[0 0 1920 1080],'Visible','off');
hB = gca;
m_proj('mercator','long',map_lon_limits,'lat',map_lat_limits);
hold on; set(gcf,'Renderer','painters');
if havePass
    m_plot(Pass7(2).Lon, Pass7(2).Lat, '-', 'linewidth',1,'Color',[.25 .25 .25]);
    m_plot(Pass7(1).Lon-0.135, Pass7(1).Lat, '-', 'linewidth',1,'Color',[.25 .25 .25]);
    m_plot(Pass7(1).Lon+0.135, Pass7(1).Lat, '-', 'linewidth',1,'Color',[.25 .25 .25]);
    m_plot(Pass7(1).Lon, Pass7(1).Lat, '-.', 'linewidth',1,'color',[248,161,90]/255);
end
m_gshhs_f('patch',[.83 .83 .83],'edgecolor','k');
m_grid('linestyle',':','FontName','Times New Roman','FontSize',20,'xtick',-68:2:-62);

hs = m_scatter(stationLon, stationLat, sz(r2_IB), r2_IB, ...
    'filled','MarkerEdgeColor',[.25 .25 .25],'LineWidth',1.5);
set(hs,'Clipping','off');
htxt = labelText(stationLon, stationLat, r2_IB, 0.70, 0.40);
colormap(hB, cmapt); caxis(hB,[0 100]);
reorderLayers(hB, hs, htxt);
title('(b)','FontName','Times New Roman','FontSize',20);
set(gca,'FontName','Times New Roman','FontSize',20);
print('-dpng', fullfile(outputDir,'Fig05b_IB.png'), '-r600');
hold off;

%% =======================================================================
%  PANEL (c): the 5 DAC series. Per station, a row of 5 mini-bubbles in
%  open ocean connected by a leader line to the real location; fixed column
%  order (Trad, CalVal, MSAS, SIROCCO, ERA5), labelled once.
%  =======================================================================
figure('Units','pixels','Position',[0 0 1920 1080],'Visible','off');
hC = gca;
m_proj('mercator','long',map_lon_limits,'lat',map_lat_limits);
hold on; set(gcf,'Renderer','painters');
if havePass
    m_plot(Pass7(2).Lon, Pass7(2).Lat, '-', 'linewidth',1,'Color',[.25 .25 .25]);
    m_plot(Pass7(1).Lon-0.135, Pass7(1).Lat, '-', 'linewidth',1,'Color',[.25 .25 .25]);
    m_plot(Pass7(1).Lon+0.135, Pass7(1).Lat, '-', 'linewidth',1,'Color',[.25 .25 .25]);
    m_plot(Pass7(1).Lon, Pass7(1).Lat, '-.', 'linewidth',1,'color',[248,161,90]/255);
end
m_gshhs_f('patch',[.83 .83 .83],'edgecolor','k');
m_grid('linestyle',':','FontName','Times New Roman','FontSize',20,'xtick',-68:2:-62);
colormap(hC, cmapt); caxis(hC,[0 100]);

% Anchor (centre of the 5-bubble row) per station, in open ocean.
% Row order = stationLon (Lander, SBE26, PDeseado, Carina, Vega)
anchorLon = [-64.30, -63.45, -63.60, -64.30, -64.30];
deltaSouth = 0.03;  % shifts the 3 southern anchors up as a block
anchorLat = [-52.10+deltaSouth, -41.18, -47.75, -52.85+deltaSouth, -53.55+deltaSouth];

% Horizontal spacing between mini-bubbles of the same station (deg lon)
dLon = 0.62;
% Size scale unified with (a)/(b): same sz(v) reduced by a constant area
% factor, so the visual comparison across panels remains valid
k_scale = 0.2535;
mini_sz = @(v) sz(v) * k_scale;

% Column offsets relative to the anchor centre
colOff = ((1:nDAC) - (nDAC+1)/2) * dLon;   % symmetric around 0

% Vertical offset of the value below its bubble (deg lat)
dLatVal = 0.16;
% Fine per-station vertical nudge [Lander, SanMatias, PDeseado, Carina, Vega]
valNudge_c = [0, -0.07, 0, 0, 0];

hScatterAll = [];
hValTxt = [];

for s = 1:nSta
    % leader line from the real location to the centre of the bubble row
    m_plot([stationLon(s) anchorLon(s)], [stationLat(s) anchorLat(s)], ...
        '-', 'Color',[.4 .4 .4], 'LineWidth',1.0);
    % real station location: filled black square
    m_plot(stationLon(s), stationLat(s), 's', 'MarkerSize',8, ...
        'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',1.0);

    for d = 1:nDAC
        blon = anchorLon(s) + colOff(d);
        blat = anchorLat(s);
        v = r2_DAC(s,d);
        h = m_scatter(blon, blat, mini_sz(v), v, ...
            'filled','MarkerEdgeColor',[.25 .25 .25],'LineWidth',1.0);
        set(h,'Clipping','off');
        hScatterAll(end+1) = h;
        % numeric value below each mini-bubble
        hValTxt(end+1) = boxedText(blon, blat-dLatVal+valNudge_c(s), sprintf('%.0f', v), 15, 'center', 'top');
    end
end

% --- Mini-legend: which column is which DAC (top left, over land) ---
legLon0 = -68.95; legLat0 = -39.8; legStep = 0.45;
for d = 1:nDAC
    haloText(legLon0, legLat0 - (d-1)*legStep, sprintf('%d = %s', d, dacNames{d}), ...
        15, 'left', 'middle', 0.018, 0.012);
end
haloText(legLon0, legLat0 + legStep*0.9, 'Column order (left to right):', ...
    15, 'left', 'middle', 0.018, 0.012, 'bold');

% Layer order: bubbles in front, then value texts, then everything else
ch = get(hC,'Children');
isBub = ismember(ch, hScatterAll(:));
isTxt = ismember(ch, hValTxt(:));
isRest = ~(isBub | isTxt);
set(hC,'Children',[ch(isBub); ch(isTxt); ch(isRest)]);

haloTitle('(c)', 20);
set(gca,'FontName','Times New Roman','FontSize',20);
print('-dpng', fullfile(outputDir,'Fig05c_DAC.png'), '-r600');
hold off;

fprintf('Panels Fig05a/b/c saved in %s\n', outputDir);

%% ---------------------- Local functions --------------------------------
function ht = labelText(lon, lat, vals, dLonLbl, dLatLbl, latNudge)
% R^2 value next to each bubble (georeferenced m_text) with a white
% double-text halo for contrast. South (Lander/Carina/Vega): label to the
% right; north (SBE26/P.Deseado): above.
% dLonLbl: eastward offset of the southern labels (deg lon).
% dLatLbl: upward offset of the northern labels (deg lat).
% latNudge: individual vertical nudge per station (deg).
    if nargin < 4 || isempty(dLonLbl),  dLonLbl  = 0.45; end
    if nargin < 5 || isempty(dLatLbl),  dLatLbl  = 0.40; end
    if nargin < 6 || isempty(latNudge), latNudge = zeros(1, numel(lon)); end
    south = [1 4 5];
    ht = gobjects(1, numel(lon));
    for k = 1:numel(lon)
        if ismember(k, south)
            ht(k) = haloText(lon(k)+dLonLbl, lat(k)+latNudge(k), sprintf('%.1f', vals(k)), ...
                18, 'left', 'middle', 0.018, 0.012);
        else
            ht(k) = haloText(lon(k), lat(k)+dLatLbl+latNudge(k), sprintf('%.1f', vals(k)), ...
                18, 'center', 'bottom', 0.018, 0.012);
        end
    end
end

function ht = boxedText(lon, lat, str, fs, ha, va)
% Panel (c) text: smaller halo radius (values sit close together over the
% ocean) with the same thickness as in panels a/b.
    ht = haloText(lon, lat, str, fs, ha, va, 0.010, 0.012);
end

function htFront = haloText(lon, lat, str, fs, ha, va, haloOff, haloThickness, fontWeight)
% Draws 'str' with a continuous white halo (concentric rings of white
% copies) and a black copy on top, using georeferenced m_text.
% haloOff: outer halo radius (deg). haloThickness: halo thickness (deg).
% Returns the handle of the front (black) copy.
    if nargin < 7 || isempty(haloOff),       haloOff       = 0.018; end
    if nargin < 8 || isempty(haloThickness), haloThickness = 0.006; end
    if nargin < 9 || isempty(fontWeight),    fontWeight    = 'normal'; end

    ringStep = 0.004;  % radial spacing between concentric rings
    nRings = max(1, round(haloThickness/ringStep));
    radii = linspace(haloOff - haloThickness, haloOff, nRings);
    radii = radii(radii > 0);

    for r = radii
        % more directions on larger rings (continuous outline)
        nDir = max(8, round(16*r/0.018));
        ang = linspace(0, 360, nDir+1); ang = ang(1:end-1);
        for a = ang
            dlon = r*cosd(a);
            dlat = r*sind(a);
            m_text(lon+dlon, lat+dlat, str, ...
                'FontName','Times New Roman','FontSize',fs,'FontWeight',fontWeight, ...
                'HorizontalAlignment',ha,'VerticalAlignment',va, ...
                'Color','w');
        end
    end

    htFront = m_text(lon, lat, str, ...
        'FontName','Times New Roman','FontSize',fs,'FontWeight',fontWeight, ...
        'HorizontalAlignment',ha,'VerticalAlignment',va, ...
        'Color','k');
end

function reorderLayers(ax, hBubbles, hTexts)
% Brings the bubbles to the front of the axes layer stack.
    ch = get(ax,'Children');
    isB = ismember(ch, hBubbles(:));
    set(ax,'Children',[ch(isB); ch(~isB)]);
end

function haloTitle(str, fs)
% title() with the same white halo as the labels; uses text() in normalized
% axes coordinates (not m_text) since it is not a map location.
    ax = gca;
    haloOff = 0.006;        % halo radius, normalized units
    haloThickness = 0.004;  % halo thickness
    ringStep = 0.0015;
    nRings = max(1, round(haloThickness/ringStep));
    radii = linspace(haloOff - haloThickness, haloOff, nRings);
    radii = radii(radii > 0);

    % Standard title() position: centred at x=0.5, just above the axes
    tx = 0.5; ty = 1.025;

    for r = radii
        nDir = max(8, round(16*r/haloOff));
        ang = linspace(0, 360, nDir+1); ang = ang(1:end-1);
        for a = ang
            dx = r*cosd(a); dy = r*sind(a);
            text(ax, tx+dx, ty+dy, str, 'Units','normalized', ...
                'FontName','Times New Roman','FontSize',fs,'FontWeight','bold', ...
                'HorizontalAlignment','center','VerticalAlignment','bottom', ...
                'Color','w');
        end
    end

    text(ax, tx, ty, str, 'Units','normalized', ...
        'FontName','Times New Roman','FontSize',fs,'FontWeight','bold', ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'Color','k');
end
