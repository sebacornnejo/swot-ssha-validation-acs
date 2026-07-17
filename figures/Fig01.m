clear all
close all
clc
addpath("../../m_map1.4f")

% Load SWOT L3 data
load('./swot_pata022v2_0_1.mat','swot');
swot022=swot; clear swot

load('./swot_pata007v2_0_1.mat','swot');
swot007=swot; clear swot

% Load SWOT L4 data
load('./swot_pataL4.mat','swot');

[lonalt,latalt]=meshgrid(swot.lon,swot.lat);

% Load conventional altimetry data
file_nc2='./DATA/ClassicAltimetry/cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.25deg_P1D.nc';
tp=[];
tp.lon=double(ncread(file_nc2,'longitude'));
tp.lat=double(ncread(file_nc2,'latitude'));
tp.adt=ncread(file_nc2,'adt');
tp.sla=ncread(file_nc2,'sla');
tp.t=datenum(1970,1,1,0,0,double(ncread(file_nc2,'time')));

[lontrad,lattrad]=meshgrid(tp.lon,tp.lat);

Pass7 = kmz2struct('Pass 7.kmz');
Pass22 = kmz2struct('Pass 22.kmz');

% Bathymetry
load('./New_SAtl_0-05_NN.mat','gebco_shn');

% Colormap (blue-white-red)
aux1 = [
    0,   0,   139;
    0,   0,   203;
    80,  80,  254;
    141, 141, 254;
    223, 223, 254;
    254, 254, 254;
    254, 254, 24;
    254, 146, 24;
    254, 75,  24;
    203, 0,   0;
    139, 0,   0
] / 255;

aux2=linspace(1,100,length(aux1(:,1)));
aux3=1:100;
R=interp1(aux2,aux1(:,1),aux3,'v5cubic');
G=interp1(aux2,aux1(:,2),aux3,'v5cubic');
B=interp1(aux2,aux1(:,3),aux3,'v5cubic');
cmap=[R' G' B'];
cmap(cmap>1)=1;
cmap(cmap<0)=0;
clear aux1 aux2 aux3 R G B

load swot007L3_processed_re_na.mat
swot_re007 = swot_re;
swot_na007 = swot_na; clear swot_na swot_re

load swot022L3_processed_re_na.mat
swot_re022 = swot_re;
swot_na022 = swot_na; clear swot_na swot_re


for tt=1:98;mask(tt,:,:)=swot007(tt).sshak;end;mask=squeeze(nanmean(mask,1)).*0+1;

ind = 93;
for i = 1:length(swot_na007.t)
    hxc = find(datenum(datestr(swot_na022.t(ind), 'dd-mmm-yyyy')) == datenum(datestr(swot_na007.t(i), 'dd-mmm-yyyy')));
    if isempty(hxc) == 0

        % --- Common data ---
        file_nc2 = './DATA/ClassicAltimetry/cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.25deg_P1D.nc';
        long = double(ncread(file_nc2, 'longitude'));
        lat = double(ncread(file_nc2, 'latitude'));
        minlon = min(long); maxlon = max(long);
        minlat = min(lat); maxlat = max(lat);

        hxc1 = find(datenum(datestr(swot_na022.t(ind), 'dd-mmm-yyyy')) == datenum(datestr(tp.t, 'dd-mmm-yyyy')));
        [lona, lata] = meshgrid(tp.lon, tp.lat);

        % --- Stations ---
        stationNames = {'SBE37 - Lander', 'SBE26 - G. San Matías', 'TG - P. Deseado', 'Radar - Hidra Norte', 'Radar - Carina', 'Radar - Vega-Pleyade'};
        stationLon = [-66.67, -63.7803, -65.917, -68.22095, -67.21960, -67.74600];
        stationLat = [-52.14, -41.1822, -47.75, -52.82061, -52.75718, -53.29900];

        stationColors = [
            [99, 49, 197];   % Purple
            [217, 95, 2];    % Orange
            [18, 191, 128];  % Teal
            [2, 235, 170];   % Turquoise
            [220, 50, 47];   % Red
            [50, 130, 250]   % Blue
        ] / 255;

        markerStyle = 's';
        markerSize = 6;
        markerFaceColor = 'w';
        lineWidth = 1.3;

        % --- Map: SWOT L3 (DAC-CalVal) ---
        figure;
        set(gcf, 'WindowState', 'maximized');
        m_proj('mercator', 'long', [-71.2 -56], 'lat', [-54 -38.4]);
        hold on;

        h1 = m_pcolorjw(lona, lata, squeeze(tp.sla(:,:,hxc1))');
        set(h1, 'FaceAlpha', 1); shading flat;
        h2 = m_pcolorjw(swot022(ind).lon-360, swot022(ind).lat, squeeze(swot_re022.sla(i,:,:)));
        set(h2, 'FaceAlpha', 1); shading flat;
        h3 = m_pcolorjw(swot007(i).lon-360, swot007(i).lat, squeeze(swot_re007.sla(i,:,:)).*mask);
        set(h3, 'FaceAlpha', 1); shading flat;

        colormap(cmap);
        hh=colorbar('southoutside');
        caxis([-.35 .35]);
        m_contour(gebco_shn.lon, gebco_shn.lat, gebco_shn.z', [-200 -200], 'k');
        m_plot(Pass22(2).Lon, Pass22(2).Lat, '-', 'markersize', 1.2, 'linewidth', 1, 'Color', [.25 .25 .25]);
        m_plot(Pass22(1).Lon-0.135, Pass22(1).Lat, '-', 'markersize', 1.2, 'linewidth', 1, 'Color', [.25 .25 .25]);
        m_plot(Pass22(1).Lon+0.135, Pass22(1).Lat, '-', 'markersize', 1.2, 'linewidth', 1, 'Color', [.25 .25 .25]);
        m_plot(Pass22(1).Lon, Pass22(1).Lat, '-.', 'linewidth', 1, 'color', [248, 161, 90]/255);
        m_plot(Pass7(2).Lon, Pass7(2).Lat, '-', 'markersize', 1.2, 'linewidth', 1, 'Color', [.25 .25 .25]);
        m_plot(Pass7(1).Lon-0.135, Pass7(1).Lat, '-', 'markersize', 1.2, 'linewidth', 1, 'Color', [.25 .25 .25]);
        m_plot(Pass7(1).Lon+0.135, Pass7(1).Lat, '-', 'markersize', 1.2, 'linewidth', 1, 'Color', [.25 .25 .25]);
        m_plot(Pass7(1).Lon, Pass7(1).Lat, '-.', 'linewidth', 1, 'color', [248, 161, 90]/255);

        for k = 1:length(stationNames)
            m_line(stationLon(k), stationLat(k), 'linestyle', 'none', 'marker', markerStyle, 'markersize', 10, ...
                   'markeredgecolor', stationColors(k,:), 'markerfacecolor', markerFaceColor, 'linewi', 2);
        end

        set(hh, ...
            'Ticks',[-.35 -.2 -.1 0 .1 .2 .35], ...
            'TickDirection','out', ...
            'FontName','Times New Roman', ...
            'FontSize',12);

        hh.Label.String   = 'Sea level anomaly [m]';
        hh.Label.FontName = 'Times New Roman';
        hh.Label.FontSize = 12;
        hh.Label.Rotation = 0;

        % Tick labels rotated 90 deg: invisible axes overlaid on the colorbar
        % (the colorbar has no tick-label rotation property)
        tickVals   = get(hh,'Ticks');
        tickLabels = compose('%g',tickVals);
        hh.TickLabels = [];

        pos = hh.Position;
        ax_overlay = axes('Position',pos, 'Color','none', 'Visible','off', ...
                          'XLim', [min(tickVals) max(tickVals)], 'YLim',[0 1], ...
                          'Box','off');
        set(ax_overlay, 'XAxisLocation','bottom', 'YAxisLocation','left');
        ax_overlay.XTick      = tickVals;
        ax_overlay.XTickLabel = tickLabels;
        xtickangle(ax_overlay, 90);
        ax_overlay.FontName = 'Times New Roman';
        ax_overlay.FontSize = 12;

        m_gshhs_f('patch', [.83 .83 .83], 'linewidth', 0.00001);
        m_grid('linestyle', ':', 'FontName', 'Times New Roman', 'FontSize', 12);

        title('(c) SWOT L3 (DAC-CalVal)');
        set(gcf, 'DefaultTextFontName', 'Times New Roman', 'DefaultTextFontSize', 12);
        set(gcf, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontSize', 12);
        ax = gca; set(ax, 'FontName', 'Times New Roman', 'FontSize', 12);

        print('-dpng', './Fig_Papers/Fig01', '-r720');
        hold off;

        clear hxc hxc1;
    end
end
