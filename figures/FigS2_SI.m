% =========================================================================
% FIGURE S2  -  HOVMOLLER OF DAC-CORRECTED SSHA (DEGRADATION)
% -------------------------------------------------------------------------
% Reuses the exact filtering + cross-track-collapse pipeline of Fig07.m
% (Hovmoller IB-only), but applied to each DAC-corrected swath field, to
% show that the northward-propagating coastally trapped waves resolved
% under IB-only lose continuity in the 2-10 day band once any DAC is applied.
%
% Output: 6 stacked Hovmoller panels (latitude vs time) sharing one colour
%   scale -> (a) IB only, (b) Trad, (c) CalVal, (d) MSAS, (e) SIROCCO,
%   (f) DAC-ERA5.  Companion to Fig. 7a.
%
% Input files: one processed .mat per correction, each holding a struct
% with fields .sla [t x ny x nx], .lat, .lon, .t (same layout as 'swot_ib'
% in swot007L3_processed_IB_hp_na.mat).
% =========================================================================
clear; close all; clc;
try, addpath("../../m_map1.4f"); catch, end

% --- correction -> (matfile, struct variable name inside the .mat) --------
% IB-only reuses the file already used by Fig07.m.
CORR = {
  'IB',       'swot007L3_processed_IB_hp_na.mat',     'swot_ib'
  'Trad',     'swot007L3_processed_Trad.mat',         'swot_dac'
  'CalVal',   'swot007L3_processed_CalVal.mat',       'swot_dac'
  'MSAS',     'swot007L3_processed_MSAS.mat',         'swot_dac'
  'SIROCCO',  'swot007L3_processed_SIROCCO.mat',      'swot_dac'
  'DAC-ERA5', 'swot007L3_processed_DACERA5.mat',      'swot_dac'
};
panel_letters = {'a','b','c','d','e','f'};

cutoff_days  = 10;          % same high-pass cutoff as Fig07.m
color_limits = [-0.3 0.3];  % shared SSHA colour scale [m]
FONT = 'Times New Roman';   FS = 16;

% Divergent colormap identical to Fig07.m
aux1 = [0 0 139; 0 0 203; 80 80 254; 141 141 254; 223 223 254; 254 254 254; ...
        254 254 24; 254 146 24; 254 75 24; 203 0 0; 139 0 0]/255;
aux2 = linspace(1,100,size(aux1,1)); aux3 = 1:100;
cmap = [interp1(aux2,aux1(:,1),aux3,'v5cubic')', ...
        interp1(aux2,aux1(:,2),aux3,'v5cubic')', ...
        interp1(aux2,aux1(:,3),aux3,'v5cubic')'];
cmap(cmap>1)=1; cmap(cmap<0)=0;

fig = figure('Units','inches','Position',[1 1 9 16],'Visible','off');
tlo = tiledlayout(numel(panel_letters),1,'TileSpacing','compact','Padding','compact');

for k = 1:size(CORR,1)
    label = CORR{k,1}; matfile = CORR{k,2}; varname = CORR{k,3};
    ax = nexttile; hold(ax,'on');

    if ~exist(matfile,'file')
        title(ax, sprintf('(%s) %s  [missing: %s]', panel_letters{k}, label, matfile), ...
              'FontName',FONT,'FontSize',FS,'Interpreter','none');
        set(ax,'FontName',FONT,'FontSize',FS); continue;
    end
    S = load(matfile, varname);
    swot_data = S.(varname);

    % --- mask pixels with >=10% NaNs (identical to Fig07.m) ---------------
    ntime = size(swot_data.sla,1);
    nanc = squeeze(sum(isnan(swot_data.sla),1));
    nanpct = (nanc/ntime)*100; mask = ones(size(nanpct)); mask(nanpct>=10)=NaN;
    for t = 1:ntime
        swot_data.sla(t,:,:) = squeeze(swot_data.sla(t,:,:)).*mask;
    end

    % --- temporal high-pass (Fig07.m methodology) -------------------------
    dt_days = nanmean(diff(swot_data.t)); fs = 1/dt_days;
    cutoff_freq = 1/cutoff_days;
    if cutoff_freq >= fs/2, cutoff_freq = 0.8*(fs/2); end
    sla_f = apply_filtering_methodology(swot_data.sla, swot_data.t, cutoff_freq, fs);

    % --- collapse cross-track -> latitude-time ----------------------------
    sla_f = nanmean(sla_f,2);
    lat   = nanmean(swot_data.lat,1); lat = lat(:);
    tvec  = swot_data.t;
    valid = find(~all(all(isnan(sla_f),3),2));
    tsub  = tvec(valid);
    treg  = floor(tsub(1)):1:ceil(tsub(end));
    [~,ti] = ismember(round(tsub),treg);
    HOV = nan(numel(lat),numel(treg));
    for ii = 1:numel(valid)
        col = squeeze(nanmean(sla_f(valid(ii),:,:),2));
        if ti(ii)>0, HOV(:,ti(ii)) = col; end
    end

    [TT,YY] = meshgrid(treg,lat);
    pcolor(ax,TT,YY,HOV); shading(ax,'interp');
    contour(ax,TT,YY,HOV,[0 0],'k:','LineWidth',0.6);
    caxis(ax,color_limits); colormap(ax,cmap);
    ylim(ax,[-53 -41.15]); xlim(ax,[treg(1) treg(end)]);
    set(ax,'Layer','top','Box','on','FontName',FONT,'FontSize',FS);
    ax.YTickLabel = abs(ax.YTick);
    ylabel(ax,['Lat (' char(176) 'S)'],'FontName',FONT,'FontSize',FS);
    title(ax,sprintf('(%s) %s',panel_letters{k},label),'FontName',FONT,'FontSize',FS);
    if k==size(CORR,1), datetick(ax,'x','dd/mm','keepticks'); xtickangle(ax,45);
    else, set(ax,'XTickLabel',{}); end
end
cb = colorbar('Position',[0.92 0.35 0.02 0.3]);
ylabel(cb,'SSHA [m]','FontName',FONT,'FontSize',FS); set(cb,'FontName',FONT,'FontSize',FS);

print(fig,'-dpng','FigS2.png','-r600','-opengl');
close(fig);
fprintf('Figure saved: FigS2.png\n');

% ===== local function (verbatim from Fig07.m) ============================
function filtered_data = apply_filtering_methodology(data, time_vector, cutoff_freq, fs)
    [T, nlon, nlat] = size(data);
    filtered_data = nan(T, nlon, nlat);
    t_days = time_vector - time_vector(1);
    for i = 1:nlon
        for j = 1:nlat
            ts = squeeze(data(:,i,j)); ts = ts(:);
            vi = ~isnan(ts);
            if sum(vi) < 20, continue; end
            tf = ts;
            if sum(vi) >= 10
                p = polyfit(t_days(vi), ts(vi), 1); tf = ts - polyval(p,t_days)';
            end
            if sum(vi) >= 15 && T >= 30
                w30 = 2*pi/30;
                A = [ones(T,1), cos(w30*t_days(:)), sin(w30*t_days(:))];
                c = A(vi,:)\tf(vi); tf = tf - A*c;
            end
            if sum(vi) >= 10
                mv = datevec(time_vector); mo = mv(:,2);
                for m = unique(mo)'
                    mi = (mo==m) & vi;
                    if sum(mi) >= 3, tf(mi) = tf(mi) - nanmean(tf(mi)); end
                end
            end
            if sum(vi) >= 30
                ti2 = tf;
                if any(~vi), ti2(~vi) = interp1(t_days(vi),tf(vi),t_days(~vi),'linear','extrap'); end
                try
                    [b,a] = butter(2, cutoff_freq/(fs/2), 'high');
                    ff = filtfilt(b,a,ti2); ff(~vi) = NaN; tf = ff;
                catch, end
            end
            filtered_data(:,i,j) = tf;
        end
    end
end
