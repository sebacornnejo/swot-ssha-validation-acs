%% =======================================================================
%  INVERTED-BAROMETER SKILL BY STATION  (SWOT L2 = REFERENCE / BASE)
%  -----------------------------------------------------------------------
%  Reads 'Barometros_Invertidos_Estaciones_SWOT_CalVal.xlsx' and condenses
%  the error metrics of every pressure product against SWOT L2 into ONE
%  space-efficient figure: four stacked GROUPED-BAR panels sharing the same
%  station axis (grouped bars avoid overlapping markers):
%     (a) Correlation r
%     (b) RMSE                       (total error)
%     (c) Centred RMSD               (random / pattern error)
%     (d) sigma_product / sigma_SWOT (variability ratio)
%  + a metrics table printed and exported to CSV/XLSX.
%
%  PRODUCTS: ERA5 MSLP | CFSv2 | ECMWF (IFS HRES) | ERA5 sp
%  (SWOT L2 = reference, not plotted as a series.)
%  The reader auto-detects columns by header keywords (PRODUCTS / REF_KEYS).
%  =======================================================================
clear; close all; clc;
%% --------------------------- USER CONFIG -------------------------------
XLSX_FILE   = './Barometros_Invertidos_SWOT/Barometros_Invertidos_Estaciones_SWOT_CalVal.xlsx';
OUT_DIR     = './Barometros_Invertidos_SWOT/IB_Intercomparison_Output';
% Station order along the x-axis (left -> right).
STATION_ORDER = {'San_Matias','Puerto_Deseado','Lander','Carina','Vega'};
PRODUCTS = {
    'ERA5'   , {'era5mslp','mslp'}             % -> IB_ERA5_MSLP_m (CDS)
%    'GFS'         , {'gfs'}
    'CFSv2'       , {'cfsv2','cfs_v2','cfs2','cfs'}
    'ECMWF'       , {'ecmwf'}                       % -> IB_ECMWF_m (single IFS HRES)
    'ERA5 (SP)'     , {'era5'}                      % -> IB_ERA5_m (first 'era5' column)
};
REF_KEYS  = {'invbarcor','inv_bar_cor','swotl2','swot'};
USE_CM = true;     % RMSE / CRMSD in cm (true) or metres (false)
FONT   = 'Times New Roman';
FS     = 20;       % base font size
if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end
if ~exist(XLSX_FILE,'file')
    error('"%s" not found. Place it next to the script or fix XLSX_FILE.', XLSX_FILE);
end
%% --------------------------- PALETTE -----------------------------------
PAL = [196 196 196;   % ERA5 MSLP     (grey)
       126 222 214;   % GFS         (teal)
       195 147 255;   % CFSv2       (violet)
       250 223  90;   % ECMWF       (yellow)
       253 154  95]/255;          % ERA5 SP   (orange)
COLOR_PROD = PAL(1:size(PRODUCTS,1), :);
prodLabels = PRODUCTS(:,1)';
nProd      = numel(prodLabels);
%% --------------------------- READ DATA ---------------------------------
sheets = string(sheetnames(XLSX_FILE));
stationData = struct('name',{},'ref',{},'prod',{});
if numel(sheets) > 1
    for k = 1:numel(sheets)
        T = readtable(XLSX_FILE,'Sheet',char(sheets(k)),'VariableNamingRule','preserve');
        sd = extract_station(T, PRODUCTS, REF_KEYS, char(sheets(k)));
        if ~isempty(sd.ref), stationData(end+1) = sd; end
    end
else
    T = readtable(XLSX_FILE,'Sheet',char(sheets(1)),'VariableNamingRule','preserve');
    vn = T.Properties.VariableNames;
    iStat = findcol(vn, {'station','estacion','estación','site','nombre','name'});
    if iStat == 0
        error(['Single sheet with no station column. Rename the station ', ...
               'column (e.g. "Station") or use one sheet per station.']);
    end
    stCol = string(T{:, iStat}); uSt = unique(stCol,'stable');
    for k = 1:numel(uSt)
        sd = extract_station(T(stCol==uSt(k),:), PRODUCTS, REF_KEYS, char(uSt(k)));
        if ~isempty(sd.ref), stationData(end+1) = sd; end
    end
end
if isempty(stationData), error('No station could be extracted. Check names/keywords.'); end
names = {stationData.name}; ordIdx = [];
for k = 1:numel(STATION_ORDER)
    j = find(strcmpi(names, STATION_ORDER{k}),1);
    if ~isempty(j), ordIdx(end+1) = j; end
end
ordIdx      = [ordIdx, setdiff(1:numel(stationData), ordIdx, 'stable')];
stationData = stationData(ordIdx);
nStat       = numel(stationData);
stNames     = cellfun(@(s) strrep(s,'_',' '), {stationData.name}, 'UniformOutput', false);
%% --------------------------- METRICS -----------------------------------
R     = nan(nStat,nProd);  STDR  = nan(nStat,nProd);
BIAS  = nan(nStat,nProd);  CRMSD = nan(nStat,nProd);  RMSE = nan(nStat,nProd);
NPTS  = nan(nStat,1);
for i = 1:nStat
    ref = stationData(i).ref(:); P = stationData(i).prod;
    mask = isfinite(ref);
    for p = 1:nProd, if ~isempty(P{p}), mask = mask & isfinite(P{p}(:)); end, end
    NPTS(i) = nnz(mask);
    o = ref(mask); mo = mean(o); so = std(o,1);
    for p = 1:nProd
        if isempty(P{p}), continue; end
        x = P{p}(:); x = x(mask);
        cc = corrcoef(x,o);
        R(i,p)     = cc(1,2);
        STDR(i,p)  = std(x,1)/so;
        BIAS(i,p)  = mean(x) - mo;
        CRMSD(i,p) = sqrt(mean(((x-mean(x))-(o-mo)).^2));
        RMSE(i,p)  = sqrt(mean((x-o).^2));
    end
end
%% --------------------------- METRICS TABLE -----------------------------
rows = {}; vals = [];
for i = 1:nStat
    for p = 1:nProd
        rows(end+1,:) = {stationData(i).name, prodLabels{p}};
        vals(end+1,:) = [R(i,p), STDR(i,p), BIAS(i,p), CRMSD(i,p), RMSE(i,p), NPTS(i)];
    end
end
Tm = table(rows(:,1),rows(:,2),vals(:,1),vals(:,2),vals(:,3),vals(:,4),vals(:,5),vals(:,6), ...
    'VariableNames',{'Station','Product','r','STD_ratio','Bias_m','CRMSD_m','RMSE_m','N'});
disp('========================  METRICS (vs SWOT L2)  ========================');
disp(Tm);
writetable(Tm, fullfile(OUT_DIR,'IB_metrics_vs_SWOT.csv'));
try, writetable(Tm, fullfile(OUT_DIR,'IB_metrics_vs_SWOT.xlsx')); catch, end
%% --------------------------- UNITS -------------------------------------
if USE_CM, sc = 100; uStr = 'cm'; else, sc = 1; uStr = 'm'; end
CRMSD_d = CRMSD*sc;  RMSE_d = RMSE*sc;
%% =======================================================================
%  SINGLE FIGURE : STACKED GROUPED-BAR PANELS, SHARED STATION AXIS
%  =======================================================================
fig = figure('Visible','off','Color','w','Units','pixels','Position',[40 40 1460 1250]);
% panel geometry (4 panels)
margL = 0.115; margR = 0.035; topY = 0.900; botY = 0.070; gap = 0.022;
panelH = (topY - botY - 3*gap)/4;
panelW = 1 - margL - margR;
yb = topY - panelH - (0:3)*(panelH+gap);    % bottoms of panels (a)..(d)
xx = 1:nStat;  xlims = [0.5 nStat+0.5];
% { matrix , ylabel , baseline-from-zero? }
panels = { R,        '\itr',                 false
           RMSE_d,   sprintf('RMSE (%s)',uStr),         true
           CRMSD_d,  sprintf('RMSD (%s)',uStr), true
           STDR,     '\sigma_{IB_{prod}} / \sigma_{IB_{SWOT}}',   false };
letters = 'abcd';
hBar = [];
for k = 1:4
    M        = panels{k,1};
    ylab     = panels{k,2};
    fromZero = panels{k,3};
    ax = axes('Parent',fig,'Position',[margL yb(k) panelW panelH]); hold(ax,'on');
    b = bar(ax, xx, M, 'grouped', 'EdgeColor','k', 'LineWidth', 0.5, 'BarWidth', 0.9);
    for p = 1:nProd, b(p).FaceColor = COLOR_PROD(p,:); end
    if k == 1, hBar = b; end
    % y-limits
    v = M(:); v = v(isfinite(v));
    if fromZero
        hi = max(v); set(ax,'YLim',[0, hi*1.12]);
    else
        lo = min(v); hi = max(v); rng = hi-lo; if rng==0, rng = 1; end
        set(ax,'YLim',[lo-0.10*rng, hi+0.12*rng]);
    end
    set(ax, 'XLim', xlims, 'XTick', 1:nStat, 'FontName', FONT, 'FontSize', FS, ...
        'Box','on', 'Layer','top', 'TickLength',[0.006 0.006]);
    ax.YGrid = 'off'; ax.GridAlpha = 0.12;
    ylabel(ax, ylab, 'FontName', FONT, 'FontSize', FS);
    text(ax, 0.010, 0.90, sprintf('(%s)', letters(k)), 'Units','normalized', ...
        'FontName', FONT, 'FontSize', FS, 'FontWeight','bold');
    if k < 4
        set(ax, 'XTickLabel', {});
    else
        set(ax, 'XTickLabel', stNames); ax.XTickLabelRotation = 0;
    end
    if k == 1
        lg = legend(ax, hBar, prodLabels, 'Orientation','horizontal', ...
            'Location','northoutside', 'Box','on');
        set(lg, 'FontName', FONT, 'FontSize', FS);
    end
end
base = fullfile(OUT_DIR,'Fig03');
exportgraphics(fig, [base '.png'], 'Resolution', 600);
close(fig);
fprintf('\nFigure + metrics saved in: %s\n', OUT_DIR);
%% =======================================================================
%  LOCAL FUNCTIONS
%  =======================================================================
function sd = extract_station(T, PRODUCTS, REF_KEYS, name)
    sd = struct('name', name, 'ref', [], 'prod', {cell(size(PRODUCTS,1),1)});
    vn = T.Properties.VariableNames;
    iref = findcol(vn, REF_KEYS);
    if iref == 0, return; end
    sd.ref = getcol(T, iref);
    for p = 1:size(PRODUCTS,1)
        ip = findcol(vn, PRODUCTS{p,2});
        if ip > 0, sd.prod{p} = getcol(T, ip); end
    end
end
function idx = findcol(varnames, keys)
    idx = 0;
    vn = lower(regexprep(varnames, '[^a-zA-Z0-9]', ''));
    for k = 1:numel(keys)
        key = lower(regexprep(keys{k}, '[^a-zA-Z0-9]', ''));
        if isempty(key), continue; end
        for i = 1:numel(vn)
            if contains(vn{i}, key), idx = i; return; end
        end
    end
end
function v = getcol(T, idx)
    c = T{:, idx};
    if isnumeric(c) || islogical(c)
        v = double(c);
    elseif iscell(c) || isstring(c) || ischar(c)
        v = str2double(string(c));
    else
        v = nan(height(T),1);
    end
    v = v(:);
end
