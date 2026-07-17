%% =======================================================================
%  FIGURE S1 - Cross-spectral assessment of the DACs against the ERA5-MSLP
%  inverse-barometer (IB) reference, per station.
%  D = IB series - DAC series. Layout 2x3: top row = IN SITU, bottom row =
%  SWOT. Columns: 1) gamma^2(D, IB series) in the 2-10 d band
%  2) coherence-weighted phase lag of D (h)  3) FVR (%).
%  Inputs : series_FES2022vsFES2022_*/<STA>_*.xlsx (SSHA pairs per DAC)
%  Outputs: ./Figures_MSLP/FigS1.png
%           ./Stats_MSLP/Master_CrossSpectral_DAC_vs_IBMSLP.xlsx
%  =======================================================================
clear; close all; clc;

input_dir  = '.';
out_stats  = './Stats_MSLP';
out_figs   = './Figures_MSLP';

if ~exist(out_stats,'dir'); mkdir(out_stats); end
if ~exist(out_figs ,'dir'); mkdir(out_figs );  end

stations = {'CARINA','LANDER','PDESEADO','SBE26','VEGA'};
dacs       = {'Trad','CalVal','MSAS','SIROCCO','ERA5'};

map_est    = containers.Map({'SBE26','PDESEADO','VEGA','CARINA','LANDER'}, ...
             {'G. San Matias','P. Deseado','Vega-Pleyade','Carina','Lander'});

config.max_homogeneity_iterations = 5;

% North -> South order for the figure rows
order_NS = {'SBE26','PDESEADO','LANDER','CARINA','VEGA'};

nS = numel(stations); nD = numel(dacs);
G2_SW   = nan(nS,nD);  G2_IN   = nan(nS,nD);
FVR_SW  = nan(nS,nD);  FVR_IN  = nan(nS,nD);
RDCTRL  = nan(nS,nD);

% --- coherence-weighted lag of D per side (hours) ---
LAGW_SW = nan(nS,nD);  LAGW_IN = nan(nS,nD);

fprintf('=== Cross-spectral DAC assessment vs IB-MSLP: coherence + weighted lag ===\n');

for i = 1:nS
    est = stations{i};
    fprintf('\n>>> Station %s\n', est);

    % --- IB reference series (ERA5 MSLP) ---
    f_ib = fullfile(input_dir,'series_FES2022vsFES2022_DACMSLPERA5', ...
                    sprintf('%s_FES2022vsFES2022_DACMSLPERA5.xlsx', est));
    [t_ib, in_ib, sw_ib] = load_and_clean(f_ib, '', false, config);
    if isempty(t_ib); warning('No IB-MSLP series for %s. Skipped.', est); continue; end

    for j = 1:nD
        dac = dacs{j};
        switch dac
            case {'Trad','CalVal','MSAS'}
                f_d = fullfile(input_dir,'Paper_2025a','FES2022vsFES2022', ...
                               sprintf('%s_FES2022vFES2022.xlsx', est));
                [t_d, in_d, sw_d] = load_and_clean(f_d, dac, true, config);
            case 'SIROCCO'
                f_d = fullfile(input_dir,'series_FES2022vsFES2022', ...
                               sprintf('%s_FES2022vsFES2022.xlsx', est));
                [t_d, in_d, sw_d] = load_and_clean(f_d, '', false, config);
            case 'ERA5'
                f_d = fullfile(input_dir,'series_FES2022vsFES2022_DACERA5', ...
                               sprintf('%s_FES2022vsFES2022_DACERA5.xlsx', est));
                [t_d, in_d, sw_d] = load_and_clean(f_d, '', false, config);
        end
        if isempty(t_d); continue; end

        key_ib = round(t_ib * 86400);
        key_d  = round(t_d  * 86400);
        [~, ia, ibx] = intersect(key_ib, key_d);

        if numel(ia) < 30 && numel(t_ib) == numel(t_d)
            ia = (1:numel(t_ib))'; ibx = ia;
        elseif numel(ia) < 30
            warning('%s-%s: incompatible time bases. Skipped.', est, dac); continue;
        end

        tA = t_ib(ia);
        M  = [in_ib(ia), sw_ib(ia), in_d(ibx), sw_d(ibx)];
        ok = all(~isnan(M), 2);
        if sum(ok) < 30; continue; end

        i1 = find(ok,1,'first'); i2 = find(ok,1,'last');
        tA = tA(i1:i2); M = M(i1:i2,:); ok = ok(i1:i2);

        dt = median(diff(tA(ok)));
        tu = (tA(1):dt:tA(end))';
        Fu = nan(numel(tu),4);
        for c = 1:4; Fu(:,c) = interp1(tA(ok), M(ok,c), tu, 'pchip'); end
        Fu = detrend(Fu);

        in_ibU = Fu(:,1); sw_ibU = Fu(:,2);
        in_dU  = Fu(:,3); sw_dU  = Fu(:,4);

        D_in = in_ibU - in_dU;
        D_sw = sw_ibU - sw_dU;

        RDCTRL(i,j) = corr(D_sw, D_in);

        fs  = 1/dt;
        wl  = round(30/dt); if wl > numel(tu); wl = numel(tu); end
        nov = round(wl/2);  win = hanning(wl);

        % --- Coherence D vs IB series in the 2-10 d band ---
        [c_sw, fq] = mscohere(D_sw, sw_ibU, win, nov, wl, fs);
        [c_in, ~ ] = mscohere(D_in, in_ibU, win, nov, wl, fs);
        pd   = 1 ./ fq;
        band = (pd >= 2) & (pd <= 10);

        G2_SW(i,j) = mean(c_sw(band));
        G2_IN(i,j) = mean(c_in(band));

        % --- Cross phase D vs IB series -> weighted lag (h) ---
        [pxy_sw, ~] = cpsd(D_sw, sw_ibU, win, nov, wl, fs);
        [pxy_in, ~] = cpsd(D_in, in_ibU, win, nov, wl, fs);

        lag_sw = (angle(pxy_sw)*(180/pi)) ./ (360*fq) * 24;   % hours
        lag_in = (angle(pxy_in)*(180/pi)) ./ (360*fq) * 24;

        w_sw = c_sw(band) / sum(c_sw(band));
        w_in = c_in(band) / sum(c_in(band));

        LAGW_SW(i,j) = sum(abs(lag_sw(band)) .* w_sw);
        LAGW_IN(i,j) = sum(abs(lag_in(band)) .* w_in);

        % --- FVR ---
        FVR_SW(i,j) = 100 * (1 - var(sw_dU) / var(sw_ibU));
        FVR_IN(i,j) = 100 * (1 - var(in_dU) / var(in_ibU));
    end
end

%% ============ FIGURE (2x3): IN-SITU row / SWOT row =====================
% Reorder rows North -> South
[~, ord] = ismember(order_NS, stations);
row_lbl  = cellfun(@(e) map_est(e), order_NS, 'UniformOutput', false);

G2i  = G2_IN(ord,:);  G2s  = G2_SW(ord,:);
Li   = LAGW_IN(ord,:); Ls   = LAGW_SW(ord,:);
Fi   = FVR_IN(ord,:);  Fs   = FVR_SW(ord,:);

FONT = 'Times New Roman'; FS = 20;

lagMax = max([Li(:); Ls(:)],[],'omitnan'); lagCax = [0 ceil(lagMax)];
fvrLo  = min([Fi(:); Fs(:)],[],'omitnan'); fvrHi = max([Fi(:); Fs(:)],[],'omitnan');
fvrCax = [floor(fvrLo/10)*10, ceil(fvrHi/10)*10];

scr = get(groot,'ScreenSize');
fig = figure('Visible','off','Color','w','Units','pixels','Position',[1 1 scr(3) scr(4)]);
tlo = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

% --- Row 1: IN SITU ---
ax = nexttile;
heat20(ax, G2i, row_lbl, dacs, '(a) \gamma^2(D, In-Situ_{IB-MSLP})  2-10 d', [0 1],   seq_cmap(),  '%.2f', FONT, FS, true);

ax = nexttile;
heat20(ax, Li,  row_lbl, dacs, '(b) Weighted phase lag, In-Situ_{IB-MSLP} (h)', lagCax, seq_cmap(), '%.1f', FONT, FS, false);

ax = nexttile;
heat20(ax, Fi,  row_lbl, dacs, '(c) FVR In-Situ_{IB-MSLP} (%)', fvrCax, seq_cmap(), '%+.0f', FONT, FS, false);

% --- Row 2: SWOT ---
ax = nexttile;
heat20(ax, G2s, row_lbl, dacs, '(d) \gamma^2(D, SWOT_{IB-MSLP})  2-10 d', [0 1],   seq_cmap(),  '%.2f', FONT, FS, true);

ax = nexttile;
heat20(ax, Ls,  row_lbl, dacs, '(e) Weighted phase lag, SWOT_{IB-MSLP} (h)', lagCax, seq_cmap(), '%.1f', FONT, FS, false);

ax = nexttile;
heat20(ax, Fs,  row_lbl, dacs, '(f) FVR SWOT_{IB-MSLP} (%)', fvrCax, seq_cmap(), '%+.0f', FONT, FS, false);

print(fig, fullfile(out_figs,'FigS1.png'), '-dpng','-r600');
fprintf('\n[FIGURE SAVED]: %s\n', fullfile(out_figs,'FigS1.png'));
close(fig);

%% ===================== BACKUP EXCEL ====================================
xls = fullfile(out_stats,'Master_CrossSpectral_DAC_vs_IBMSLP.xlsx');
if exist(xls,'file'); delete(xls); end

rl = cellfun(@(e) map_est(e), stations, 'UniformOutput', false);
mats = {G2_SW,'g2_D_vs_SWOT'; G2_IN,'g2_D_vs_InSitu'; ...
        LAGW_SW,'LagPond_SWOT_h'; LAGW_IN,'LagPond_InSitu_h'; ...
        FVR_SW,'FVR_SWOT_pct'; FVR_IN,'FVR_InSitu_pct'; ...
        RDCTRL,'Control_corr_Dsw_Din'};

for m = 1:size(mats,1)
    Tm = array2table(mats{m,1},'VariableNames',dacs,'RowNames',rl);
    writetable(Tm, xls, 'Sheet', mats{m,2}, 'WriteRowNames', true);
end
fprintf('[EXCEL SAVED]: %s\n', xls);

%% ========================= LOCAL FUNCTIONS =============================
function heat20(ax, M, rows, cols, ttl, cl, cmap, fmt, FONT, FS, showY)
    imagesc(ax, M, 'AlphaData', ~isnan(M)); hold(ax,'on');
    set(ax,'Color',[0.94 0.94 0.94]);
    colormap(ax, cmap); caxis(ax, cl); cb = colorbar(ax);
    cb.FontName = FONT; cb.FontSize = FS;

    if showY, yl = rows; else, yl = repmat({''},1,numel(rows)); end

    set(ax,'YDir','reverse','XTick',1:numel(cols),'XTickLabel',cols, ...
        'XTickLabelRotation',0, ...
        'YTick',1:numel(rows),'YTickLabel',yl,'TickLabelInterpreter','none', ...
        'FontName',FONT,'FontSize',FS,'Layer','top');

    xlim(ax,[0.5 numel(cols)+0.5]); ylim(ax,[0.5 numel(rows)+0.5]);

    for r = 1:size(M,1)
        for c = 1:size(M,2)
            if isnan(M(r,c)); continue; end
            v  = (M(r,c) - cl(1)) / (cl(2) - cl(1));
            tc = 'k'; if v > 0.72 || v < 0.22; tc = 'w'; end
            text(ax, c, r, sprintf(fmt, M(r,c)), 'HorizontalAlignment','center', ...
                 'FontName',FONT,'FontSize',FS,'FontWeight','bold','Color',tc);
        end
    end
    title(ax, ttl, 'FontWeight','normal','FontName',FONT,'FontSize',FS);
end

function [t, insitu, swot] = load_and_clean(fpath, sheet, multis, config)
    t = []; insitu = []; swot = [];
    if ~exist(fpath,'file'); warning('Not found: %s', fpath); return; end
    try
        if multis
            opts = detectImportOptions(fpath,'Sheet',sheet,'ReadVariableNames',true);
        else
            opts = detectImportOptions(fpath,'ReadVariableNames',true);
        end
        for k = 2:numel(opts.VariableNames)
            opts = setvartype(opts, opts.VariableNames{k}, 'double');
        end
        T = readtable(fpath, opts);
    catch
        warning('Could not read: %s', fpath); return;
    end

    [T, ~] = iterative_outlier_removal(T, [2 3], config);

    m  = ~isnan(T{:,2}) & ~isnan(T{:,3});
    i1 = find(m, 1, 'first');
    if ~isempty(i1); T{i1,2} = NaN; T{i1,3} = NaN; end

    m = ~isnan(T{:,2}) & ~isnan(T{:,3});
    if sum(m) > 2
        T{:,2} = T{:,2} - mean(T{m,2});
        T{:,3} = T{:,3} - mean(T{m,3});
    end

    traw = T{:,1};
    if isdatetime(traw); t = datenum(traw); else; t = traw; end
    insitu = T{:,2}; swot = T{:,3};
end

function cm = seq_cmap()
    a = [68 1 84; 59 82 139; 33 145 140; 94 201 98; 253 231 37]/255;
    x = linspace(0,1,size(a,1)); xi = linspace(0,1,256);
    cm = [interp1(x,a(:,1),xi)' interp1(x,a(:,2),xi)' interp1(x,a(:,3),xi)'];
end

function [dataTable, outlier_info_table] = iterative_outlier_removal(dataTable, cols_to_check, config)
    dataTable_cleaned  = dataTable;
    outlier_candidates = table('Size',[0,3],'VariableTypes',{'double','double','double'}, ...
                               'VariableNames',{'RowIndex','ColIndex','DeviationValue'});
    last_removed_deviation = inf; halt_due_to_homogeneity = false;
    outlier_info_table = table();
    var_names = dataTable.Properties.VariableNames;

    for iter = 1:config.max_homogeneity_iterations + 1
        max_dev_this_pass = -1; candidate_row_idx = -1; candidate_col_idx = -1;
        for k = 1:length(cols_to_check)
            current_col  = cols_to_check(k);
            current_data = dataTable_cleaned{:, current_col};
            valid_mask   = ~isnan(current_data);
            if sum(valid_mask) < 3; continue; end
            data_clean   = current_data(valid_mask);
            data_mean    = mean(data_clean); data_std = std(data_clean);
            deviations   = abs(current_data - data_mean);
            [current_max_dev, current_idx] = max(deviations);
            if current_max_dev > (2*data_std) && current_max_dev > max_dev_this_pass
                max_dev_this_pass = current_max_dev;
                candidate_row_idx = current_idx; candidate_col_idx = current_col;
            end
        end
        if candidate_row_idx == -1; break; end
        if iter > 1 && (max_dev_this_pass <= (2/3)*last_removed_deviation); break; end
        if iter > config.max_homogeneity_iterations; halt_due_to_homogeneity = true; break; end

        outlier_candidates = [outlier_candidates; {candidate_row_idx, candidate_col_idx, max_dev_this_pass}];
        last_removed_deviation = max_dev_this_pass;
        dataTable_cleaned{candidate_row_idx, candidate_col_idx} = NaN;
    end

    if ~halt_due_to_homogeneity && ~isempty(outlier_candidates)
        outlier_info_table = table('Size',[0,3],'VariableTypes',{'string','string','double'}, ...
                                   'VariableNames',{'Date','Variable','Deviation'});
        for k = 1:height(outlier_candidates)
            row_to_remove = outlier_candidates.RowIndex(k);
            col_to_remove = outlier_candidates.ColIndex(k);
            dev_value     = outlier_candidates.DeviationValue(k);
            outlier_date  = dataTable{row_to_remove, 1};
            variable_name = var_names{col_to_remove};
            if isnumeric(outlier_date)
                date_str = datestr(outlier_date, 'dd-mm-yyyy');
            else
                date_str = char(outlier_date);
            end
            outlier_info_table = [outlier_info_table; {date_str, variable_name, dev_value}];
            dataTable{row_to_remove, col_to_remove} = NaN;
        end
    end
end
