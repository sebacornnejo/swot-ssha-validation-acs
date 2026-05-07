%% Spectral Analysis of SWOT and In-Situ Sea Surface Height Anomalies
% Description:
%   Computes the Power Spectral Density (PSD) of Sea Surface Height 
%   Anomalies (SSHA) for individual monitoring stations and regional 
%   composites over the Argentine Continental Shelf. The time series are 
%   detrended, and evaluated using a tapered periodogram approach to 
%   maximize frequency resolution over the short (~100-day) observational 
%   window of the SWOT Cal/Val phase.

clear; close all; clc;

%% Configuration and Directories
dirs.input_shared = './FES2022vsFES2022';
dirs.input_hr     = './Series_insitu_FES2022b_SIROCCO'; 
dirs.output       = './Spectrals/Combined_Analysis_FES2022b';

if ~exist(dirs.output, 'dir'); mkdir(dirs.output); end
dirs.out_indiv    = fullfile(dirs.output, 'Individual_Stations');
dirs.out_regional = fullfile(dirs.output, 'Regional_Averages');
if ~exist(dirs.out_indiv, 'dir'); mkdir(dirs.out_indiv); end
if ~exist(dirs.out_regional, 'dir'); mkdir(dirs.out_regional); end

% Quality control parameters
params.max_variance_iterations = 5;

% Spectral estimation parameters
params.spec.min_valid_fraction   = 0.70;   
params.spec.max_missing_days     = 5;      
params.spec.force_dt_days        = 1.0;    
params.spec.window_type          = 'tukey';
params.spec.tukey_alpha          = 0.5;    
params.spec.min_segment_days     = 16;     
params.spec.segment_fraction_ind = 0.90; % Evaluation window fraction (Individual)
params.spec.segment_fraction_reg = 0.90; % Evaluation window fraction (Regional)
params.spec.cutoff_period_days   = 4 * 30.4; 

% Plotting aesthetics
params.colors.insitu    = [18, 191, 128] / 255;
params.colors.swot      = [99, 49, 197] / 255;
params.colors.red_noise = [0.9, 0.2, 0.2];
params.colors.swot_base = [0.6, 0.6, 0.6]; 

% Regional composite exclusion configuration
params.exclude_stations = {'San Matias'};       

%% Station Definitions
locations = struct(...
    'name',           {'Lander', 'Carina', 'Vega-Pleyade', 'P. Deseado', 'San Matias'}, ...
    'shared_filename',{'LANDER_FES2022vFES2022.xlsx', ...
                       'CARINA_FES2022vFES2022.xlsx', ...
                       'VEGA_FES2022vFES2022.xlsx', ...
                       'PDESEADO_FES2022vFES2022.xlsx', ...
                       'SBE26_FES2022vFES2022.xlsx'}, ...
    'hr_filename',    {'Serie_Lander_SLA_FES2022b_SIROCCO_CalVal.xlsx', ...
                       'Serie_CARINA_SLA_FES2022b_SIROCCO_CalVal.xlsx', ...
                       'Serie_VEGA_SLA_FES2022b_SIROCCO_CalVal.xlsx', ...
                       'Serie_PDESEADO_SLA_FES2022b_SIROCCO_CalVal.xlsx', ...
                       'Serie_SBE26_SLA_FES2022b_SIROCCO_CalVal.xlsx'}, ...
    'lat',            {-52.14, -52.75718, -53.299, -47.75, -41.18223}, ...
    'lon',            {-66.67, -67.21960, -67.746, -65.917, -63.78029}, ...
    'ylim_power',     {1, 1, 1, 2, 3} ...
);

%% 1. Data Assimilation and Quality Control
fprintf('------------------------------------------------------------\n');
fprintf('Part 1: Data Assimilation and Quality Control\n');
fprintf('------------------------------------------------------------\n');

dataset = struct('name', {}, 'shared_dates', {}, 'insitu_sir', {}, 'swot_sir', {}, 'swot_trad', {}, 'hr_dates', {}, 'hr_sla', {});

for loc_idx = 1:length(locations)
    current_loc = locations(loc_idx);
    fprintf('Assimilating station: %s\n', current_loc.name);
    dataset(loc_idx).name = current_loc.name;
    
    % Shared quasi-daily arrays
    shared_path = fullfile(dirs.input_shared, current_loc.shared_filename);
    if exist(shared_path, 'file')
        [dates_sir, insitu_sir, swot_sir] = preprocess_sheet_data(shared_path, 'SIROCCO', params);
        [dates_trad, ~, swot_trad]        = preprocess_sheet_data(shared_path, 'Trad', params);
        
        if isempty(dates_sir) || isempty(dates_trad) || length(dates_sir) ~= length(dates_trad)
            fprintf('  -> Warning: Data integrity failure for shared records. Skipping.\n');
            continue;
        end
        dataset(loc_idx).shared_dates = dates_sir;
        dataset(loc_idx).insitu_sir   = insitu_sir;
        dataset(loc_idx).swot_sir     = swot_sir;
        dataset(loc_idx).swot_trad    = swot_trad;
    else
        fprintf('  -> Warning: Shared file not found: %s\n', shared_path);
        continue;
    end
    
    % High-resolution arrays
    hr_path = fullfile(dirs.input_hr, current_loc.hr_filename);
    if exist(hr_path, 'file')
        try
            hr_opts  = detectImportOptions(hr_path);
            hr_data  = readtable(hr_path, hr_opts);
            dataset(loc_idx).hr_dates = hr_data{:, 2}; 
            dataset(loc_idx).hr_sla   = hr_data{:, 3};
        catch
            fprintf('  -> Warning: Failed to parse high-resolution file.\n');
        end
    end
end

%% 2. Individual Station Spectral Analysis
fprintf('\n------------------------------------------------------------\n');
fprintf('Part 2: Individual Station Spectra\n');
fprintf('------------------------------------------------------------\n');

for loc_idx = 1:length(dataset)
    if isempty(dataset(loc_idx).shared_dates); continue; end
    fprintf('Evaluating spectra for station: %s\n', dataset(loc_idx).name);
    
    dnum_all   = dataset(loc_idx).shared_dates;
    insitu_sir = dataset(loc_idx).insitu_sir;
    swot_sir   = dataset(loc_idx).swot_sir;
    swot_trad  = dataset(loc_idx).swot_trad;
    
    [grid_dates_num, dt_grid] = build_time_grid(dnum_all, insitu_sir, swot_sir, swot_trad, params);
    if isempty(grid_dates_num); continue; end
    
    xg_insitu = resample_to_grid(dnum_all, insitu_sir, grid_dates_num, params, dt_grid);
    xg_swot   = resample_to_grid(dnum_all, swot_sir,   grid_dates_num, params, dt_grid);
    xg_trad   = resample_to_grid(dnum_all, swot_trad,  grid_dates_num, params, dt_grid);
    
    has_hr = ~isempty(dataset(loc_idx).hr_dates);
    
    fig = figure('Units', 'inches', 'Position', [1, 1, 8, 5], 'Visible', 'off');
    ax = axes(fig);
    title_str = sprintf('Station: %s', dataset(loc_idx).name);
    ylim_max = 10^locations(loc_idx).ylim_power;
    
    plot_power_spectrum(ax, grid_dates_num, xg_insitu, xg_swot, xg_trad, ...
        has_hr, dataset(loc_idx).hr_dates, dataset(loc_idx).hr_sla, ...
        title_str, params, true, true, locations(loc_idx), ylim_max, params.spec.segment_fraction_ind);
        
    fig_filename = fullfile(dirs.out_indiv, sprintf('Spectrum_%s.png', dataset(loc_idx).name));
    print(fig, fig_filename, '-dpng', '-r600');
    close(fig);
end

%% 3. Regional Composite Spectral Analysis
fprintf('\n------------------------------------------------------------\n');
fprintf('Part 3: Regional Composite Spectra\n');
fprintf('------------------------------------------------------------\n');

mask_all = true(1, numel(dataset));
mask_no_sanmatias = true(1, numel(dataset));

for k = 1:numel(dataset)
    if isempty(dataset(k).shared_dates)
        mask_all(k) = false;
        mask_no_sanmatias(k) = false;
    end
    if any(strcmpi(dataset(k).name, params.exclude_stations))
        mask_no_sanmatias(k) = false;
    end
end

regional_cases = struct('tag', {'ALL_STATIONS', 'NO_SAN_MATIAS'}, 'mask', {mask_all, mask_no_sanmatias});

for cc = 1:numel(regional_cases)
    case_tag  = regional_cases(cc).tag;
    case_mask = regional_cases(cc).mask;
    
    fprintf('Evaluating composite: %s (%d stations)\n', case_tag, nnz(case_mask));
    if nnz(case_mask) == 0; continue; end

    station_indices = find(case_mask);
    [grid_dates_num, dt_grid] = build_global_grid(dataset, params, station_indices);
    if isempty(grid_dates_num); continue; end
    
    num_stations = numel(station_indices);
    insitu_mat = NaN(numel(grid_dates_num), num_stations);
    swot_mat_s = NaN(numel(grid_dates_num), num_stations);
    swot_mat_t = NaN(numel(grid_dates_num), num_stations);
    
    for kk = 1:num_stations
        k = station_indices(kk);
        dnum = dataset(k).shared_dates;
        insitu_mat(:,kk) = resample_to_grid(dnum, dataset(k).insitu_sir, grid_dates_num, params, dt_grid);
        swot_mat_s(:,kk) = resample_to_grid(dnum, dataset(k).swot_sir,   grid_dates_num, params, dt_grid);
        swot_mat_t(:,kk) = resample_to_grid(dnum, dataset(k).swot_trad,  grid_dates_num, params, dt_grid);
    end
    
    insitu_avg = mean(insitu_mat, 2, 'omitnan');
    swot_avg_s = mean(swot_mat_s, 2, 'omitnan');
    swot_avg_t = mean(swot_mat_t, 2, 'omitnan');
    
    [f_i, psd_i, dof_i, dt_i, x_i_used] = compute_tapered_periodogram(grid_dates_num, insitu_avg, params, false, params.spec.segment_fraction_reg);
    [f_s, psd_s]                        = compute_tapered_periodogram(grid_dates_num, swot_avg_s, params, false, params.spec.segment_fraction_reg);
    [f_t, psd_t]                        = compute_tapered_periodogram(grid_dates_num, swot_avg_t, params, false, params.spec.segment_fraction_reg);
    
    psd_hr_avg = [];
    if ~isempty(f_i)
        psd_hr_mat = NaN(length(f_i), num_stations);
        for kk = 1:num_stations
            k = station_indices(kk);
            if ~isempty(dataset(k).hr_dates) && ~isempty(dataset(k).hr_sla)
                [f_hr, p_hr] = compute_tapered_periodogram(dataset(k).hr_dates, dataset(k).hr_sla, params, true, params.spec.segment_fraction_reg);
                if ~isempty(f_hr)
                    psd_hr_mat(:, kk) = interp1(f_hr, p_hr, f_i, 'linear', NaN);
                end
            end
        end
        psd_hr_avg = mean(psd_hr_mat, 2, 'omitnan');
    end
    
    fig = figure('Units', 'inches', 'Position', [1, 1, 8, 5], 'Visible', 'off');
    ax = axes(fig);
    title_str = sprintf('Regional Composite: %s', strrep(case_tag, '_', ' '));
    ylim_max = 10^2; 
    
    plot_regional_spectrum(ax, f_i, psd_i, f_s, psd_s, f_t, psd_t, psd_hr_avg, dof_i, dt_i, x_i_used, title_str, params, true, ylim_max, cc);
    
    fig_filename = fullfile(dirs.out_regional, sprintf('Spectrum_Regional_%s.png', case_tag));
    print(fig, fig_filename, '-dpng', '-r600');
    close(fig);
end

fprintf('\n================== SPECTRAL ANALYSIS COMPLETED ==================\n');

%% Local Functions

function [dates_num, insitu_clean, swot_clean] = preprocess_sheet_data(file_path, target_sheet, params)
    dates_num = []; insitu_clean = []; swot_clean = [];
    try
        opts = detectImportOptions(file_path, 'Sheet', target_sheet, 'ReadVariableNames', true);
        varNames = opts.VariableNames;
        for k = 2:length(varNames); opts = setvartype(opts, varNames{k}, 'double'); end
        opts.MissingRule = 'fill'; opts.ImportErrorRule = 'fill';
        dataTable = readtable(file_path, opts);
    catch
        return;
    end
    
    if isdatetime(dataTable{:,1})
        dates_num = datenum(dataTable{:,1});
    elseif iscellstr(dataTable{:,1}) || isstring(dataTable{:,1})
        dates_num = datenum(dataTable{:,1});
    else
        dates_num = dataTable{:,1};
    end

    cols_to_check = [2, 3]; 
    dataTable_cleaned = apply_variance_qc(dataTable, cols_to_check, params);
    
    insitu_col = 2; swot_col = 3;
    
    mask_common_initial = ~isnan(dataTable_cleaned{:, insitu_col}) & ~isnan(dataTable_cleaned{:, swot_col});
    idx_first_common = find(mask_common_initial, 1, 'first');
    if ~isempty(idx_first_common)
        dataTable_cleaned{idx_first_common, insitu_col} = NaN;
        dataTable_cleaned{idx_first_common, swot_col} = NaN;
    end
    
    mask_common = ~isnan(dataTable_cleaned{:, insitu_col}) & ~isnan(dataTable_cleaned{:, swot_col});
    if sum(mask_common) > 0
        mean_in = mean(dataTable_cleaned{mask_common, insitu_col});
        mean_sw = mean(dataTable_cleaned{mask_common, swot_col});
        dataTable_cleaned{:, insitu_col} = dataTable_cleaned{:, insitu_col} - mean_in;
        dataTable_cleaned{:, swot_col}   = dataTable_cleaned{:, swot_col} - mean_sw;
    end
    
    insitu_clean = dataTable_cleaned{:, insitu_col};
    swot_clean   = dataTable_cleaned{:, swot_col};
end

function dataTable_qc = apply_variance_qc(dataTable, cols_to_check, params)
    dataTable_qc = dataTable;
    last_dev = inf;
    for iter = 1:params.max_variance_iterations + 1
        max_dev = -1; row_idx = -1; col_idx = -1;
        for k = 1:length(cols_to_check)
            c_col = cols_to_check(k);
            c_data = dataTable_qc{:, c_col};
            v_mask = ~isnan(c_data);
            if sum(v_mask) < 3; continue; end
            d_clean = c_data(v_mask);
            devs = abs(c_data - mean(d_clean));
            [c_max_dev, c_idx] = max(devs);
            if c_max_dev > (2 * std(d_clean)) && c_max_dev > max_dev
                max_dev = c_max_dev; row_idx = c_idx; col_idx = c_col;
            end
        end
        if row_idx == -1; break; end
        if iter > 1 && (max_dev <= (2/3) * last_dev); break; end
        if iter > params.max_variance_iterations; break; end
        last_dev = max_dev;
        dataTable_qc{row_idx, col_idx} = NaN;
    end
end

function [grid_dates, dt_grid] = build_time_grid(dnum_all, x1, x2, x3, params)
    dnum_all = dnum_all(:);
    mask_any = isfinite(dnum_all) & (~isnan(x1(:)) | ~isnan(x2(:)) | ~isnan(x3(:)));
    d = unique(sort(dnum_all(mask_any)));
    if numel(d) < 2
        grid_dates = []; dt_grid = []; return;
    end
    if ~isempty(params.spec.force_dt_days) && isfinite(params.spec.force_dt_days)
        dt_grid = params.spec.force_dt_days;
    else
        dt_grid = median(diff(d));
        if ~isfinite(dt_grid) || dt_grid <= 0; dt_grid = 1.0; end
        if abs(dt_grid - 1.0) < 0.20; dt_grid = 1.0; end
    end
    grid_dates = (d(1):dt_grid:d(end)).';
end

function [grid_dates, dt_grid] = build_global_grid(dataset, params, station_indices)
    all_dates = [];
    for kk = 1:numel(station_indices)
        k = station_indices(kk);
        dnum = dataset(k).shared_dates;
        valid_mask = ~isnan(dataset(k).insitu_sir) | ~isnan(dataset(k).swot_sir) | ~isnan(dataset(k).swot_trad);
        all_dates = [all_dates; dnum(valid_mask)];
    end
    all_dates = unique(sort(all_dates));
    if numel(all_dates) < 2
        grid_dates = []; dt_grid = []; return;
    end
    if ~isempty(params.spec.force_dt_days) && isfinite(params.spec.force_dt_days)
        dt_grid = params.spec.force_dt_days;
    else
        dt_grid = median(diff(all_dates));
        if ~isfinite(dt_grid) || dt_grid <= 0; dt_grid = 1.0; end
        if abs(dt_grid - 1.0) < 0.20; dt_grid = 1.0; end
    end
    grid_dates = (all_dates(1):dt_grid:all_dates(end)).';
end

function xg = resample_to_grid(dates_num, x, grid_dates_num, params, dt_grid)
    xg = NaN(numel(grid_dates_num), 1);
    valid = isfinite(dates_num(:)) & ~isnan(x(:));
    if nnz(valid) < max(2, params.spec.min_segment_days); return; end
    d = dates_num(valid); x_v = x(valid);
    [d, order] = sort(d); x_v = x_v(order);
    [dates_unique, ~, ic] = unique(d);
    if numel(dates_unique) < numel(d)
        x_v = accumarray(ic, x_v, [], @mean);
        d = dates_unique;
    end
    dt_local = median(diff(d));
    if ~isfinite(dt_local) || dt_local <= 0; dt_local = dt_grid; end
    if abs(dt_local - 1.0) < 0.20; dt_local = 1.0; end
    missing_samples = max(0, round(diff(d)./dt_local) - 1);
    if any(missing_samples > params.spec.max_missing_days); return; end
    span_mask = grid_dates_num >= d(1) & grid_dates_num <= d(end);
    n_span = nnz(span_mask);
    if n_span <= 0 || (numel(x_v) / n_span) < params.spec.min_valid_fraction; return; end
    xg(span_mask) = interp1(d, x_v, grid_dates_num(span_mask), 'linear');
end

function [freq, psd, dof, dt, x_used] = compute_tapered_periodogram(dates_num, data, params, is_highres, fraction)
    freq = []; psd = []; dof = []; dt = []; x_used = [];
    valid = isfinite(dates_num(:)) & isfinite(data(:));
    if nnz(valid) < max(8, params.spec.min_segment_days); return; end
    
    d = dates_num(valid); x = data(valid);
    [d, idx] = sort(d); x = x(idx);
    
    dt_local = median(diff(d));
    if ~isfinite(dt_local) || dt_local <= 0; return; end
    
    if ~is_highres && ~isempty(params.spec.force_dt_days) && isfinite(params.spec.force_dt_days)
        dt = params.spec.force_dt_days;
    elseif ~is_highres && abs(dt_local - 1.0) < 0.20
        dt = 1.0;
    else
        dt = dt_local;
    end
    
    [tgrid, xgrid] = enforce_uniform_sampling(d, x, dt, params, is_highres);
    if isempty(tgrid) || isempty(xgrid); return; end
    
    fs = 1 / dt; n = numel(xgrid);
    xgrid = detrend(xgrid, 'linear');
    x_used = xgrid;
    
    if is_highres
        duration_days = tgrid(end) - tgrid(1);
        eval_window_days = max(params.spec.min_segment_days, fraction * duration_days);
        window_length = round(eval_window_days / dt);
        window_length = min(window_length, n);
    else
        min_window_len = round(params.spec.min_segment_days / dt);
        window_length = max(min_window_len, min(n, round(fraction * n)));
    end
    
    if window_length < 4 || window_length > n; return; end
    
    noverlap = round(0.50 * window_length);
    noverlap = min(noverlap, window_length - 1);
    
    if strcmpi(params.spec.window_type, 'tukey')
        win = tukeywin(window_length, params.spec.tukey_alpha);
    else
        win = hann(window_length);
    end
    
    nfft = 2^nextpow2(window_length);
    
    [psd_calc, f_calc] = pwelch(xgrid, win, noverlap, nfft, fs, 'psd');
    dof = compute_effective_dof(win, window_length, noverlap, n);
    
    maskf = f_calc > 0 & isfinite(psd_calc);
    freq  = f_calc(maskf);
    psd   = psd_calc(maskf);
end

function dof = compute_effective_dof(win, window_length, noverlap, n)
    step = window_length - noverlap;
    if step <= 0 || n < window_length; dof = 2; return; end
    K = max(1, 1 + floor((n - window_length) / step));
    if K == 1 || noverlap == 0; dof = max(2, 2*K); return; end
    E = sum(win.^2);
    if ~isfinite(E) || E <= 0; dof = max(2, 2*K); return; end
    sum_term = 0;
    for m = 1:(K-1)
        lag = m * step; if lag >= window_length; break; end
        c = sum(win(1:end-lag) .* win(1+lag:end)) / E;
        sum_term = sum_term + (1 - m/K) * (c^2);
    end
    Keff = max(1, min(K, K / (1 + 2*sum_term)));
    dof = max(2, 2*Keff);
end

function [tgrid, xgrid] = enforce_uniform_sampling(d, x, dt, params, is_highres)
    tgrid = []; xgrid = [];
    if numel(d) < 2; return; end
    max_missing_pts = ceil(params.spec.max_missing_days / dt);
    missing = max(0, round(diff(d)./dt) - 1);
    if any(missing > max_missing_pts); return; end
    tgrid = (d(1):dt:d(end)).';
    if numel(tgrid) < 2; return; end
    if ~is_highres && (numel(x) / numel(tgrid)) < params.spec.min_valid_fraction; return; end
    xgrid = interp1(d, x, tgrid, 'linear');
    if any(isnan(xgrid))
        xgrid = fillmissing(xgrid, 'linear', 'EndValues', 'nearest');
    end
end

function h_group = plot_red_noise_confidence(ax, freq, x_used, dof, dt, color)
    h_group = hggroup('Parent', ax);
    if isempty(dof) || dof < 2 || isempty(freq) || isempty(x_used); return; end
    x = x_used(:) - mean(x_used, 'omitnan');
    if numel(x) < 3; return; end
    
    den = sum(x(1:end-1).^2);
    if den <= 0 || ~isfinite(den); return; end
    r1 = max(-0.99, min(0.99, sum(x(1:end-1) .* x(2:end)) / den));
    
    sigma2 = var(x, 1);
    red_noise = (2*dt) * sigma2 * (1 - r1^2) ./ (1 - 2*r1*cos(2*pi*freq*dt) + r1^2);
    
    qL = chi2inv_local(0.025, dof);
    qU = chi2inv_local(0.975, dof);
    
    loglog(ax, freq, red_noise * (qL / dof), ':', 'Color', color, 'LineWidth', 0.8, 'Parent', h_group);
    loglog(ax, freq, red_noise * (qU / dof), ':', 'Color', color, 'LineWidth', 0.8, 'Parent', h_group);
    loglog(ax, freq, red_noise, '--', 'Color', color, 'LineWidth', 1.0, 'Parent', h_group);
    set(get(get(h_group,'Annotation'),'LegendInformation'), 'IconDisplayStyle', 'on');
    uistack(h_group, 'bottom');
end

function x = chi2inv_local(p, v)
    if exist('chi2inv', 'file') == 2
        x = chi2inv(p, v);
    else
        x = 2 * gammaincinv(p, v/2);
    end
end

function format_frequency_axis(ax)
    periods_days  = [30.4, 21, 14, 7, 5, 4, 3];
    period_labels = {'1m', '3w', '2w', '1w', '5d', '4d', '3d'};
    tick_freqs = 1 ./ periods_days;
    
    current_xlim = xlim(ax);
    visible = (tick_freqs >= current_xlim(1)) & (tick_freqs <= current_xlim(2));
    
    if any(visible)
        set(ax, 'XTick', tick_freqs(visible), 'XTickLabel', period_labels(visible));
        xlabel(ax, 'Period', 'FontWeight', 'bold', 'FontName', 'Times New Roman', 'FontSize', 20);
    else
        xlabel(ax, 'Frequency (cpd)', 'FontWeight', 'bold', 'FontName', 'Times New Roman', 'FontSize', 20);
    end
    set(ax, 'XMinorTick', 'on', 'FontName', 'Times New Roman', 'FontSize', 20);
end

function plot_power_spectrum(ax, dates_grid_num, insitu_shared, swot_shared, trad_shared, ...
    has_highres, dates_hr, sla_hr, title_str, params, show_legend, show_station_info, loc_info, ylim_max, seg_fraction)
    
    hold(ax, 'on');
    fmin = 1 / params.spec.cutoff_period_days;
    
    [fi, Pii, dof_i, dt_i, x_used] = compute_tapered_periodogram(dates_grid_num, insitu_shared, params, false, seg_fraction);
    if isempty(fi) || isempty(Pii)
        text(ax, 0.5, 0.5, 'Insufficient valid data for spectrum', 'Units','normalized', ...
            'HorizontalAlignment','center', 'FontName','Times New Roman', 'FontSize', 20);
        set(ax, 'XScale','log', 'YScale','log'); grid(ax,'on'); box(ax,'on');
        title(ax, title_str, 'FontWeight','bold', 'FontName','Times New Roman', 'FontSize', 20);
        return;
    end
    
    mask_i = fi >= fmin; fi_k = fi(mask_i); Pii_k = Pii(mask_i);
    [ft, Ptt] = compute_tapered_periodogram(dates_grid_num, trad_shared, params, false, seg_fraction);
    mask_t = ft >= fmin; ft_k = ft(mask_t); Ptt_k = Ptt(mask_t);
    [fs, Pss] = compute_tapered_periodogram(dates_grid_num, swot_shared, params, false, seg_fraction);
    mask_s = fs >= fmin; fs_k = fs(mask_s); Pss_k = Pss(mask_s);

    f_hr_k = []; P_hr_k = [];
    if has_highres
        [f_hr, P_hr] = compute_tapered_periodogram(dates_hr, sla_hr, params, true, seg_fraction);
        if ~isempty(f_hr)
            mask_hr = f_hr >= fmin; f_hr_k = f_hr(mask_hr); P_hr_k = P_hr(mask_hr);
        end
    end

    h_trad = gobjects(0); h_hr = gobjects(0); h_insitu = gobjects(0); h_swot = gobjects(0);
    h_red_group = gobjects(0);
    
    if ~isempty(ft_k); h_trad = loglog(ax, ft_k, Ptt_k, ':', 'Color', params.colors.swot_base, 'LineWidth', 2.0); end
    if ~isempty(f_hr_k); h_hr = loglog(ax, f_hr_k, P_hr_k, ':', 'Color', params.colors.insitu, 'LineWidth', 2.0); end
    if ~isempty(fi_k); h_insitu = loglog(ax, fi_k, Pii_k, '-', 'Color', params.colors.insitu, 'LineWidth', 2.0); end
    if ~isempty(fs_k); h_swot = loglog(ax, fs_k, Pss_k, '-', 'Color', params.colors.swot, 'LineWidth', 2.0); end
    if ~isempty(fi_k)
        h_red_group = plot_red_noise_confidence(ax, fi_k, x_used, dof_i, dt_i, params.colors.red_noise);
    end

    set(ax, 'XScale', 'log', 'YScale', 'log', 'YAxisLocation', 'right');
    grid(ax, 'on'); box(ax, 'on');
    ylabel(ax, 'PSD [m^2/cpd]', 'FontWeight', 'bold', 'FontName', 'Times New Roman', 'FontSize', 20);
    title(ax, title_str, 'FontWeight', 'bold', 'FontName', 'Times New Roman', 'FontSize', 20);
    xlim(ax, [1/30.4, 1/3]);
    ylim(ax, [1e-5, ylim_max]);
    set(ax, 'FontName', 'Times New Roman', 'FontSize', 20);
    format_frequency_axis(ax);
    
    if show_station_info && nargin >= 13
        info_text = sprintf('%.2f°S, %.2f°W', abs(loc_info.lat), abs(loc_info.lon));
        text(ax, 0.98, 0.95, info_text, 'Units', 'normalized', ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
            'BackgroundColor', 'w', 'EdgeColor', 'k', 'FontWeight', 'bold', ...
            'FontName', 'Times New Roman', 'FontSize', 16);
    end
    
    if show_legend
        hlist = []; llist = {};
        if isgraphics(h_trad),      hlist(end+1) = h_trad;      llist{end+1} = 'SWOT (DAC Trad.)'; end 
        if isgraphics(h_hr),        hlist(end+1) = h_hr;        llist{end+1} = 'In-situ Original Sampling'; end
        if isgraphics(h_insitu),    hlist(end+1) = h_insitu;    llist{end+1} = 'In-situ (DAC SIROCCO)'; end 
        if isgraphics(h_swot),      hlist(end+1) = h_swot;      llist{end+1} = 'SWOT (DAC SIROCCO)'; end 
        if isgraphics(h_red_group), hlist(end+1) = h_red_group; llist{end+1} = 'Red Noise (95% CI)'; end 
        if ~isempty(hlist)
            lgd = legend(ax, hlist, llist, 'Location', 'southwest');
            set(lgd, 'FontName', 'Times New Roman', 'FontSize', 14); 
        end
    end
end

function plot_regional_spectrum(ax, f_i, psd_i, f_s, psd_s, f_t, psd_t, psd_hr_avg, dof_i, dt_i, x_i_used, title_str, params, show_legend, ylim_max, case_index)
    hold(ax, 'on');
    if isempty(f_i) || isempty(psd_i)
        text(ax, 0.5, 0.5, 'Insufficient valid data for composite', 'Units','normalized', ...
            'HorizontalAlignment','center', 'FontName','Times New Roman', 'FontSize', 20);
        set(ax, 'XScale','log', 'YScale','log'); grid(ax,'on'); box(ax,'on');
        title(ax, title_str, 'FontWeight','bold', 'FontName','Times New Roman', 'FontSize', 20);
        return;
    end
    
    h_red_group = plot_red_noise_confidence(ax, f_i, x_i_used, dof_i, dt_i, params.colors.red_noise);
    
    h_trad = gobjects(0); h_hr = gobjects(0); h_swot = gobjects(0);
    if ~isempty(f_t) && ~isempty(psd_t)
        h_trad = loglog(ax, f_t, psd_t, ':', 'Color', params.colors.swot_base, 'LineWidth', 2.0);
    end
    if ~isempty(psd_hr_avg)
        h_hr = loglog(ax, f_i, psd_hr_avg, ':', 'Color', params.colors.insitu, 'LineWidth', 2.0);
    end
    h_insitu = loglog(ax, f_i, psd_i, '-', 'Color', params.colors.insitu, 'LineWidth', 2.0);
    if ~isempty(f_s) && ~isempty(psd_s)
        h_swot = loglog(ax, f_s, psd_s, '-', 'Color', params.colors.swot, 'LineWidth', 2.0);
    end
    
    if show_legend
        hlist = []; llist = {};
        if isgraphics(h_trad),      hlist(end+1) = h_trad;      llist{end+1} = 'SWOT Avg (DAC Trad.)'; end
        if isgraphics(h_hr),        hlist(end+1) = h_hr;        llist{end+1} = 'In-situ Avg High-Res'; end
        if isgraphics(h_insitu),    hlist(end+1) = h_insitu;    llist{end+1} = 'In-situ Avg (DAC SIROCCO)'; end
        if isgraphics(h_swot),      hlist(end+1) = h_swot;      llist{end+1} = 'SWOT Avg (DAC SIROCCO)'; end
        if isgraphics(h_red_group), hlist(end+1) = h_red_group; llist{end+1} = 'Red Noise (95% CI)'; end
        if ~isempty(hlist)
            lgd = legend(ax, hlist, llist, 'Location', 'southwest');
            set(lgd, 'FontName', 'Times New Roman', 'FontSize', 16);
        end
    end
    
    set(ax, 'XScale', 'log', 'YScale', 'log', 'YAxisLocation', 'right');
    grid(ax, 'on'); box(ax, 'on');
    ylabel(ax, 'PSD [m^2/cpd]', 'FontWeight', 'bold', 'FontName', 'Times New Roman', 'FontSize', 20);
    title(ax, title_str, 'FontWeight', 'bold', 'FontName', 'Times New Roman', 'FontSize', 20);
    xlim(ax, [1/30.4, 1/3]);
    
    if case_index == 1
        ylim(ax, [10^(-5), ylim_max]);      
    else
        ylim(ax, [10^(-5.5), ylim_max]);    
    end
    
    set(ax, 'FontName', 'Times New Roman', 'FontSize', 20);
    format_frequency_axis(ax);
end