%% SWOT L3 Spatiotemporal Analysis: CEOF Propagation & Hovmöller Diagrams
% Description:
% This script performs a complete analysis of SWOT L3 Sea Surface Height 
% Anomalies (SSHA) for Track 007 during the Cal/Val phase. 
% It compares two Dynamic Atmospheric Corrections (DAC): 
%   1) Default DAC (Traditional)
%   2) SIROCCO DAC (Regional high-resolution)
%
% THE WORKFLOW INCLUDES:
%   1. High-pass filtering (< 10 days) to isolate synoptic-scale waves,
%      following adapted methodology from Dinapoli et al. (2025).
%   2. Complex Empirical Orthogonal Function (CEOF) analysis to identify
%      propagating features in the Southwestern Atlantic.
%   3. Construction of Hovmöller diagrams (Latitude vs. Time) for Mode 1
%      reconstructions and their differences.

clear; close all; clc;

% --- Path Setup (Optional) ---
% addpath("../../m_map1.4f"); 

%% 1. DATA LOADING AND INITIALIZATION
% =========================================================================
fprintf('Step 1: Loading processed SWOT L3 data...\n');
input_file = 'swot007L3_processed_reflorent_na.mat';

if ~exist(input_file, 'file')
    error('Input file %s not found in path.', input_file);
end

load(input_file); % Expected variables: swot_na (Default), swot_re (SIROCCO)

% Filtering Parameters
cutoff_days = 10; % High-pass filter threshold
dt_days = nanmean(diff(swot_na.t)); 
fs = 1/dt_days; 
cutoff_freq = 1/cutoff_days;

datasets = {'na', 're'};
dataset_labels = {'Default Correction', 'SIROCCO Correction'};

%% 2. SIGNAL PROCESSING (DINAPOLI ET AL. 2025 ADAPTATION)
% =========================================================================
% We remove low-frequency components (seasonal cycles and interannual drift)
% to focus on the sub-10-day variability.
fprintf('Step 2: Applying adapted high-pass filtering methodology...\n');

filtered_structs = struct();

for d = 1:numel(datasets)
    ds_name = datasets{d};
    fprintf('  -> Processing %s dataset...\n', ds_name);
    
    current_ds = eval(['swot_' ds_name]);
    
    % Apply multi-step filtering (Detrending -> Mean removal -> 10-day High-pass)
    filtered_sla = apply_filtering_methodology(current_ds.sla, current_ds.t, cutoff_freq, fs);
    
    % Store in a temporary structure
    filtered_structs.(ds_name).sla = filtered_sla;
    filtered_structs.(ds_name).t   = current_ds.t;
    filtered_structs.(ds_name).lat = current_ds.lat;
    filtered_structs.(ds_name).lon = current_ds.lon;
end

%% 3. COMPLEX EOF (CEOF) ANALYSIS
% =========================================================================
% CEOF identifies propagating modes by using the Hilbert Transform to 
% construct an analytic signal: Z(t,x) = SSHA(t,x) + i*H(SSHA(t,x)).
fprintf('Step 3: Performing CEOF analysis and reconstructing Mode 1...\n');

reconstructed_fields = struct();
variance_explained = nan(numel(datasets), 1);

for d = 1:numel(datasets)
    ds_name = datasets{d};
    curr = filtered_structs.(ds_name);
    [nt, nlon, nlat] = size(curr.sla);
    
    % Flatten 3D spatial dimensions to 2D for EOF computation
    % Identify valid ocean points (non-NaN)
    mean_map = squeeze(nanmean(curr.sla, 1));
    ocean_idx = find(~isnan(mean_map));
    
    if isempty(ocean_idx)
        error('No valid ocean points found for %s.', ds_name);
    end
    
    % Reshape valid points: [Time x Space]
    M_real = nan(nt, length(ocean_idx));
    for i = 1:length(ocean_idx)
        [r, c] = ind2sub([nlon, nlat], ocean_idx(i));
        M_real(:, i) = curr.sla(:, r, c);
    end
    
    % Fill sporadic NaNs with zero (neutral for EOF) and remove mean
    M_real(isnan(M_real)) = 0;
    M_real_anom = M_real - mean(M_real, 1);
    
    % Construct analytic signal via Hilbert Transform
    M_complex = hilbert(M_real_anom);
    
    % Call the EOF function (Keeping the user-specified function name)
    n_modes = 1; 
    [L, var, CEOFs, PCs, ~] = EOF_(M_complex, n_modes);
    
    variance_explained(d) = var(1);
    fprintf('  -> %s Mode 1 Variance: %.2f%%\n', dataset_labels{d}, var(1));
    
    % Reconstruct the 3D SSHA field using only Mode 1
    % Field = Real part of (PC1 * Spatial_Mode1')
    recon_2d = real(PCs(:, 1) * conj(CEOFs(:, 1))');
    
    % Map back to 3D grid
    sla_recon = nan(nt, nlon, nlat);
    for i = 1:length(ocean_idx)
        [r, c] = ind2sub([nlon, nlat], ocean_idx(i));
        sla_recon(:, r, c) = recon_2d(:, i);
    end
    
    reconstructed_fields.(ds_name) = sla_recon;
end

%% 4. HOVMÖLLER DIAGRAM CONSTRUCTION
% =========================================================================
fprintf('Step 4: Building Hovmöller matrices (Latitude vs. Time)...\n');

HOV_matrices = cell(2, 1);
time_grids = cell(2, 1);
lat_vectors = cell(2, 1);

for d = 1:numel(datasets)
    ds_name = datasets{d};
    curr_filt = filtered_structs.(ds_name);
    curr_recon = reconstructed_fields.(ds_name);
    
    % Define the Zonal Average slice (Mean across Longitude)
    lat_slice = squeeze(nanmean(curr_filt.lat, 1));
    
    % Create regular daily time grid
    t_start = floor(min(curr_filt.t));
    t_end = ceil(max(curr_filt.t));
    reg_time = t_start:1:t_end;
    
    % Compute zonal mean of the reconstructed Mode 1
    zonal_mean_recon = squeeze(nanmean(curr_recon, 2)); % [Time x Latitude]
    
    % Map irregular SWOT sampling to daily grid
    hov_daily = nan(length(lat_slice), length(reg_time));
    [~, time_loc] = ismember(round(curr_filt.t), reg_time);
    
    for ii = 1:length(time_loc)
        if time_loc(ii) > 0
            hov_daily(:, time_loc(ii)) = zonal_mean_recon(ii, :);
        end
    end
    
    HOV_matrices{d} = hov_daily;
    time_grids{d}   = reg_time;
    lat_vectors{d}  = lat_slice;
end

% --- Difference Matrix (SIROCCO - Default) ---
% We align both datasets to the same Latitude/Time grid (using 'na' as reference)
common_lat = lat_vectors{1};
t_common_start = min(time_grids{1}(1), time_grids{2}(1));
t_common_end   = max(time_grids{1}(end), time_grids{2}(end));
common_time = t_common_start:1:t_common_end;

% Align 'na' (A)
A = nan(length(common_lat), length(common_time));
[~, idxA] = ismember(time_grids{1}, common_time);
A(:, idxA) = HOV_matrices{1};

% Align 're' (B)
B_raw = HOV_matrices{2};
% Interpolate B to A's latitude if they differ
if ~isequal(lat_vectors{2}, common_lat)
    B_raw = interp_lat_matrix(B_raw, lat_vectors{2}, common_lat);
end
B = nan(length(common_lat), length(common_time));
[~, idxB] = ismember(time_grids{2}, common_time);
B(:, idxB) = B_raw;

% Final Difference (Physical Difference)
D = B - A;

%% 5. VISUALIZATION
% =========================================================================
fprintf('Step 5: Generating final figures...\n');

% Custom Colormap Construction
cmap_colors = [0,0,139; 0,0,203; 80,80,254; 141,141,254; 223,223,254; 254,254,254; ...
               254,254,24; 254,146,24; 254,75,24; 203,0,0; 139,0,0] / 255;
cmap = interp1(linspace(1,100,size(cmap_colors,1)), cmap_colors, 1:100, 'v5cubic');
cmap = max(0, min(1, cmap));

fig = figure('Units','inches','Position',[1 1 10 18],'Color','w');

% Plot Panels
panel_titles = {sprintf('%s CEOF Mode 1 (%.1f%%)', dataset_labels{1}, variance_explained(1)), ...
                sprintf('%s CEOF Mode 1 (%.1f%%)', dataset_labels{2}, variance_explained(2)), ...
                'Difference (SIROCCO - Default)'};

for p = 1:3
    subplot(3,1,p);
    if p == 1, data_plot = A; t_plot = common_time; clim_val = [-0.015 0.015]; end
    if p == 2, data_plot = B; t_plot = common_time; clim_val = [-0.015 0.015]; end
    if p == 3, data_plot = D; t_plot = common_time; clim_val = [-0.01 0.01]; end
    
    [TT, YY] = meshgrid(t_plot, common_lat);
    pcolor(TT, YY, data_plot); shading interp;
    hold on;
    contour(TT, YY, data_plot, [0 0], 'k:', 'LineWidth', 0.8);
    
    title(panel_titles{p}, 'FontSize', 18, 'FontName', 'Times New Roman');
    ylabel(['Latitude (' char(176) 'S)'], 'FontSize', 16);
    caxis(clim_val);
    colormap(gca, cmap);
    cb = colorbar; ylabel(cb, 'SLA [m]');
    
    ylim([-53 -41.15]);
    datetick('x', 'dd/mm', 'keepticks');
    set(gca, 'Layer', 'top', 'TickDir', 'in', 'FontSize', 14, 'FontName', 'Times New Roman');
    
    % Format Y-axis to show absolute values
    yticks = get(gca, 'YTick');
    set(gca, 'YTickLabel', abs(yticks));
end

fprintf('Process complete. Hovmöller diagrams generated.\n');

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function filtered_data = apply_filtering_methodology(data, time_vector, cutoff_freq, fs)
    [T, nlon, nlat] = size(data);
    filtered_data = nan(T, nlon, nlat);
    t_days = time_vector - time_vector(1);
    
    for i = 1:nlon
        for j = 1:nlat
            ts = squeeze(data(:, i, j));
            ts = ts(:);
            valid_idx = ~isnan(ts);
            
            if sum(valid_idx) < 30, continue; end
            
            % 1. Linear Detrend
            p = polyfit(t_days(valid_idx), ts(valid_idx), 1);
            ts_filt = ts - polyval(p, t_days);
            
            % 2. Adapted Seasonal/Low-frequency removal (Period > 30 days)
            A = [ones(length(t_days), 1), cos(2*pi/30 * t_days(:)), sin(2*pi/30 * t_days(:))];
            coef = A(valid_idx, :) \ ts_filt(valid_idx);
            ts_filt = ts_filt - A * coef;
            
            % 3. Remove Monthly Means
            [~, months] = datevec(time_vector);
            unique_months = unique(months);
            for m = unique_months'
                m_idx = (months == m) & valid_idx;
                if sum(m_idx) >= 3
                    ts_filt(m_idx) = ts_filt(m_idx) - nanmean(ts_filt(m_idx));
                end
            end
            
            % 4. 2nd Order Butterworth High-pass Filter
            ts_interp = ts_filt;
            ts_interp(~valid_idx) = interp1(t_days(valid_idx), ts_filt(valid_idx), t_days(~valid_idx), 'linear', 'extrap');
            [b, a] = butter(2, cutoff_freq/(fs/2), 'high');
            ts_final = filtfilt(b, a, ts_interp);
            ts_final(~valid_idx) = NaN;
            
            filtered_data(:, i, j) = ts_final;
        end
    end
end

function data_out = interp_lat_matrix(data_in, lat_in, lat_out)
    lat_in = lat_in(:); lat_out = lat_out(:);
    ok = ~isnan(lat_in);
    [lat_u, ia] = unique(lat_in(ok));
    data_u = data_in(ia, :);
    data_out = nan(length(lat_out), size(data_in, 2));
    for k = 1:size(data_in, 2)
        data_out(:, k) = interp1(lat_u, data_u(:, k), lat_out, 'linear', NaN);
    end
end