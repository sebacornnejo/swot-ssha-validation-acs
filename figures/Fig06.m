% =========================================================================
% TEMPORAL FILTERING + SWATH COLLAPSE + SLA HOVMOLLER (IB)
% + phase velocity and characteristic period (Lomb-Scargle, 2-10 day band)
% Manual axes and manual date ticks (no datetick).
% =========================================================================
clear all; close all; clc;
try
    addpath("../../m_map1.4f");
catch
end
%% 1. LOAD AND FILTER DATA (IB)
load('swot007L3_processed_IB_hp_na.mat'); % contains swot_na, swot_ib, swot_ib_hp
swot_data = swot_ib;
clear swot_na swot_ib swot_ib_hp;
fprintf('swot_data.sla: %s | %d timesteps\n', mat2str(size(swot_data.sla)), length(swot_data.t));
T = length(swot_data.t);
dt_days = nanmean(diff(swot_data.t));
fs = 1/dt_days;
cutoff_days = 10;
cutoff_freq = 1/cutoff_days;
nyquist_freq = fs/2;
if cutoff_freq >= nyquist_freq
    warning('Cutoff frequency too high for this time series');
    cutoff_freq = nyquist_freq * 0.8;
    fprintf('Cutoff frequency adjusted to: %.4f [1/day]\n', cutoff_freq);
end
%% MASKING (based on the percentage of NaNs in the IB field)
ntime = size(swot_data.sla, 1);
for t = 1:ntime
    if t == 1
        nan_count = squeeze(isnan(swot_data.sla(t,:,:)));   % [ny x nx]
    else
        nan_count = nan_count + squeeze(isnan(swot_data.sla(t,:,:)));
    end
end
nan_pct = (nan_count / ntime) * 100;
nan_pct(nan_pct>=10) = NaN; nan_pct = nan_pct.*0+1;
for t = 1:ntime
    swot_data.sla(t,:,:) = squeeze(swot_data.sla(t,:,:)).*nan_pct;
end
swot_filtered = swot_data;
swot_filtered.sla = apply_filtering_methodology(swot_data.sla, swot_data.t, cutoff_freq, fs);
%% 2. SWATH COLLAPSE (CROSS-TRACK AVERAGE)
swot_filtered.sla = nanmean(swot_filtered.sla, 2);
swot_filtered.lat = nanmean(swot_filtered.lat, 1);
swot_filtered.lon = nanmean(swot_filtered.lon, 1);
output_filename_filt = 'swot007L3_processed_IB_filtered_track_line.mat';
fprintf('Saving filtered track to: %s\n', output_filename_filt);
save(output_filename_filt, 'swot_filtered', 'cutoff_freq', 'fs', 'dt_days', '-v7.3');
%% 3. HOVMOLLER + PHASE VELOCITY
% Point pairs (time, latitude) digitized along the crests
x_coords = [7.389850541271989e+05, 7.389860284167795e+05,7.389941641566266e+05, 7.389948418674698e+05, 7.389972816265060e+05, ...
            7.389986370481928e+05, 7.390014834337350e+05, 7.390027033132530e+05, ...
            7.390105647590362e+05, 7.390115135542168e+05, 7.390116490963856e+05, ...
            7.390127334337350e+05, 7.390151731927710e+05, 7.390163930722892e+05, ...
            7.390172063253012e+05, 7.390184262048192e+05, 7.390220858433734e+05, ...
            7.390231701807228e+05, 7.390329292168674e+05, 7.390349623493976e+05, ...
            7.390365888554216e+05, 7.390375376506024e+05, 7.390392996987952e+05, ...
            7.390411972891566e+05, 7.390455346385542e+05, 7.390467545180722e+05, ...
            7.390547515060240e+05, 7.390557003012048e+05, 7.390563780120482e+05, ...
            7.390575978915662e+05, 7.390588177710844e+05, 7.390605798192772e+05, ...
            7.390609864457832e+05, 7.390624774096386e+05, 7.390630195783132e+05, ...
            7.390639683734940e+05, 7.390650527108434e+05, 7.390662725903614e+05];
y_coords = [-52.9126, -41.2374, -52.8407, -41.1819, -52.8407, -41.3730, -52.8407, -41.3093, ...
            -52.8407, -41.3730, -52.7770, -41.3730, -52.8407, -41.3730, ...
            -52.8407, -41.3093, -52.9044, -41.1819, -52.8407, -41.3093, ...
            -52.8407, -41.5004, -52.8407, -41.3093, -52.9044, -41.3730, ...
            -52.8407, -41.3093, -52.7770, -41.3093, -52.8407, -41.3093, ...
            -52.7133, -41.3093, -52.8407, -41.3730, -52.8407, -41.3093];
x_pairs = reshape(x_coords, 2, [])';
y_pairs = reshape(y_coords, 2, [])';
R_earth = 6371000.0;
delta_y = y_pairs(:, 2) - y_pairs(:, 1);
delta_d = delta_y * (pi / 180.0) * R_earth;
delta_x = x_pairs(:, 2) - x_pairs(:, 1);
delta_t = delta_x * 86400.0;
v_phase = delta_d ./ delta_t;
mean_v  = mean(v_phase);
std_v   = std(v_phase);
mean_x  = mean(x_pairs, 2);
% --- Style and colormap ---
color_limits_sla = [-0.3 0.3];
aux1 = [0 0 139; 0 0 203; 80 80 254; 141 141 254; 223 223 254; 254 254 254; 254 254 24; 254 146 24; 254 75 24; 203 0 0; 139 0 0] / 255;
aux2 = linspace(1,100,size(aux1,1)); aux3 = 1:100;
R = interp1(aux2,aux1(:,1),aux3,'v5cubic'); G = interp1(aux2,aux1(:,2),aux3,'v5cubic'); B = interp1(aux2,aux1(:,3),aux3,'v5cubic');
cmap = [R' G' B']; cmap(cmap>1)=1; cmap(cmap<0)=0;
clear aux1 aux2 aux3 R G B
sla_filt = swot_filtered.sla;
lat = swot_filtered.lat;
time_vec = swot_filtered.t;
valid_idx = find(~all(all(isnan(sla_filt),3),2));
lat_slice = squeeze(nanmean(lat,1));
time_subset = time_vec(valid_idx);
regular_time_vec = floor(time_subset(1)):1:ceil(time_subset(end));
[~, time_indices] = ismember(round(time_subset), regular_time_vec);
hov_irreg = nan(length(lat_slice), length(time_subset));
for ii = 1:length(valid_idx)
    t = valid_idx(ii);
    hov_irreg(:,ii) = squeeze(nanmean(sla_filt(t,:,:),2));
end
HOV_disp = nan(length(lat_slice), length(regular_time_vec));
for ii = 1:length(time_indices)
    if time_indices(ii) > 0
        HOV_disp(:, time_indices(ii)) = hov_irreg(:, ii);
    end
end
common_lat = lat_slice(:);
% =========================================================================
% 3b. Characteristic period within the 2-10 day band (Lomb-Scargle, filtered)
% =========================================================================
band_lo_days = 2;  band_hi_days = 10;
f_band = linspace(1/band_hi_days, 1/band_lo_days, 400);
fprintf('\nAnalysis band: %.1f - %.1f days\n', band_lo_days, band_hi_days);
[P_lat_filt, Pxx_mean_filt] = band_periods_by_latitude(HOV_disp, regular_time_vec, f_band);
valid_P  = P_lat_filt(~isnan(P_lat_filt));
T_mean   = mean(valid_P);
T_median = median(valid_P);
T_std    = std(valid_P);
T_q1     = manual_percentile(valid_P, 25);
T_q3     = manual_percentile(valid_P, 75);
T_iqr    = T_q3 - T_q1;
T_from_fmean = 1 / mean(1 ./ valid_P);
[~, idx_avg]  = max(Pxx_mean_filt);  T_avgspec = 1/f_band(idx_avg);
fprintf('  Latitudes analysed                    : %d\n', numel(valid_P));
fprintf('  Dominant period per latitude (mean)   : %.2f days\n', T_mean);
fprintf('  Dominant period per latitude (median) : %.2f days\n', T_median);
fprintf('  Dominant period per latitude (std)    : %.2f days\n', T_std);
fprintf('  IQR: %.2f days  [Q1=%.2f, Q3=%.2f]\n', T_iqr, T_q1, T_q3);
fprintf('  Period of the mean frequency (1/mean(f)): %.2f days\n', T_from_fmean);
fprintf('  Peak of the latitude-averaged spectrum: %.2f days\n', T_avgspec);
fprintf('  Regional characteristic period (median): %.2f days\n', T_median);
Pmat_filt = per_latitude_peak_periods(HOV_disp, regular_time_vec, f_band, 5);
print_peak_rank_table(Pmat_filt, 'Median per peak rank (2-10 d band, filtered, per latitude)');
T_char_display = T_median;
% =========================================================================
% 3c. Spectral analysis on unfiltered data (console only)
% =========================================================================
hov_unf_irreg = nan(length(common_lat), length(time_subset));
for ii = 1:length(valid_idx)
    tt = valid_idx(ii);
    hov_unf_irreg(:,ii) = squeeze(nanmean(swot_data.sla(tt,:,:),2));
end
HOV_unfilt = nan(length(common_lat), length(regular_time_vec));
for ii = 1:length(time_indices)
    if time_indices(ii) > 0
        HOV_unfilt(:, time_indices(ii)) = hov_unf_irreg(:, ii);
    end
end
T_record = regular_time_vec(end) - regular_time_vec(1);
f_full   = linspace(1/T_record, 0.5, 800);
char_rel = 0.15;
[Pmat_unf_full, frac_band_unf, frac_char_unf] = per_latitude_spectral_diagnostics( ...
        HOV_unfilt, regular_time_vec, f_full, 5, ...
        1/band_hi_days, 1/band_lo_days, T_median, char_rel);
Pmat_unf_band = per_latitude_peak_periods(HOV_unfilt, regular_time_vec, f_band, 5);
print_peak_rank_table(Pmat_unf_full, 'Full spectrum (unfiltered), median per peak rank');
print_peak_rank_table(Pmat_unf_band, '2-10 d band (unfiltered), median per peak rank');
% =========================================================================
% 3d. Spectral energy partition of the unfiltered signal (console only)
% =========================================================================
fb = frac_band_unf(~isnan(frac_band_unf));
fc = frac_char_unf(~isnan(frac_char_unf));
fprintf('  Energy in the 2-10 d band / total resolvable energy (variance):\n');
fprintf('    median = %.1f %%   (IQR %.1f-%.1f %%, n latitudes = %d)\n', ...
        100*median(fb), 100*manual_percentile(fb,25), 100*manual_percentile(fb,75), numel(fb));
fprintf('  Energy at the characteristic period (%.1f d, +-%.0f%% window) / total energy:\n', ...
        T_median, 100*char_rel);
fprintf('    median = %.1f %%   (IQR %.1f-%.1f %%, n latitudes = %d)\n', ...
        100*median(fc), 100*manual_percentile(fc,25), 100*manual_percentile(fc,75), numel(fc));

% Manual date ticks: datetick(...,'keepticks') in R2026a regenerates the
% labels incorrectly; ticks every 10 days with datestr labels (plain text).
% Both panels share the same ticks; labels only on the bottom panel.
date_ticks = regular_time_vec(1):10:regular_time_vec(end);
if date_ticks(end) < regular_time_vec(end)
    date_ticks = [date_ticks, regular_time_vec(end)];
end
date_labels = cellstr(datestr(date_ticks, 'dd/mm'));

% =========================================================================
% FIGURE -- MANUAL AXES (no tiledlayout)
% =========================================================================
FONT = 'Times New Roman';
FS   = 16;
LX  = 0.105;  PW  = 0.775;  PH  = 0.355;
Y2  = 0.115;  Y1  = 0.575;
CBX = LX + PW + 0.015;  CBW = 0.020;

fig = figure('Units','inches','Position',[1 1 9 6.4],'Visible','off');

%% PANEL (a): IB Hovmoller with phase-velocity lines
ax1 = axes('Position',[LX Y1 PW PH]); hold(ax1,'on');
[T1, Y1m] = meshgrid(regular_time_vec, common_lat);
pcolor(ax1, T1, Y1m, HOV_disp); shading(ax1,'interp');
contour(ax1, T1, Y1m, HOV_disp, [0 0], 'k:', 'LineWidth', 0.6);
for i = 1:size(x_pairs, 1)
    plot(ax1, [x_pairs(i,1), x_pairs(i,2)], [y_pairs(i,1), y_pairs(i,2)], ...
        '-o', 'LineWidth', 2.0, 'MarkerSize', 4, 'MarkerFaceColor',[1 1 1],'Color',[1 1 1]);
    plot(ax1, [x_pairs(i,1), x_pairs(i,2)], [y_pairs(i,1), y_pairs(i,2)], ...
         '-o', 'LineWidth', 1.2, 'MarkerSize', 3, 'MarkerFaceColor',[26, 224, 217]/255,'Color',[26, 224, 217]/255);
end
x_txt = regular_time_vec(1) + 0.015*(regular_time_vec(end) - regular_time_vec(1));
txt_str = sprintf('T \\approx %.1fd', T_char_display);
text(ax1, x_txt, -41.45, txt_str, ...
    'FontName',FONT,'FontSize',FS-2, ...
    'HorizontalAlignment','left','VerticalAlignment','top', ...
    'BackgroundColor','w','EdgeColor','k','LineWidth',0.5,'Margin',1);
hold(ax1,'off');
title(ax1, '(a)', 'FontName',FONT,'FontSize',FS);
ylabel(ax1, ['Latitude (' char(176) 'S)'], 'FontName',FONT,'FontSize',FS);
xlim(ax1, [regular_time_vec(1) regular_time_vec(end)]);
ylim(ax1, [-53 -41.15]);
caxis(ax1, color_limits_sla); colormap(ax1, cmap);
set_axes_style(ax1, FONT, FS);
set(ax1, 'XTick', date_ticks, 'XTickLabel', {});   % ticks yes, dates no
h1 = colorbar(ax1, 'Position',[CBX Y1 CBW PH]);
set(h1,'FontName',FONT,'FontSize',FS);
ylabel(h1,'SSHA [m]', 'FontName',FONT,'FontSize',FS);

%% PANEL (b): phase-velocity bars aligned in time
ax2 = axes('Position',[LX Y2 PW PH]); hold(ax2, 'on');
x_shade = [regular_time_vec(1), regular_time_vec(end)];
fill_h = fill(ax2, [x_shade, fliplr(x_shade)], [mean_v + std_v, mean_v + std_v, mean_v - std_v, mean_v - std_v], ...
              [0 0 0], 'FaceAlpha', 0.16, 'EdgeColor', 'none');
mean_h = yline(ax2, mean_v, '--', 'LineWidth', 1.5, 'Color', [0 0 0]);
bar_h = bar(ax2, mean_x, v_phase, 1.5, 'FaceColor', [26, 224, 217]/255, 'EdgeColor', 'k');
bar_h.FaceAlpha = 0.85;
hold(ax2, 'off');
title(ax2, '(b)', 'FontName',FONT,'FontSize',FS);
ylabel(ax2, 'Phase Velocity [m/s]', 'FontName',FONT,'FontSize',FS);
xlim(ax2, [regular_time_vec(1) regular_time_vec(end)]);
ylim(ax2, [0 25]);
set_axes_style(ax2, FONT, FS);
set(ax2, 'XTick', date_ticks, 'XTickLabel', date_labels);   % manual date labels
xtickangle(ax2, 90);
lgd = legend(ax2, [bar_h, mean_h, fill_h], ...
             'Group phase velocity', ...
             sprintf('Mean: %.2f m/s', mean_v), ...
             sprintf('Std. Dev.: \\pm%.2f m/s', std_v), ...
             'Location', 'northeast');
set(lgd, 'FontName', FONT, 'FontSize', FS-5);
drawnow;

% -opengl: forces print to use the exact screen buffer (avoids the R2026a
% hardcopy font-relocation bug)
fig_filename = 'Fig06.png';
print(fig, '-dpng', fig_filename, '-r600', '-opengl');
close(fig);
fprintf('Figure saved: %s\n', fig_filename);
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

            if sum(valid_idx) < 20
                continue;
            end
            ts_filtered = ts;

            % 1. Linear trend
            if sum(valid_idx) >= 10
                p = polyfit(t_days(valid_idx), ts(valid_idx), 1);
                trend = polyval(p, t_days);
                trend = trend(:);
                ts_filtered = ts - trend;
            end

            % 2. Low frequency
            if sum(valid_idx) >= 15
                A = [ones(T, 1)];
                if T >= 30
                    omega_30 = 2*pi/30;
                    A = [A, cos(omega_30 * t_days(:)), sin(omega_30 * t_days(:))];
                end
                if ~isempty(A)
                    A_valid = A(valid_idx, :);
                    ts_valid = ts_filtered(valid_idx);
                    if size(A_valid, 1) > size(A_valid, 2)
                        coef = A_valid \ ts_valid;
                        seasonal_component = A * coef;
                        seasonal_component = seasonal_component(:);
                        ts_filtered = ts_filtered - seasonal_component;
                    end
                end
            end

            % 3. Monthly means
            if sum(valid_idx) >= 10
                date_vec = datevec(time_vector);
                months = date_vec(:, 2);
                unique_months = unique(months);
                for m = unique_months'
                    month_idx = (months == m) & valid_idx;
                    if sum(month_idx) >= 3
                        monthly_mean = nanmean(ts_filtered(month_idx));
                        ts_filtered(month_idx) = ts_filtered(month_idx) - monthly_mean;
                    end
                end
            end

            % 4. Temporal high-pass
            if sum(valid_idx) >= 30
                ts_interp = ts_filtered;
                if any(~valid_idx)
                    ts_interp(~valid_idx) = interp1(t_days(valid_idx), ...
                        ts_filtered(valid_idx), t_days(~valid_idx), 'linear', 'extrap');
                end
                try
                    [b, a] = butter(2, cutoff_freq/(fs/2), 'high');
                    ts_filtered_final = filtfilt(b, a, ts_interp);
                    ts_filtered_final = ts_filtered_final(:);
                    ts_filtered_final(~valid_idx) = NaN;
                    ts_filtered = ts_filtered_final;
                catch
                end
            end

            filtered_data(:, i, j) = ts_filtered;
        end
    end
end
function set_axes_style(ax, fontname, fontsize)
    ax.Layer   = 'top';
    ax.TickDir = 'in';
    box(ax, 'on');
    ax.FontName = fontname;
    ax.FontSize = fontsize;
    current_yticks = get(ax, 'YTick');
    y_limits = get(ax, 'YLim');
    ax.YTick = unique([y_limits(1), current_yticks, y_limits(2)]);
    ax.YTickLabel = abs(ax.YTick);
end
% =========================================================================
% LOCAL FUNCTIONS: SPECTRAL ANALYSIS (LOMB-SCARGLE) AND ROBUST STATISTICS
% =========================================================================
function [P_per_lat, Pxx_mean] = band_periods_by_latitude(HOV, tvec, f)
    nlat      = size(HOV, 1);
    P_per_lat = nan(nlat, 1);
    Pxx_acc   = zeros(numel(f), 1);
    cnt       = 0;
    for r = 1:nlat
        x = HOV(r, :);
        [P_per_lat(r), Pxx] = ls_dominant_period(tvec, x, f);
        if all(~isnan(Pxx)) && max(Pxx) > 0
            Pxx_acc = Pxx_acc + Pxx / max(Pxx);
            cnt     = cnt + 1;
        end
    end
    Pxx_mean = Pxx_acc / max(cnt, 1);
end
function [P_peak, Pxx] = ls_dominant_period(t, x, f)
    good = ~isnan(x) & ~isnan(t);
    t = t(good); x = x(good);
    if numel(x) < 8
        P_peak = NaN; Pxx = nan(numel(f), 1); return;
    end
    t = t(:) - t(1);
    x = x(:) - mean(x);
    Pxx = zeros(numel(f), 1);
    for k = 1:numel(f)
        w   = 2*pi*f(k);
        tau = atan2(sum(sin(2*w*t)), sum(cos(2*w*t))) / (2*w);
        wt  = w*(t - tau);
        cwt = cos(wt); swt = sin(wt);
        cc  = sum(cwt.^2); ss = sum(swt.^2);
        t1 = 0; t2 = 0;
        if cc > eps, t1 = (sum(x.*cwt))^2 / cc; end
        if ss > eps, t2 = (sum(x.*swt))^2 / ss; end
        Pxx(k) = 0.5*(t1 + t2);
    end
    [~, idx] = max(Pxx);
    P_peak = 1 / f(idx);
end
function [Pmat, frac_band, frac_char] = per_latitude_spectral_diagnostics( ...
                HOV, tvec, f_full, nPeaks, fb_lo, fb_hi, Tc, rel)
    nlat = size(HOV, 1);
    Pmat = nan(nlat, nPeaks);
    frac_band = nan(nlat, 1);
    frac_char = nan(nlat, 1);
    in_band = (f_full >= fb_lo) & (f_full <= fb_hi);
    fc_lo = 1/(Tc*(1+rel));  fc_hi = 1/(Tc*(1-rel));
    in_char = (f_full >= fc_lo) & (f_full <= fc_hi);
    for r = 1:nlat
        x = HOV(r, :);
        [~, Pxx] = ls_dominant_period(tvec, x, f_full);
        if any(isnan(Pxx)) || sum(Pxx) <= 0
            continue;
        end
        tot = sum(Pxx);
        frac_band(r) = sum(Pxx(in_band)) / tot;
        frac_char(r) = sum(Pxx(in_char)) / tot;
        locs = local_maxima_indices(Pxx);
        if ~isempty(locs)
            [~, order] = sort(Pxx(locs), 'descend');
            locs = locs(order);
            nTake = min(nPeaks, numel(locs));
            Pmat(r, 1:nTake) = 1 ./ f_full(locs(1:nTake));
        end
    end
end
function Pmat = per_latitude_peak_periods(HOV, tvec, f, nPeaks)
    nlat = size(HOV, 1);
    Pmat = nan(nlat, nPeaks);
    for r = 1:nlat
        x = HOV(r, :);
        [~, Pxx] = ls_dominant_period(tvec, x, f);
        if any(isnan(Pxx)) || max(Pxx) <= 0
            continue;
        end
        locs = local_maxima_indices(Pxx);
        if isempty(locs)
            continue;
        end
        [~, order] = sort(Pxx(locs), 'descend');
        locs = locs(order);
        nTake = min(nPeaks, numel(locs));
        Pmat(r, 1:nTake) = 1 ./ f(locs(1:nTake));
    end
end
function locs = local_maxima_indices(Pxx)
    Pxx = Pxx(:);
    n = numel(Pxx);
    locs = [];
    for i = 1:n
        leftOK  = (i == 1) || (Pxx(i) >  Pxx(i-1));
        rightOK = (i == n) || (Pxx(i) >= Pxx(i+1));
        if leftOK && rightOK
            locs(end+1) = i;
        end
    end
end
function print_peak_rank_table(Pmat, label)
    fprintf('  -- %s --\n', label);
    nP = size(Pmat, 2);
    for k = 1:nP
        col = Pmat(~isnan(Pmat(:, k)), k);
        if isempty(col)
            fprintf('    Peak #%d: (not enough maxima across latitudes)\n', k);
        else
            fprintf('    Peak #%d: median T = %6.2f days   (IQR %.2f-%.2f d, n latitudes = %d)\n', ...
                k, median(col), manual_percentile(col,25), manual_percentile(col,75), numel(col));
        end
    end
end
function y = manual_percentile(x, p)
    % Percentile without the Statistics Toolbox (linear interpolation,
    % equivalent to the default method of prctile).
    x = sort(x(~isnan(x)));
    n = numel(x);
    if n == 0
        y = NaN; return;
    end
    if n == 1
        y = x(1); return;
    end
    pos = p/100*(n-1) + 1;
    lo  = floor(pos); hi = ceil(pos);
    lo  = max(min(lo, n), 1);
    hi  = max(min(hi, n), 1);
    frac = pos - lo;
    y = x(lo) + frac*(x(hi) - x(lo));
end
