%% Tidal Harmonic Validation: Root Sum Square (RSS) Misfit
% Evaluates global and regional tidal models against in-situ harmonic
% constituents derived from pyFES. Uses the common constituents across all
% datasets at each station, computes the total RSS misfit, and exports a
% grouped bar chart plus a summary table.

clear; close all; clc;

%% 1. Configuration and Parameters
dirs.input  = '.'; 
dirs.output = './RSSmisfit';
if ~exist(dirs.output, 'dir'); mkdir(dirs.output); end

params.fig_width_cm  = 15;
params.fig_height_cm = 12;
params.dpi           = 600;
params.font_name     = 'Times New Roman';
params.font_size     = 12;

% Model colors
colors = [
    195, 195, 195;  % FES2022
     26, 224, 217;  % TPXO10v2
    209, 148, 255;  % EOT20
    254, 222,  76;  % GOT5.5
    253, 154,  95   % MSAS (Regional)
] / 255;

model_names = {'FES2022', 'TPXO10v2', 'EOT20', 'GOT5.5', 'MSAS'};

%% 2. Station Definitions
% Filenames and CSV search keys for each station
stations = struct(...
    'name', {'Carina', 'Lander', 'P. Deseado', 'Vega-Pleyade', 'San Matias'}, ...
    'csv_key', {'carina', 'lander', 'puerto deseado', 'vega-pleyade', 'sbe26 san matias'}, ...
    'pyfes', {'tidal_pyFES_components_carina.xlsx', ...
              'tidal_pyFES_components_Lander.xlsx', ...
              'tidal_pyFES_components_PDeseado.xlsx', ...
              'tidal_pyFES_components_Vega.xlsx', ...
              'tidal_pyFES_components_sbe26_sm.xlsx'}, ...
    'fes',   {'tides_constituents_Carina_FES2022.xlsx', ...
              'tides_constituents_Lander_FES2022.xlsx', ...
              'tides_constituents_PDeseado_FES2022v2.xlsx', ... % Note the v2
              'tides_constituents_Vega_FES2022.xlsx', ...
              'tides_constituents_sbe26SM_FES2022.xlsx'}, ...
    'tpxo',  {'tides_constituents_Carina_TPXO10v2.xlsx', ...
              'tides_constituents_Lander_TPXO10v2.xlsx', ...
              'tides_constituents_PDeseado_TPXO10v2.xlsx', ...
              'tides_constituents_Vega_TPXO10v2.xlsx', ...
              'tides_constituents_sbe26SM_TPXO10v2.xlsx'}, ...
    'eot',   {'tides_constituents_Carina_EOT20.xlsx', ...
              'tides_constituents_Lander_EOT20.xlsx', ...
              'tides_constituents_PDeseado_EOT20.xlsx', ...
              'tides_constituents_Vega_EOT20.xlsx', ...
              'tides_constituents_sbe26SM_EOT20.xlsx'}, ...
    'got',   {'tides_constituents_Carina_GOT5_5.xlsx', ...
              'tides_constituents_Lander_GOT5_5.xlsx', ...
              'tides_constituents_PDeseado_GOT5_5.xlsx', ...
              'tides_constituents_Vega_GOT5_5.xlsx', ...
              'tides_constituents_sbe26SM_GOT5_5.xlsx'} ...
);

num_stations = numel(stations);
num_models   = numel(model_names);
rss_results  = nan(num_stations, num_models);

% Analytical formulation for Root Sum Square (RSS) misfit
calc_rss = @(H1, g1, H2, g2) sqrt(0.5 * sum((H1.*cos(g1) - H2.*cos(g2)).^2 + ...
                                            (H1.*sin(g1) - H2.*sin(g2)).^2));
deg2rad_local = @(x) x * (pi/180);

%% 3. Load Regional Model (CSV)
csv_file = fullfile(dirs.input, 'armonicos_estaciones_modelo.csv');
fprintf('Loading regional model database: %s\n', csv_file);
T_reg_full = readtable(csv_file, 'VariableNamingRule', 'preserve');
reg_vars = T_reg_full.Properties.VariableNames;
reg_constituents = lower(string(T_reg_full{:, 1}));

%% 4. Iterative Processing per Station
fprintf('\nComputing RSS misfits...\n');

for s = 1:num_stations
    st = stations(s);
    fprintf('Processing station: %s\n', st.name);
    
    % 4.1 Load external data
    T_pyfes = readtable(fullfile(dirs.input, st.pyfes), 'Sheet', 'Complete');
    T_fes   = readtable(fullfile(dirs.input, st.fes));
    T_tpxo  = readtable(fullfile(dirs.input, st.tpxo));
    T_eot   = readtable(fullfile(dirs.input, st.eot));
    T_got   = readtable(fullfile(dirs.input, st.got));
    
    % 4.2 Extract target columns from regional CSV
    csv_cols = find(contains(lower(reg_vars), st.csv_key));
    if numel(csv_cols) < 2
        error('Insufficient columns found for station "%s" in regional CSV.', st.name);
    end
    
    med_vals = median(abs(T_reg_full{:, csv_cols}), 1, 'omitnan');
    amp_idx = csv_cols(find(med_vals < 30, 1, 'first'));
    pha_idx = csv_cols(find(~ismember(1:numel(csv_cols), find(csv_cols == amp_idx)), 1, 'first'));
    if isempty(amp_idx) || isempty(pha_idx)
        amp_idx = csv_cols(1); pha_idx = csv_cols(2);
    end
    
    amp_reg_all = T_reg_full{:, amp_idx} * 100; % Convert m to cm
    pha_reg_all = wrapTo360(T_reg_full{:, pha_idx});
    
    % 4.3 Determine intersection of valid constituents
    set_pyfes = lower(string(T_pyfes.constituent));
    set_fes   = lower(string(T_fes.constituent));
    set_tpxo  = lower(string(T_tpxo.constituent));
    set_eot   = lower(string(T_eot.constituent));
    set_got   = lower(string(T_got.constituent));
    
    common_set = intersect(set_pyfes, set_fes);
    common_set = intersect(common_set, set_tpxo);
    common_set = intersect(common_set, set_eot);
    common_set = intersect(common_set, set_got);
    common_set = intersect(common_set, reg_constituents);
    
    num_common = numel(common_set);
    
    % Allocate arrays for common constituents
    amp_ref = nan(num_common, 1); pha_ref = nan(num_common, 1);
    amp_mod = nan(num_common, num_models); pha_mod = nan(num_common, num_models);
    
    % 4.4 Extract amplitude and phase for matched constituents
    for i = 1:num_common
        c = common_set(i);
        
        idx = strcmpi(T_pyfes.constituent, c);
        amp_ref(i) = T_pyfes.amplitude(idx) * 100;
        pha_ref(i) = wrapTo360(T_pyfes.phase_deg(idx));
        
        idx = strcmpi(T_fes.constituent, c);
        amp_mod(i, 1) = T_fes.amplitude(idx) * 100;
        pha_mod(i, 1) = wrapTo360(T_fes.phase_deg(idx));
        
        idx = strcmpi(T_tpxo.constituent, c);
        amp_mod(i, 2) = T_tpxo.amplitude(idx) * 100;
        pha_mod(i, 2) = wrapTo360(T_tpxo.phase_deg(idx));
        
        idx = strcmpi(T_eot.constituent, c);
        amp_mod(i, 3) = T_eot.amplitude(idx) * 100;
        pha_mod(i, 3) = wrapTo360(T_eot.phase_deg(idx));
        
        idx = strcmpi(T_got.constituent, c);
        amp_mod(i, 4) = T_got.amplitude(idx) * 100;
        pha_mod(i, 4) = wrapTo360(T_got.phase_deg(idx));
        
        idx = strcmpi(reg_constituents, c);
        amp_mod(i, 5) = amp_reg_all(idx);
        pha_mod(i, 5) = pha_reg_all(idx);
    end
    
    % 4.5 Compute Total RSS
    rad_ref = deg2rad_local(pha_ref);
    rad_mod = deg2rad_local(pha_mod);
    
    for m = 1:num_models
        rss_results(s, m) = calc_rss(amp_mod(:, m), rad_mod(:, m), amp_ref, rad_ref);
    end
end

%% 5. Global Visualization
fprintf('\nGenerating figure...\n');

fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', ...
             'Position', [1, 1, params.fig_width_cm, params.fig_height_cm]);
set(fig, 'PaperUnits', 'centimeters', 'PaperSize', [params.fig_width_cm, params.fig_height_cm]);
set(fig, 'PaperPositionMode', 'manual', 'PaperPosition', [0 0 params.fig_width_cm params.fig_height_cm]);

ax = axes('Parent', fig);
b = bar(ax, 1:num_stations, rss_results, 'grouped');

for k = 1:num_models
    b(k).FaceColor = colors(k, :);
end

set(ax, 'XTick', 1:num_stations, 'XTickLabel', {stations.name}, ...
    'FontName', params.font_name, 'FontSize', params.font_size, ...
    'TickDir', 'in', 'Box', 'on');
ylabel(ax, 'Total RSS Misfit (cm)', 'FontName', params.font_name, 'FontSize', params.font_size, 'FontWeight', 'bold');
title(ax, 'Total Tidal Harmonic RSS Misfit vs pyFES', 'FontName', params.font_name, 'FontSize', params.font_size + 2, 'FontWeight', 'bold');

lgd = legend(ax, model_names, 'Location', 'northeast', 'Orientation', 'horizontal');
set(lgd, 'FontName', params.font_name, 'FontSize', params.font_size - 1, 'Box', 'off');

ylim_max = max(rss_results(:)) * 1.15;
ylim(ax, [0, ylim_max]);

% Annotate bar values
for s = 1:num_stations
    for m = 1:num_models
        text(b(m).XEndPoints(s), b(m).YEndPoints(s) + ylim_max*0.02, ...
            sprintf('%.2f', rss_results(s, m)), 'HorizontalAlignment', 'center', ...
            'FontName', params.font_name, 'FontSize', params.font_size - 2, 'Rotation', 90);
    end
end

% Export figure
fig_filename = fullfile(dirs.output, 'Fig02.png');
print(fig, fig_filename, '-dpng', sprintf('-r%d', params.dpi));
fprintf('Figure saved: %s\n', fig_filename);
close(fig);

%% 6. Tabular Export
T_export = array2table(rss_results, 'VariableNames', model_names);
T_export = [table({stations.name}', 'VariableNames', {'Station'}), T_export];
excel_filename = fullfile(dirs.output, 'Table_Global_Total_RSS_Misfit.xlsx');
writetable(T_export, excel_filename);
fprintf('Table saved: %s\n', excel_filename);