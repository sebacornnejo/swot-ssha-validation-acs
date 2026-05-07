%% SWOT and In-Situ Sea Level Anomaly Validation
% Description:
%   Evaluates Sea Surface Height Anomalies (SSHA) from in-situ platforms 
%   and SWOT altimetry. The input Sea Surface Height (SSH) records for both 
%   datasets have been previously corrected for tides using reconstructed 
%   tides from the FES2022b model. This script specifically evaluates the 
%   performance of different Dynamic Atmospheric Corrections (DACs). It 
%   centers the anomalies over the common temporal overlap and computes 
%   standard validation metrics (r, RMSE, R-squared) for the SSHA corrected 
%   by each respective DAC.

clear; close all; clc;

%% Configuration and Directories
input_dir   = '.';
qc_data_dir = './QC_Data';
stats_dir   = './Statistics';
figs_dir    = './Figures';

if ~exist(qc_data_dir, 'dir'); mkdir(qc_data_dir); end
if ~exist(stats_dir, 'dir'); mkdir(stats_dir); end
if ~exist(figs_dir, 'dir'); mkdir(figs_dir);  end

% Quality control thresholds
config.max_variance_iterations = 5;

% Plotting and export parameters
FIG_WIDTH_CM   = 15;      
FIG_HEIGHT_CM  = 12;      
FIG_DPI        = 600;     
FONT_NAME      = 'Times New Roman';
FONT_SIZE      = 10;      
COLOR_INSITU   = [157, 157, 157] / 255;   
COLOR_SWOT     = [126, 79, 255] / 255;   

% Input files and correction versions
input_files = { ...
    './CARINA_FES2022vFES2022.xlsx', ...
    './LANDER_FES2022vFES2022.xlsx', ...
    './PDESEADO_FES2022vFES2022.xlsx', ...
    './SBE26_FES2022vFES2022.xlsx', ...
    './VEGA_FES2022vFES2022.xlsx'};
dac_versions = {'Trad', 'CalVal', 'MSAS', 'SIROCCO'};
stats_file = fullfile(stats_dir, 'Validation_Statistics_Summary.xlsx');

%% Main Processing Loop
for i_file = 1:length(input_files)
    current_file = input_files{i_file};
    full_path    = fullfile(input_dir, current_file);
    
    [~, station_name, ~] = fileparts(current_file);
    station_name = strrep(station_name, '_FES2022vFES2022', '');
    
    fprintf('\n------------------------------------------------------------\n');
    fprintf('Processing station: %s (%d of %d)\n', upper(station_name), i_file, length(input_files));
    fprintf('------------------------------------------------------------\n');
    
    qc_output_file = fullfile(qc_data_dir, sprintf('%s_CLEANED.xlsx', station_name));
    
    station_results = table('Size', [0, 4], ...
        'VariableTypes', {'string', 'double', 'double', 'double'}, ...
        'VariableNames', {'DAC', 'Corr_r', 'RMSE', 'R2_percent'});
    
    % Initialize figure with precise dimensions
    fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', ...
                 'Position', [1, 1, FIG_WIDTH_CM, FIG_HEIGHT_CM]);
             
    set(fig, 'PaperUnits', 'centimeters', 'PaperSize', [FIG_WIDTH_CM, FIG_HEIGHT_CM]);
    set(fig, 'PaperPositionMode', 'manual', 'PaperPosition', [0 0 FIG_WIDTH_CM FIG_HEIGHT_CM]);
    
    % Define specific panel layout coordinates
    m_left   = 0.09;  
    m_right  = 0.02;  
    m_bottom = 0.18; 
    m_top    = 0.12; 
    gap_x    = 0.07; 
    gap_y    = 0.09; 
    
    w_ax = (1 - m_left - m_right - gap_x) / 2;
    h_ax = (1 - m_bottom - m_top - gap_y) / 2;
    
    panel_positions = {
        [m_left, m_bottom + h_ax + gap_y, w_ax, h_ax], ...               
        [m_left + w_ax + gap_x, m_bottom + h_ax + gap_y, w_ax, h_ax], ...
        [m_left, m_bottom, w_ax, h_ax], ...                              
        [m_left + w_ax + gap_x, m_bottom, w_ax, h_ax]                    
    };
    
    global_max_abs = 0;
    ax_handles     = gobjects(1, length(dac_versions));
    
    % Impose bounding temporal limits per station
    if i_file == 2
        xlim_station = [datenum(2023,4,11,0,0,0), datenum(2023,7,11,23,59,59)];
    else
        xlim_station = [datenum(2023,3,27,0,0,0), datenum(2023,7,11,23,59,59)];
    end
    
    for i_dac = 1:length(dac_versions)
        current_dac = dac_versions{i_dac};
        
        % 1. Data loading
        try
            opts = detectImportOptions(full_path, 'Sheet', current_dac, 'ReadVariableNames', true);
            varNames = opts.VariableNames;
            for k = 2:length(varNames); opts = setvartype(opts, varNames{k}, 'double'); end
            opts.MissingRule = 'fill'; opts.ImportErrorRule = 'fill';
            dataTable = readtable(full_path, opts);
        catch
            continue;
        end
        
        % 2. Variance-based quality control (spike removal)
        cols_to_check = [2, 3];
        [dataTable_qc, ~] = apply_variance_qc(dataTable, cols_to_check, config);
        
        % 2.1 Discard initial transient common observation
        insitu_col = 2; swot_col = 3;
        mask_common_initial = ~isnan(dataTable_qc{:, insitu_col}) & ~isnan(dataTable_qc{:, swot_col});
        idx_first_common    = find(mask_common_initial, 1, 'first');
        if ~isempty(idx_first_common)
            dataTable_qc{idx_first_common, insitu_col} = NaN;
            dataTable_qc{idx_first_common, swot_col}   = NaN;
        end
        
        % 2.2 Sea Level Anomaly (SLA) centering based on common period
        mask_common = ~isnan(dataTable_qc{:, insitu_col}) & ~isnan(dataTable_qc{:, swot_col});
        if sum(mask_common) > 0
            mean_insitu_common = mean(dataTable_qc{mask_common, insitu_col});
            mean_swot_common   = mean(dataTable_qc{mask_common, swot_col});
            dataTable_qc{:, insitu_col} = dataTable_qc{:, insitu_col} - mean_insitu_common;
            dataTable_qc{:, swot_col}   = dataTable_qc{:, swot_col}   - mean_swot_common;
        end
        
        % 3. Error metrics computation
        valid_mask   = ~isnan(dataTable_qc{:, insitu_col}) & ~isnan(dataTable_qc{:, swot_col});
        insitu_clean = dataTable_qc{valid_mask, insitu_col};
        swot_clean   = dataTable_qc{valid_mask, swot_col};
        
        if numel(insitu_clean) > 2
            r      = corr(insitu_clean, swot_clean);
            rmse   = sqrt(mean((swot_clean - insitu_clean).^2));
            r2_pct = (r^2) * 100;
        else
            r = NaN; rmse = NaN; r2_pct = NaN;
        end
        station_results = [station_results; {current_dac, r, rmse, r2_pct}]; 
        
        % 4. Export quality-controlled data
        writetable(dataTable_qc, qc_output_file, 'Sheet', current_dac);
        
        % 5. Construct time-series panels
        ax_handles(i_dac) = axes('Parent', fig, 'Position', panel_positions{i_dac});
        
        time_vector = dataTable_qc{:, 1};
        insitu_plot = dataTable_qc{:, insitu_col};
        swot_plot   = dataTable_qc{:, swot_col};
        
        idx_valid_insitu = ~isnan(insitu_plot);
        idx_valid_swot   = ~isnan(swot_plot);
        
        max_in = max(abs(insitu_plot(idx_valid_insitu)));
        max_sw = max(abs(swot_plot(idx_valid_swot)));
        if isempty(max_in); max_in = 0; end
        if isempty(max_sw); max_sw = 0; end
        global_max_abs = max([global_max_abs, max_in, max_sw]);
        
        ax = ax_handles(i_dac);

        x_grid = min(time_vector)-7 : max(time_vector)+7;
        grid_offsets = [0, 0.5, -0.5, 1, -1, 1.5, -1.5];
        
        hold(ax, 'on')
        for k = 1:numel(grid_offsets)
            plot(ax, x_grid, grid_offsets(k)*ones(size(x_grid)), ...
                '-', 'Color', [.9 .9 .9], 'LineWidth', 0.4);
        end

        h_in = plot(ax_handles(i_dac), time_vector(idx_valid_insitu), insitu_plot(idx_valid_insitu), ...
             '-', 'Color', COLOR_INSITU, 'LineWidth', 0.5, 'DisplayName', 'In Situ');
        hold(ax_handles(i_dac), 'on');
        h_sw = plot(ax_handles(i_dac), time_vector(idx_valid_swot), swot_plot(idx_valid_swot), ...
             '-o', 'Color', COLOR_SWOT, 'MarkerSize', 1.8, 'MarkerFaceColor', COLOR_SWOT, ...
             'LineWidth', 0.5, 'DisplayName', 'SWOT');
         
        set(ax_handles(i_dac), 'FontName', FONT_NAME, 'FontSize', FONT_SIZE, 'Box', 'on', 'TickDir', 'in');
        
        labels = {'(a)', '(b)', '(c)', '(d)'};
        prefix = labels{i_dac};

        if strcmpi(current_dac, 'Trad')
            plot_label = 'Trad.';
        else
            plot_label = current_dac;
        end

        if ~isnan(r)
            title_str = sprintf('\\bf%s \\rm DAC %s\\rm\n{\\itr} = %.2f | RMSE = %.3f m', ...
                prefix, plot_label, r, rmse);
        else
            title_str = sprintf('\\bf%s\\rm\nInsufficient data', plot_label);
        end

        t = title(ax_handles(i_dac), title_str, ...
            'FontName', FONT_NAME, ...
            'FontSize', FONT_SIZE, ...
            'Interpreter', 'tex');
        
        t.Units = 'normalized';          
        pos = t.Position;
        pos(2) = pos(2) + 0.004;         
        t.Position = pos;
        
        if i_dac == 1 || i_dac == 2
            xticklabels(ax_handles(i_dac), {});
        end
    end
    
    % 5.5 Apply uniform axes limits and date ticks
    ylim_val = ceil(global_max_abs * 10) / 10;
    if ylim_val == 0; ylim_val = 0.1; end
    xtick_vals = compute_date_ticks(xlim_station(1), xlim_station(2));
    
    for k = 1:length(ax_handles)
        if ~isgraphics(ax_handles(k)); continue; end
        ax_k = ax_handles(k);
        ylim(ax_k, [-ylim_val, ylim_val]);
        xlim(ax_k, xlim_station);
        set(ax_k, 'XTick', xtick_vals);
        
        if k == 3 || k == 4
            set(ax_k, 'XTickLabel', datestr(xtick_vals, 'dd-mmm-yy'));
            xtickangle(ax_k, 45);
        else
            set(ax_k, 'XTickLabel', []);
        end
    end
    
    % 6. Global figure aesthetics
    ax_ghost = axes('Position', [0 0 1 1], 'Visible', 'off', 'HitTest', 'off');
    
    text(ax_ghost, 0.02, 0.5, 'SSHA (m)', 'Rotation', 90, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontName', FONT_NAME, 'FontSize', FONT_SIZE + 1);
        
    station_map = containers.Map( ...
        {'SBE26','PDESEADO','VEGA'}, ...
        {'G. San Matias','P. Deseado','Vega-Pleyade'} );
    
    map_key = upper(station_name);
    
    if isKey(station_map, map_key)
        plot_name = station_map(map_key);
    else
        plot_name = station_name; 
    end
    
    text(ax_ghost, 0.5, 1, sprintf('SSHA Comparison — Station: %s', upper(plot_name)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
        'FontName', FONT_NAME, 'FontSize', FONT_SIZE + 2, 'FontWeight', 'bold');
        
    lgd = legend([h_in, h_sw], 'Orientation', 'horizontal');
    set(lgd, 'Position', [0.40, 0.005, 0.20, 0.05], 'Units', 'normalized', ...
             'FontName', FONT_NAME, 'FontSize', FONT_SIZE, 'Box', 'off');
    
    % 7. Export figure
    fig_filename = fullfile(figs_dir, sprintf('%s_Series.png', station_name));
    print(fig, fig_filename, '-dpng', sprintf('-r%d', FIG_DPI));
    fprintf('Figure saved: %s\n', fig_filename);
    close(fig);
    
    if ~isempty(station_results)
        writetable(station_results, stats_file, 'Sheet', station_name);
    end
end
fprintf('\n================== PROCESSING COMPLETED ==================\n');

%% Local Functions
function tick_vals = compute_date_ticks(t_min, t_max)
    dv_min = datevec(t_min); dv_max = datevec(t_max);
    yr1 = dv_min(1); yr2 = dv_max(1);
    mo1 = dv_min(2); mo2 = dv_max(2);
    
    interior = [];
    for y = yr1:yr2
        m_start = 1; m_end = 12;
        if y == yr1; m_start = mo1; end
        if y == yr2; m_end = mo2; end
        for m = m_start:m_end
            interior(end+1) = datenum(y, m, 1);  
            interior(end+1) = datenum(y, m, 15); 
        end
    end
    
    min_dist_days = 8;
    interior_filtered = interior( (interior > t_min + min_dist_days) & (interior < t_max - min_dist_days) );
    tick_vals = unique([t_min, interior_filtered, t_max]);
    tick_vals = sort(tick_vals);
end

function [dataTable_qc, outlier_info_table] = apply_variance_qc(dataTable, cols_to_check, config)
    dataTable_qc = dataTable;
    outlier_candidates = table('Size', [0, 3], 'VariableTypes', {'double', 'double', 'double'}, 'VariableNames', {'RowIndex', 'ColIndex', 'DeviationValue'});
    last_removed_deviation  = inf;
    halt_due_to_homogeneity = false;
    outlier_info_table      = table();
    var_names = dataTable.Properties.VariableNames;
    
    for iter = 1:config.max_variance_iterations + 1
        max_dev_this_pass  = -1; candidate_row_idx  = -1; candidate_col_idx  = -1;
        for k = 1:length(cols_to_check)
            current_col = cols_to_check(k);
            current_data = dataTable_qc{:, current_col};
            valid_data_mask = ~isnan(current_data);
            if sum(valid_data_mask) < 3; continue; end
            
            data_clean  = current_data(valid_data_mask);
            data_mean   = mean(data_clean);
            data_std    = std(data_clean);
            deviations  = abs(current_data - data_mean);
            
            [current_max_dev, current_idx] = max(deviations);
            if current_max_dev > (2 * data_std) && current_max_dev > max_dev_this_pass
                max_dev_this_pass = current_max_dev;
                candidate_row_idx = current_idx; candidate_col_idx = current_col;
            end
        end
        
        if candidate_row_idx == -1; break; end
        if iter > 1 && (max_dev_this_pass <= (2/3) * last_removed_deviation); break; end
        if iter > config.max_variance_iterations; halt_due_to_homogeneity = true; break; end
        
        outlier_candidates = [outlier_candidates; {candidate_row_idx, candidate_col_idx, max_dev_this_pass}]; 
        last_removed_deviation = max_dev_this_pass;
        dataTable_qc{candidate_row_idx, candidate_col_idx} = NaN;
    end
    
    if ~halt_due_to_homogeneity && ~isempty(outlier_candidates)
        outlier_info_table = table('Size', [0, 3], 'VariableTypes', {'string', 'string', 'double'}, 'VariableNames', {'Date', 'Variable', 'Deviation'});
        for k = 1:height(outlier_candidates)
            row_to_remove = outlier_candidates.RowIndex(k);
            col_to_remove = outlier_candidates.ColIndex(k);
            dev_value     = outlier_candidates.DeviationValue(k);
            outlier_date  = dataTable{row_to_remove, 1};
            variable_name = var_names{col_to_remove};
            if isnumeric(outlier_date); date_str = datestr(outlier_date, 'dd-mmm-yyyy'); else; date_str = char(outlier_date); end
            outlier_info_table = [outlier_info_table; {date_str, variable_name, dev_value}]; 
            dataTable{row_to_remove, col_to_remove} = NaN;
        end
    end
end