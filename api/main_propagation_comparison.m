% MAIN_PROPAGATION_COMPARISON
% Master Harness for Unified Radio Propagation Computational Framework.
%
% Objectives:
%   1. Load terrain profile.
%   2. Execute multiple propagation models via standardized wrappers.
%   3. Compare results via plots and metrics (Reference: EFIE).
%
% Pre-requisites:
%   - Run setup_project.m first.
%   - Ensure 'data/terrain_profile.txt' exists.

clear; clc; close all;

%% 1. Configuration & Scenario Setup

% --- Simulation Parameters ---
Scenario = struct();
Scenario.Frequency_MHz = 970;       % Frequency
Scenario.Tx_Height_m = 52.0;        % Tx Antenna Height (AGL)
Scenario.Rx_Height_m = 2.4;         % Rx Antenna Height (AGL)
Scenario.Environment = 'suburban';  % For Hata model
Scenario.Run_EFIE = true;           % EFIE is computationally expensive (set false for quick tests)

% --- File Paths ---
Scenario.Terrain_File_Path = fullfile('data', 'terrain_profile.txt');
Scenario.Output_Dir = 'output';

% --- Model Selection Flags ---
DO_ITM = true;
DO_HATA = true;
DO_BULLINGTON = true;
DO_DEYGOUT = true;

% Ensure output directory exists
if ~exist(Scenario.Output_Dir, 'dir'), mkdir(Scenario.Output_Dir); end

fprintf('========================================================\n');
fprintf('   Unified Propagation Model Comparison Harness\n');
fprintf('========================================================\n');
fprintf('Frequency:      %.1f MHz\n', Scenario.Frequency_MHz);
fprintf('Tx/Rx Height:   %.1f m / %.1f m\n', Scenario.Tx_Height_m, Scenario.Rx_Height_m);
fprintf('Environment:    %s\n', Scenario.Environment);
fprintf('Ref Model:      EFIE (Full-Wave)\n');
fprintf('========================================================\n');

%% 2. Data Loading

if ~isfile(Scenario.Terrain_File_Path)
    % Fallback for legacy setups
    if isfile('X.txt')
        fprintf('Warning: Standard data file missing. Using legacy "X.txt".\n');
        Scenario.Terrain_File_Path = 'X.txt';
    else
        error('Terrain file not found at "%s". Run dem_profile_api.py first.', Scenario.Terrain_File_Path);
    end
end

% Load Data (Format: Distance(m) | Elevation(m))
raw_data = readmatrix(Scenario.Terrain_File_Path);

% Populate Terrain Profile Struct
TerrainProfile = struct();
TerrainProfile.Distance_m = raw_data(:, 1);
TerrainProfile.Elevation_m = raw_data(:, 2);

% Sanitize Data (Remove NaNs)
valid_mask = ~isnan(TerrainProfile.Elevation_m) & ~isnan(TerrainProfile.Distance_m);
TerrainProfile.Distance_m = TerrainProfile.Distance_m(valid_mask);
TerrainProfile.Elevation_m = TerrainProfile.Elevation_m(valid_mask);

dist_km = TerrainProfile.Distance_m / 1000.0;
num_points = length(dist_km);

fprintf('Loaded Terrain: %d points, %.2f km total length.\n', num_points, dist_km(end));

%% 3. Model Execution

Results = struct();
Results.Distance_km = dist_km;

% Add wrappers to path if strictly necessary (usually handled by addpath in setup) 
addpath(genpath('models')); 

% --- A. EFIE (Reference Model) ---
if Scenario.Run_EFIE
    fprintf('Running EFIE Model (Reference)... ');
    tic;
    Results.EFIE_dB = efie_wrapper(Scenario, TerrainProfile);
    t = toc;
    fprintf('Done (%.2f s).\n', t);
else
    fprintf('Skipping EFIE (Reference unavailable).\n');
    Results.EFIE_dB = nan(size(dist_km));
end

% --- B. Longley-Rice (ITM) ---
if DO_ITM
    fprintf('Running Longley-Rice (ITM)...     ');
    tic;
    Results.ITM_dB = itm_wrapper(Scenario, TerrainProfile);
    t = toc;
    fprintf('Done (%.2f s).\n', t);
end

% --- C. Hata / COST 231 ---
if DO_HATA
    fprintf('Running Hata/COST231...           ');
    tic;
    Results.Hata_dB = hata_wrapper(Scenario, TerrainProfile);
    t = toc;
    fprintf('Done (%.2f s).\n', t);
end

% --- D. Bullington Diffraction ---
if DO_BULLINGTON
    fprintf('Running Bullington...             ');
    tic;
    Results.Bullington_dB = bullington_wrapper(Scenario, TerrainProfile);
    t = toc;
    fprintf('Done (%.2f s).\n', t);
end

% --- E. Deygout Diffraction ---
if DO_DEYGOUT
    fprintf('Running Deygout...                ');
    tic;
    Results.Deygout_dB = deygout_wrapper(Scenario, TerrainProfile);
    t = toc;
    fprintf('Done (%.2f s).\n', t);
end

%% 4. Visualization

fprintf('\nGenerating Plots...\n');

% --- Figure 1: Path Loss Comparison ---
fig1 = figure('Name', 'Path Loss Comparison', 'Color', 'w', 'Position', [100, 100, 1000, 600]);
hold on; grid on; box on;

% Plot order: Reference last (or prominent)
if DO_HATA, plot(dist_km, Results.Hata_dB, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Hata/COST231'); end
if DO_ITM, plot(dist_km, Results.ITM_dB, 'k-', 'LineWidth', 2, 'DisplayName', 'ITM (Longley-Rice)'); end
if DO_BULLINGTON, plot(dist_km, Results.Bullington_dB, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Bullington'); end
if DO_DEYGOUT, plot(dist_km, Results.Deygout_dB, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Deygout'); end
if Scenario.Run_EFIE, plot(dist_km, Results.EFIE_dB, 'r:', 'LineWidth', 2.5, 'DisplayName', 'EFIE (Full-Wave)'); end

xlabel('Distance (km)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Path Loss (dB)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Propagation Model Comparison (f=%.1f MHz)', Scenario.Frequency_MHz), 'FontSize', 14);
legend('Location', 'Best', 'FontSize', 10);

% Smart Y-Axis Limits
all_losses = [];
if DO_ITM, all_losses = [all_losses; Results.ITM_dB]; end
if Scenario.Run_EFIE, all_losses = [all_losses; Results.EFIE_dB]; end
if ~isempty(all_losses)
    valid_vals = all_losses(~isnan(all_losses) & all_losses > 0);
    if ~isempty(valid_vals)
        ylim([min(valid_vals)-10, max(valid_vals)+20]);
    end
end

saveas(fig1, fullfile(Scenario.Output_Dir, 'comparison_plot.png'));

% --- Figure 2: Terrain Geometry ---
fig2 = figure('Name', 'Terrain Geometry', 'Color', 'w', 'Position', [150, 150, 1000, 400]);
area(dist_km, TerrainProfile.Elevation_m, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'k', 'DisplayName', 'Terrain');
hold on; grid on; box on;

% Tx/Rx Markers
tx_x = dist_km(1); tx_y = TerrainProfile.Elevation_m(1);
rx_x = dist_km(end); rx_y = TerrainProfile.Elevation_m(end);
tx_ant_h = tx_y + Scenario.Tx_Height_m;
rx_ant_h = rx_y + Scenario.Rx_Height_m;

plot([tx_x tx_x], [tx_y tx_ant_h], 'k-', 'LineWidth', 2);
plot([rx_x rx_x], [rx_y rx_ant_h], 'k-', 'LineWidth', 2);
plot(tx_x, tx_ant_h, 'b^', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'DisplayName', 'Tx');
plot(rx_x, rx_ant_h, 'rv', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Rx');
plot([tx_x rx_x], [tx_ant_h rx_ant_h], 'g--', 'LineWidth', 1.5, 'DisplayName', 'LOS');

xlabel('Distance (km)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Elevation (m)', 'FontSize', 12, 'FontWeight', 'bold');
title('Terrain Profile & Geometry', 'FontSize', 14);
legend('Location', 'Best');

saveas(fig2, fullfile(Scenario.Output_Dir, 'terrain_geometry.png'));

%% 5. Metrics Analysis (Reference: EFIE)

if Scenario.Run_EFIE
    fprintf('\nCalculating Metrics (Reference: EFIE).\n');
    
    % Metric Functions
    rmse = @(est, ref) sqrt(nanmean((est - ref).^2));
    bias = @(est, ref) nanmean(est - ref); % Positive = Model overestimates loss
    
    % Table Collection
    model_list = {};
    rmse_list = [];
    bias_list = [];
    
    check_model = @(name, data) (~isempty(data) && ~all(isnan(data)));
    
    if DO_ITM && check_model('ITM', Results.ITM_dB)
        model_list{end+1} = 'ITM';
        rmse_list(end+1) = rmse(Results.ITM_dB, Results.EFIE_dB);
        bias_list(end+1) = bias(Results.ITM_dB, Results.EFIE_dB);
    end
    
    if DO_HATA && check_model('Hata', Results.Hata_dB)
        model_list{end+1} = 'Hata';
        rmse_list(end+1) = rmse(Results.Hata_dB, Results.EFIE_dB);
        bias_list(end+1) = bias(Results.Hata_dB, Results.EFIE_dB);
    end
    
    if DO_BULLINGTON && check_model('Bullington', Results.Bullington_dB)
        model_list{end+1} = 'Bullington';
        rmse_list(end+1) = rmse(Results.Bullington_dB, Results.EFIE_dB);
        bias_list(end+1) = bias(Results.Bullington_dB, Results.EFIE_dB);
    end
    
    if DO_DEYGOUT && check_model('Deygout', Results.Deygout_dB)
        model_list{end+1} = 'Deygout';
        rmse_list(end+1) = rmse(Results.Deygout_dB, Results.EFIE_dB);
        bias_list(end+1) = bias(Results.Deygout_dB, Results.EFIE_dB);
    end
    
    % Display Results
    fprintf('\n%-12s | %-10s | %-10s\n', 'Model', 'RMSE (dB)', 'Bias (dB)');
    fprintf('------------------------------------------\n');
    for i = 1:length(model_list)
        fprintf('%-12s | %-10.2f | %-10.2f\n', model_list{i}, rmse_list(i), bias_list(i));
    end
    
    % Save CSV
    T = table(model_list', rmse_list', bias_list', 'VariableNames', {'Model', 'RMSE_dB', 'Bias_dB'});
    writetable(T, fullfile(Scenario.Output_Dir, 'metrics.csv'));
    fprintf('\nMetrics saved to "%s/metrics.csv".\n', Scenario.Output_Dir);
    
else
    fprintf('\nMetrics skipped (EFIE Reference missing).\n');
end

fprintf('\nComparison Complete.\n');