% EFIE_SOLVER_RUNNER
% Standalone runner for the EFIE Full-Wave Solver.
%
% This script sets up a single scenario and calls the shared implementation
% in 'efie_wrapper.m'.

clc; clear; close all;

%% 1. Configuration
Scenario = struct();
Scenario.Frequency_MHz = 970;
Scenario.Tx_Height_m = 0.0; % Note: Solver often puts Source ON terrain, or above. Wrapper adds this to Y(1).
Scenario.Rx_Height_m = 2.4;

% Input File
if isfile(fullfile('data', 'terrain_profile.txt'))
    filename = fullfile('data', 'terrain_profile.txt');
elseif isfile('X.txt')
    filename = 'X.txt';
else
    error('No terrain file found.');
end

%% 2. Load Data
data = load(filename);
TerrainProfile.Distance_m = data(:, 1);
TerrainProfile.Elevation_m = data(:, 2);

%% 3. Execute Model via Wrapper
fprintf('Running EFIE Model via Wrapper...\n');
tic;
PL_dB = efie_wrapper(Scenario, TerrainProfile);
t = toc;
fprintf('Calculation complete in %.2f seconds.\n', t);

%% 4. Visualization
figure('Name', 'EFIE Path Loss');
plot(TerrainProfile.Distance_m, PL_dB, 'r-', 'LineWidth', 2);
xlabel('Distance (m)');
ylabel('Path Loss (dB)');
title(['EFIE Full-Wave Solution (f = ' num2str(Scenario.Frequency_MHz) ' MHz)']);
grid on;

fprintf('Final Path Loss: %.2f dB\n', PL_dB(end));