% BULLINGTON_PATHLOSS_RUNNER
% Standalone runner for the Bullington Diffraction Model.
%
% This script sets up a single scenario and calls the shared implementation
% in 'bullington_wrapper.m'.
%
% Usage: Ensure 'data/terrain_profile.txt' or 'X.txt' exists.

clc;
clear;
close all;

%% 1. Configuration
% Modify these parameters to change the simulation
Scenario = struct();
Scenario.Frequency_MHz = 970;
Scenario.Tx_Height_m = 52.0;
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
fprintf('Running Bullington Model...\n');
PL_dB = bullington_wrapper(Scenario, TerrainProfile);

%% 4. Visualization
figure('Name', 'Bullington Path Loss');
plot(TerrainProfile.Distance_m, PL_dB, 'b-', 'LineWidth', 2);
xlabel('Distance (m)');
ylabel('Path Loss (dB)');
title(['Bullington Model (f = ' num2str(Scenario.Frequency_MHz) ' MHz)']);
grid on;

fprintf('Final Path Loss: %.2f dB\n', PL_dB(end));