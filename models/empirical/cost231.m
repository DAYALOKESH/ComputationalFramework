% COST 231 Hata Model Implementation and Plotting
% Implements the COST 231 Hata model for Urban, Suburban, and Rural environments.
% Plots the Path Loss vs Distance.

clear; clc; close all;

tic; % Start timer for performance analysis

%% --- Parameters ---
f = 1500;           % Frequency in MHz
h_B = 52;           % Transmitter (Base Station) height in meters
h_R = 2.4;          % Receiver (Mobile Station) height in meters
d = 1:0.1:20;       % Distance vector from 1 to 20 km with 0.1 km step

%% --- Calculations ---
% Calculate Path Loss for each environment
PL_urban = cost231_hata(f, h_B, h_R, d, 'urban');
PL_suburban = cost231_hata(f, h_B, h_R, d, 'suburban');
PL_rural = cost231_hata(f, h_B, h_R, d, 'rural');

%% --- Plotting ---
figure('Name', 'COST 231 Hata Path Loss', 'NumberTitle', 'off');
plot(d, PL_urban, 'r-', 'LineWidth', 2, 'DisplayName', 'Urban');
hold on;
plot(d, PL_suburban, 'b--', 'LineWidth', 2, 'DisplayName', 'Suburban');
plot(d, PL_rural, 'g-.', 'LineWidth', 2, 'DisplayName', 'Rural');

% Graph formatting
grid on;
xlabel('Distance (km)', 'FontSize', 12);
ylabel('Path Loss (dB)', 'FontSize', 12);
title(sprintf('COST 231 Hata Path Loss\nf = %d MHz, h_B = %d m, h_R = %.1f m', f, h_B, h_R), 'FontSize', 14);
legend('show', 'Location', 'best', 'FontSize', 10);
xlim([1 20]);

%% --- Save Plot ---
output_filename = 'cost231_pathloss.png';
saveas(gcf, output_filename);
fprintf('Plot saved as %s\n', output_filename);

%% --- Performance Analysis ---
elapsed_time = toc;
fprintf('Elapsed Time: %.4f seconds\n', elapsed_time);

% Memory analysis
user_variables = whos;
total_memory_bytes = sum([user_variables.bytes]);
total_memory_kb = total_memory_bytes / 1024;
total_memory_mb = total_memory_kb / 1024;

fprintf('Total Workspace Memory Used: %.2f KB (%.4f MB)\n', total_memory_kb, total_memory_mb);


%% ========================================================================
%  LOCAL FUNCTIONS
% ========================================================================

function PL = cost231_hata(f, h_B, h_R, d, environment)
%COST231_HATA Calculate path loss using COST 231 Hata model
%
%   PL = cost231_hata(f, h_B, h_R, d, environment)
%
%   Input Parameters:
%   ----------------
%   f           : Frequency in MHz (scalar or vector)
%   h_B         : Base station antenna height in meters (scalar or vector)
%   h_R         : Mobile antenna height in meters (scalar or vector)
%   d           : Distance in kilometers (scalar or vector)
%   environment : Type of environment ('urban', 'suburban', 'rural')
%
%   Output:
%   -------
%   PL          : Path loss in dB

    % Input validation
    validateInput(f, h_B, h_R, d, environment);

    % Ensure all inputs are column vectors for consistent broadcasting
    f = f(:);
    h_B = h_B(:);
    h_R = h_R(:);
    d = d(:);

    % Handle Scalar vs Vector broadcasting
    sizes = [numel(f), numel(h_B), numel(h_R), numel(d)];
    max_size = max(sizes);

    if numel(f) == 1, f = repmat(f, max_size, 1); end
    if numel(h_B) == 1, h_B = repmat(h_B, max_size, 1); end
    if numel(h_R) == 1, h_R = repmat(h_R, max_size, 1); end
    if numel(d) == 1, d = repmat(d, max_size, 1); end

    env_lower = lower(string(environment));

    % --- Step 1: Calculate Basic Urban Path Loss ---
    % Determine 'a(h_R)' correction factor based on environment
    % Note: For Suburban and Rural, we use the standard (small/medium city) correction
    % For Urban, COST 231 typically refers to metropolitan areas (Cm=3) or uses specific a(h_R)

    if env_lower == "urban"
        % Urban: Use C_m = 3 dB (Metropolitan) and Urban a(h_R)
        C_m = 3;
        % a(h_R) for Urban (f > 400 MHz)
        a_hR = 3.20 * (log10(11.75 * h_R)).^2 - 4.97;
    else
        % For calculation of Suburban/Rural, we start with the basic Reference Path Loss (Urban Medium/Small)
        % This implies C_m = 0 and the standard a(h_R)
        C_m = 0;
        % a(h_R) for Small/Medium City (used as base for Suburban/Rural)
        a_hR = (1.1 * log10(f) - 0.7) .* h_R - (1.56 * log10(f) - 0.8);
    end

    % Standard COST 231 Hata Formula (Urban Basis)
    PL_base = 46.3 + ...
              33.9 * log10(f) - ...
              13.82 * log10(h_B) - ...
              a_hR + ...
              (44.9 - 6.55 * log10(h_B)) .* log10(d) + ...
              C_m;

    % --- Step 2: Apply Environment Specific Corrections ---

    if env_lower == "urban"
        PL = PL_base;

    elseif env_lower == "suburban"
        % Suburban Correction (relative to Urban Small/Medium City)
        % PL_suburban = PL_urban - 2*[log10(f/28)]^2 - 5.4
        correction = 2 * (log10(f / 28)).^2 + 5.4;
        PL = PL_base - correction;

    elseif env_lower == "rural"
        % Rural / Open Area Correction
        % PL_rural = PL_urban - 4.78*[log10(f)]^2 + 18.33*log10(f) - 40.94
        correction = 4.78 * (log10(f)).^2 - 18.33 * log10(f) + 40.94;
        PL = PL_base - correction;

    else
        error('Unknown environment type.');
    end

    % Ensure output is double
    PL = double(PL);

end

function validateInput(f, h_B, h_R, d, environment)
%VALIDATEINPUT Validate input parameters

    % Check numeric inputs
    if ~isnumeric(f) || any(f <= 0)
        error('Frequency f must be positive numeric value(s) in MHz');
    end

    if ~isnumeric(h_B) || any(h_B <= 0)
        error('Base station height h_B must be positive numeric value(s) in meters');
    end

    if ~isnumeric(h_R) || any(h_R <= 0)
        error('Mobile antenna height h_R must be positive numeric value(s) in meters');
    end

    if ~isnumeric(d) || any(d <= 0)
        error('Distance d must be positive numeric value(s) in kilometers');
    end

    % Check environment
    env_str = string(environment);
    valid_envs = ["urban", "suburban", "rural"];

    if ~any(lower(env_str) == valid_envs)
        error('Environment must be "urban", "suburban", or "rural"');
    end

end
