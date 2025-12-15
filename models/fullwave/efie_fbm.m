%% Electric Field Integral Equation (EFIE) with Forward-Backward Method (FBM)
% This script calculates the field scattered from an undulating terrain surface
% which has been excited by an infinite line source.
% 
% Based on Balanis EM Theory (Page 680) and the Forward-Backward Method.
%
% Author: EFIE_FBM Implementation
% Date: December 2024
%
% Outputs:
%   1. Surface current magnitude vs distance along terrain
%   2. Electric field vs distance along terrain
%   3. Data files containing numerical results
%   4. Computation time and memory usage statistics

clear; close all; clc;

%% ========================================================================
%  START TIMING AND MEMORY MEASUREMENT
% =========================================================================
% Record start time for computation (excludes plotting)
computation_start_time = tic;

% Get initial memory usage
try
    % MATLAB memory function
    [user_mem, sys_mem] = memory;
    initial_memory_bytes = user_mem.MemUsedMATLAB;
    memory_tracking_available = true;
catch
    % Octave - use whos to estimate memory
    initial_memory_bytes = 0;
    memory_tracking_available = false;
end

%% ========================================================================
%  PHYSICAL CONSTANTS AND PARAMETERS
% =========================================================================
fprintf('=========================================================\n');
fprintf('  EFIE - Forward-Backward Method for Terrain Scattering\n');
fprintf('=========================================================\n\n');

% Fundamental constants
EXP_E = exp(1);                      % Euler's number e = 2.718281828...
PI = pi;                              % Pi constant

% Electromagnetic constants
Epsilon_0 = 8.854e-12;               % Permittivity of free space (F/m)
Mu_0 = 4*PI*1e-7;                    % Permeability of free space (H/m)
c = 1.0 / sqrt(Mu_0 * Epsilon_0);    % Speed of light in free space (m/s)

% Operating parameters
f = 970e6;                            % Frequency (Hz) - 970 MHz
Lambda = c / f;                       % Wavelength (m)
Omega = 2*PI*f;                       % Angular frequency (rad/s)

% Derived parameters
Beta_0 = Omega * sqrt(Mu_0 * Epsilon_0);  % Free-space wavenumber (rad/m)
Eta_0 = sqrt(Mu_0 / Epsilon_0);           % Free-space intrinsic impedance (Ohms)

fprintf('Operating Frequency: %.3f MHz\n', f/1e6);
fprintf('Wavelength: %.4f m\n', Lambda);
fprintf('Wavenumber (Beta_0): %.4f rad/m\n', Beta_0);
fprintf('Intrinsic Impedance: %.2f Ohms\n\n', Eta_0);

% Discretization parameters
GrossStep = 10.0;                     % Step size in terrain data (m)
NOOFSTEPS = 70;                       % Number of steps to process (1 step = 10m). Set to control distance.

% Choose discretization level:
% Lambda/4 is standard for high accuracy, but Lambda/2 or Lambda gives faster computation
% User can adjust this parameter based on accuracy vs. speed requirements
% Recommended values:
%   4.0 = lambda/4 (highest accuracy, slowest)
%   2.0 = lambda/2 (good accuracy, balanced speed)
%   1.0 = lambda (moderate accuracy, fast)
%   0.5 = 2*lambda (lower accuracy, very fast - for testing)
%
% Correction: Lambda/2 (2.0) is too coarse for pulse basis functions (leads to large phase errors).
% Changed to Lambda/4 (4.0) to ensure physical correctness and reasonable accuracy.
discretization_factor = 4.0;          % Improved accuracy (Lambda/4)
DeltaX = Lambda / discretization_factor;  % Segment length for MoM (m)

% Source parameters
Xsource = 0.0;                        % Source x-coordinate (m)
Ysource = 442.0;                      % Source y-coordinate (height above terrain) (m)
I_source = 1.0;                       % Source current amplitude (A)

% Observation height above terrain surface
ObsHeight = 2.4;                      % Observation height (m)

% Impedance coefficient (used throughout calculations)
Z_coeff = (Beta_0^2) / (4.0 * Omega * Epsilon_0);

%% ========================================================================
%  READ TERRAIN PROFILE DATA FROM X.txt
% =========================================================================
fprintf('Reading terrain data from X.txt...\n');

% Read terrain data
terrain_data = load('X.txt');
X_terrain = terrain_data(:, 1);       % X coordinates (m)
Y_terrain = terrain_data(:, 2);       % Y coordinates (height) (m)

NumTerrainPts = length(X_terrain);
TerrainLength = X_terrain(end);

fprintf('  Number of terrain points: %d\n', NumTerrainPts);
fprintf('  Terrain length: %.1f m\n', TerrainLength);
fprintf('  Terrain height range: %.2f to %.2f m\n', min(Y_terrain), max(Y_terrain));

%% ========================================================================
%  CALCULATE NUMBER OF LINE SEGMENTS (MoM Elements)
% =========================================================================
% Calculate maximum available steps from file data
MaxFileSteps = floor(TerrainLength / GrossStep);

% Apply user limit if NOOFSTEPS is defined and valid
if exist('NOOFSTEPS', 'var') && ~isempty(NOOFSTEPS) && NOOFSTEPS > 0
    if NOOFSTEPS < MaxFileSteps
        GrossNoSteps = NOOFSTEPS;
        TerrainLength = GrossNoSteps * GrossStep; % Update effective terrain length
        fprintf('User Override: Processing limited to %d steps (%.1f m)\n', GrossNoSteps, TerrainLength);
    else
        GrossNoSteps = MaxFileSteps;
        fprintf('User Override: Requested steps (%d) >= available steps (%d). Using full terrain.\n', NOOFSTEPS, MaxFileSteps);
    end
else
    GrossNoSteps = MaxFileSteps;
end

NoLinesubs = floor(TerrainLength / DeltaX);

fprintf('\nDiscretization:\n');
fprintf('  Segment length (DeltaX): %.4f m (lambda/%.1f)\n', DeltaX, discretization_factor);
fprintf('  Number of segments: %d\n', NoLinesubs);

%% ========================================================================
%  PRE-COMPUTE SEGMENT POSITIONS
% =========================================================================
fprintf('\nComputing segment positions...\n');
x_seg = zeros(NoLinesubs, 1);
y_seg = zeros(NoLinesubs, 1);

for idx = 1:NoLinesubs
    % x-coordinate for segment index (MATLAB 1-based indexing)
    x_seg(idx) = (idx-1) * DeltaX;
    
    % y-coordinate by interpolation from terrain data
    x_pos = x_seg(idx);
    Temp = x_pos / GrossStep;
    Index = floor(Temp) + 1;  % MATLAB 1-based indexing
    
    % Boundary check
    if Index < 1
        Index = 1;
    end

    Prop = Temp - floor(Temp);

    % Fixed interpolation boundary logic:
    % If we reach the end of the terrain data, Index can be NumTerrainPts.
    % The original code would clamp Index to NumTerrainPts-1 but Prop was 0 (relative to NumTerrainPts).
    % This would use Y(N-1) instead of Y(N).
    % Correct logic: If we are past N-1, use the last segment for interpolation (N-1 to N).
    if Index >= NumTerrainPts
        Index = NumTerrainPts - 1;
        % Calculate Prop relative to the clamped index (N-1)
        % x_pos is absolute position. (Index-1)*GrossStep is position of Y(Index).
        % We want fraction between Index and Index+1.
        % Prop = (x_pos - x_start_of_segment) / GrossStep
        Prop = (x_pos - (Index-1)*GrossStep) / GrossStep;
    end
    
    y_seg(idx) = Y_terrain(Index) + Prop * (Y_terrain(Index+1) - Y_terrain(Index));
end

fprintf('  Segment x-range: %.2f to %.2f m\n', x_seg(1), x_seg(end));
fprintf('  Segment y-range: %.2f to %.2f m\n', min(y_seg), max(y_seg));

%% ========================================================================
%  PRE-COMPUTE DISTANCES AND SEGMENT LENGTHS
% =========================================================================
fprintf('\nPre-computing distances...\n');

% Distance from source to each segment
R_source_seg = zeros(NoLinesubs, 1);
for p = 1:NoLinesubs
    R_source_seg(p) = sqrt((Xsource - x_seg(p))^2 + (Ysource - y_seg(p))^2);
end

% Distance from source to observation point (2.4m above terrain)
R_source_obs = zeros(NoLinesubs, 1);
for p = 1:NoLinesubs
    R_source_obs(p) = sqrt((Xsource - x_seg(p))^2 + (Ysource - y_seg(p) - ObsHeight)^2);
end

% Segment lengths (distance between consecutive segments)
seg_length = zeros(NoLinesubs, 1);
for p = 1:(NoLinesubs-1)
    seg_length(p) = sqrt((x_seg(p+1) - x_seg(p))^2 + (y_seg(p+1) - y_seg(p))^2);
end
seg_length(NoLinesubs) = DeltaX;  % Last segment uses nominal length

fprintf('  Pre-computation completed.\n');

%% ========================================================================
%  HANKEL FUNCTION H0^(2) AND FIELD CALCULATIONS
% =========================================================================
% H0^(2)(x) = J0(x) - j*Y0(x)
% Inline function for Hankel function of second kind, order zero
calc_H02 = @(Arg) besselj(0, Arg) - 1j * bessely(0, Arg);

% Incident field from line source
calc_EiRad = @(dist) -Z_coeff * calc_H02(Beta_0 * dist);

% Impedance element for given distance R
calc_Z = @(R) Z_coeff * calc_H02(Beta_0 * R);

%% ========================================================================
%  FORWARD-BACKWARD METHOD (FBM) TO SOLVE FOR SURFACE CURRENT
% =========================================================================
fprintf('\n--- Solving for Surface Current using Forward-Backward Method ---\n');

% -------------------------------------------------------------------------
%  SELF-IMPEDANCE CALCULATION (Corrected)
% -------------------------------------------------------------------------
% The original code used a small-argument approximation for Zself which is
% invalid for segments larger than approx lambda/10.
% We now calculate Zself using numerical integration to accurately handle
% the singularity of the Hankel function and the segment size (Lambda/4).
% Formula: Zself = (k*eta/4) * Integral(-L/2 to L/2) of H0^(2)(k*|x|) dx
%                = 2 * Z_coeff * Integral(0 to L/2) of (J0(k*x) - j*Y0(k*x)) dx

fprintf('  Computing self-impedances using numerical integration...\n');

% Define integrands (vectorized for quadgk)
fun_real = @(x) besselj(0, Beta_0 * x);
fun_imag = @(x) bessely(0, Beta_0 * x);

Zself = zeros(NoLinesubs, 1);

for i = 1:NoLinesubs
    L_half = seg_length(i) / 2.0;

    % Real part: Integral of J0 (smooth function)
    % Imag part: Integral of Y0 (logarithmic singularity at 0)
    % quadgk handles endpoint singularities robustly.

    val_real = quadgk(fun_real, 0, L_half);
    val_imag = quadgk(fun_imag, 0, L_half);

    % Combine results
    Zself(i) = 2.0 * Z_coeff * (val_real - 1j * val_imag);
end

% Pre-compute incident fields (vectorized)
Ei = -Z_coeff * (besselj(0, Beta_0 * R_source_seg) - 1j * bessely(0, Beta_0 * R_source_seg));

% Initialize current array
J = zeros(NoLinesubs, 1) + 1j*zeros(NoLinesubs, 1);

% Determine progress update interval based on number of segments
progress_interval = max(100, floor(NoLinesubs / 20));

% -------------------------------------------------------------------------
%  FORWARD SWEEP (Vectorized inner loop)
% -------------------------------------------------------------------------
fprintf('  Forward sweep...\n');

% Initial condition
J(1) = Ei(1) / Zself(1);

for p = 2:NoLinesubs
    % Vectorized distance calculation for all q < p
    q_idx = 1:(p-1);
    dist_pq = sqrt((x_seg(q_idx) - x_seg(p)).^2 + (y_seg(q_idx) - y_seg(p)).^2);
    
    % Vectorized impedance calculation
    Z_pq = Z_coeff * (besselj(0, Beta_0 * dist_pq) - 1j * bessely(0, Beta_0 * dist_pq));
    
    % Vectorized sum
    SUM = sum(seg_length(q_idx) .* Z_pq .* J(q_idx));
    
    J(p) = (Ei(p) - SUM) / Zself(p);
    
    % Progress indicator
    if mod(p, progress_interval) == 0
        fprintf('    Processed %d/%d segments (%.1f%%)\n', p, NoLinesubs, 100*p/NoLinesubs);
    end
end

fprintf('  Forward sweep completed.\n');

% -------------------------------------------------------------------------
%  BACKWARD SWEEP (Vectorized inner loop)
% -------------------------------------------------------------------------
fprintf('  Backward sweep...\n');

for p = (NoLinesubs-1):-1:1
    % Vectorized distance calculation for all q > p
    q_idx = (p+1):NoLinesubs;
    dist_pq = sqrt((x_seg(q_idx) - x_seg(p)).^2 + (y_seg(q_idx) - y_seg(p)).^2);
    
    % Vectorized impedance calculation
    Z_pq = Z_coeff * (besselj(0, Beta_0 * dist_pq) - 1j * bessely(0, Beta_0 * dist_pq));
    
    % Vectorized sum
    SUM = sum(seg_length(q_idx) .* Z_pq .* J(q_idx));
    
    J(p) = J(p) + (-1.0 * SUM) / Zself(p);
    
    % Progress indicator
    if mod(p, progress_interval) == 0 || p == 1
        fprintf('    Processed segment %d (%.1f%%)\n', p, 100*(NoLinesubs-p)/NoLinesubs);
    end
end

fprintf('  Backward sweep completed.\n');

%% ========================================================================
%  CALCULATE TOTAL ELECTRIC FIELD ABOVE THE SURFACE
% =========================================================================
fprintf('\n--- Calculating Total Electric Field ---\n');

% Pre-compute incident fields at observation points
Ei_obs = -Z_coeff * (besselj(0, Beta_0 * R_source_obs) - 1j * bessely(0, Beta_0 * R_source_obs));

Et = zeros(NoLinesubs, 1) + 1j*zeros(NoLinesubs, 1);

for idx = 1:NoLinesubs
    % Vectorized calculation for all segments up to idx
    n_idx = 1:idx;
    
    % Distance from segments to observation point at idx (2.4m above terrain)
    dist_n_idx = sqrt((x_seg(idx) - x_seg(n_idx)).^2 + ((y_seg(idx) + ObsHeight) - y_seg(n_idx)).^2);
    
    % Vectorized impedance calculation
    Z_n = Z_coeff * (besselj(0, Beta_0 * dist_n_idx) - 1j * bessely(0, Beta_0 * dist_n_idx));
    
    % Vectorized scattered field calculation
    E_scattered = sum(J(n_idx) .* seg_length(n_idx) .* Z_n);
    
    % Total field = Incident field - Scattered field
    Et(idx) = Ei_obs(idx) - E_scattered;
    
    % Progress indicator
    if mod(idx, progress_interval) == 0
        fprintf('  Processed %d/%d observation points (%.1f%%)\n', idx, NoLinesubs, 100*idx/NoLinesubs);
    end
end

fprintf('Electric field calculation completed.\n');

%% ========================================================================
%  COMPUTE MAGNITUDES
% =========================================================================
fprintf('\n--- Computing Field Magnitudes ---\n');

% Surface current magnitude
ModJ = abs(J);

% Electric field magnitude (vectorized computation)
ModEt_linear = abs(Et);

% Electric field in dB, normalized by sqrt of distance (vectorized)
valid_idx = (R_source_obs > 0) & (ModEt_linear > 0);
ModEt_dB = -200 * ones(NoLinesubs, 1);  % Initialize with small value
ModEt_dB(valid_idx) = 20.0 * log10(ModEt_linear(valid_idx) ./ sqrt(R_source_obs(valid_idx)));

%% ========================================================================
%  END COMPUTATION TIMING AND MEMORY MEASUREMENT
% =========================================================================
% Stop computation timer (before plotting)
computation_time = toc(computation_start_time);

% Get final memory usage
try
    % MATLAB memory function
    [user_mem, sys_mem] = memory;
    final_memory_bytes = user_mem.MemUsedMATLAB;
    memory_used_bytes = final_memory_bytes - initial_memory_bytes;
catch
    % Octave - estimate memory from workspace variables
    vars_info = whos;
    memory_used_bytes = sum([vars_info.bytes]);
    memory_tracking_available = false;
end

% Convert memory to appropriate units
if memory_used_bytes > 1e9
    memory_used_str = sprintf('%.2f GB', memory_used_bytes / 1e9);
elseif memory_used_bytes > 1e6
    memory_used_str = sprintf('%.2f MB', memory_used_bytes / 1e6);
else
    memory_used_str = sprintf('%.2f KB', memory_used_bytes / 1e3);
end

% Convert time to appropriate units
if computation_time > 3600
    time_str = sprintf('%.2f hours (%.2f seconds)', computation_time / 3600, computation_time);
elseif computation_time > 60
    time_str = sprintf('%.2f minutes (%.2f seconds)', computation_time / 60, computation_time);
else
    time_str = sprintf('%.2f seconds', computation_time);
end

fprintf('\n=========================================================\n');
fprintf('         COMPUTATION PERFORMANCE METRICS                  \n');
fprintf('=========================================================\n');
fprintf('  Computation Time: %s\n', time_str);
fprintf('  Memory Used: %s\n', memory_used_str);
if ~memory_tracking_available
    fprintf('  (Note: Memory estimate based on workspace variables)\n');
end
fprintf('=========================================================\n');

%% ========================================================================
%  SAVE OUTPUT DATA
% =========================================================================
fprintf('\n--- Saving Output Data ---\n');

% Save surface current data
J_output = [x_seg, ModJ];
save('J_current.dat', 'J_output', '-ascii', '-double');
fprintf('  Surface current saved to: J_current.dat\n');

% Save electric field data
E_output = [x_seg, ModEt_dB];
save('E_field.dat', 'E_output', '-ascii', '-double');
fprintf('  Electric field saved to: E_field.dat\n');

% Save linear electric field data
E_linear_output = [x_seg, ModEt_linear];
save('E_field_linear.dat', 'E_linear_output', '-ascii', '-double');
fprintf('  Electric field (linear) saved to: E_field_linear.dat\n');

% Save performance metrics
perf_file = fopen('performance_metrics.txt', 'w');
fprintf(perf_file, 'EFIE-FBM Computation Performance Metrics\n');
fprintf(perf_file, '=========================================\n');
fprintf(perf_file, 'Number of segments: %d\n', NoLinesubs);
fprintf(perf_file, 'Computation Time: %s\n', time_str);
fprintf(perf_file, 'Computation Time (seconds): %.6f\n', computation_time);
fprintf(perf_file, 'Memory Used: %s\n', memory_used_str);
fprintf(perf_file, 'Memory Used (bytes): %.0f\n', memory_used_bytes);
fclose(perf_file);
fprintf('  Performance metrics saved to: performance_metrics.txt\n');

%% ========================================================================
%  GENERATE PLOTS
% =========================================================================
fprintf('\n--- Generating Plots ---\n');

% Set up for headless plotting (no display required)
% This makes the script work in environments without a display (servers, CI)
try
    graphics_toolkit('gnuplot');  % Octave - use gnuplot for file output
catch
    % MATLAB or Octave without gnuplot - will use default
end

% -------------------------------------------------------------------------
%  Plot 1: Terrain Profile
% -------------------------------------------------------------------------
hFig1 = figure('Position', [100, 100, 800, 600], 'Color', 'white', 'Visible', 'off');
plot(X_terrain, Y_terrain, 'b-', 'LineWidth', 1.5);
hold on;
plot(Xsource, Ysource, 'r*', 'MarkerSize', 15, 'LineWidth', 2);
xlabel('Distance along terrain (m)', 'FontSize', 12);
ylabel('Terrain Height (m)', 'FontSize', 12);
title('Terrain Profile with Source Location', 'FontSize', 14, 'FontWeight', 'bold');
legend('Terrain', 'Line Source', 'Location', 'best');
grid on;
xlim([0, TerrainLength]);
print(hFig1, '-dpng', '-r150', 'EFIE_FBM_Terrain_Profile.png');
fprintf('  Figure saved to: EFIE_FBM_Terrain_Profile.png\n');
close(hFig1);

% -------------------------------------------------------------------------
%  Plot 2: Surface Current Magnitude vs Distance
% -------------------------------------------------------------------------
hFig2 = figure('Position', [100, 100, 800, 600], 'Color', 'white', 'Visible', 'off');
plot(x_seg, ModJ, 'r-', 'LineWidth', 1.5);
xlabel('Distance along terrain (m)', 'FontSize', 12);
ylabel('Surface Current Magnitude |J| (A/m)', 'FontSize', 12);
title('Surface Current Distribution on Terrain', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
xlim([0, max(x_seg)]);
print(hFig2, '-dpng', '-r150', 'EFIE_FBM_Surface_Current.png');
fprintf('  Figure saved to: EFIE_FBM_Surface_Current.png\n');
close(hFig2);

% -------------------------------------------------------------------------
%  Plot 3: Electric Field vs Distance (dB scale)
% -------------------------------------------------------------------------
hFig3 = figure('Position', [100, 100, 800, 600], 'Color', 'white', 'Visible', 'off');
plot(x_seg, ModEt_dB, 'g-', 'LineWidth', 1.5);
xlabel('Distance along terrain (m)', 'FontSize', 12);
ylabel('Electric Field (dB rel. to 1/sqrt(R))', 'FontSize', 12);
title('Total Electric Field vs Distance', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
xlim([0, max(x_seg)]);
print(hFig3, '-dpng', '-r150', 'EFIE_FBM_Electric_Field_dB.png');
fprintf('  Figure saved to: EFIE_FBM_Electric_Field_dB.png\n');
close(hFig3);

% -------------------------------------------------------------------------
%  Plot 4: Electric Field vs Distance (Linear scale)
% -------------------------------------------------------------------------
hFig4 = figure('Position', [100, 100, 800, 600], 'Color', 'white', 'Visible', 'off');
semilogy(x_seg, ModEt_linear, 'm-', 'LineWidth', 1.5);
xlabel('Distance along terrain (m)', 'FontSize', 12);
ylabel('Electric Field Magnitude |E| (V/m)', 'FontSize', 12);
title('Total Electric Field vs Distance (Log Scale)', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
xlim([0, max(x_seg)]);
print(hFig4, '-dpng', '-r150', 'EFIE_FBM_Electric_Field_Linear.png');
fprintf('  Figure saved to: EFIE_FBM_Electric_Field_Linear.png\n');
close(hFig4);

% -------------------------------------------------------------------------
%  Plot 5: Surface Current with terrain overlay
% -------------------------------------------------------------------------
hFig5 = figure('Position', [100, 100, 800, 600], 'Color', 'white', 'Visible', 'off');
% Use plotyy for compatibility with Octave
[ax, h1, h2] = plotyy(x_seg, y_seg, x_seg, ModJ);
set(h1, 'LineWidth', 1, 'Color', 'b');
set(h2, 'LineWidth', 1.5, 'Color', 'r');
ylabel(ax(1), 'Terrain Height (m)', 'FontSize', 12);
ylabel(ax(2), 'Surface Current |J| (A/m)', 'FontSize', 12);
xlabel('Distance along terrain (m)', 'FontSize', 12);
title('Surface Current over Terrain', 'FontSize', 14, 'FontWeight', 'bold');
legend([h1; h2], 'Terrain', 'Surface Current', 'Location', 'best');
grid on;
xlim(ax(1), [0, max(x_seg)]);
xlim(ax(2), [0, max(x_seg)]);
print(hFig5, '-dpng', '-r150', 'EFIE_FBM_Overlay_Current.png');
fprintf('  Figure saved to: EFIE_FBM_Overlay_Current.png\n');
close(hFig5);

% -------------------------------------------------------------------------
%  Plot 6: Electric Field with terrain overlay
% -------------------------------------------------------------------------
hFig6 = figure('Position', [100, 100, 800, 600], 'Color', 'white', 'Visible', 'off');
[ax, h1, h2] = plotyy(x_seg, y_seg, x_seg, ModEt_dB);
set(h1, 'LineWidth', 1, 'Color', 'b');
set(h2, 'LineWidth', 1.5, 'Color', 'g');
ylabel(ax(1), 'Terrain Height (m)', 'FontSize', 12);
ylabel(ax(2), 'Electric Field (dB)', 'FontSize', 12);
xlabel('Distance along terrain (m)', 'FontSize', 12);
title('Electric Field over Terrain', 'FontSize', 14, 'FontWeight', 'bold');
legend([h1; h2], 'Terrain', 'Electric Field', 'Location', 'best');
grid on;
xlim(ax(1), [0, max(x_seg)]);
xlim(ax(2), [0, max(x_seg)]);
print(hFig6, '-dpng', '-r150', 'EFIE_FBM_Overlay_Field.png');
fprintf('  Figure saved to: EFIE_FBM_Overlay_Field.png\n');
close(hFig6);

%% ========================================================================
%  SUMMARY STATISTICS
% =========================================================================
fprintf('\n=========================================================\n');
fprintf('                    SUMMARY RESULTS                       \n');
fprintf('=========================================================\n');
fprintf('Surface Current Statistics:\n');
fprintf('  Maximum |J|: %.6e A/m\n', max(ModJ));
fprintf('  Minimum |J|: %.6e A/m\n', min(ModJ));
fprintf('  Mean |J|: %.6e A/m\n', mean(ModJ));
fprintf('\nElectric Field Statistics:\n');
fprintf('  Maximum |E|: %.6e V/m\n', max(ModEt_linear));
fprintf('  Minimum |E|: %.6e V/m\n', min(ModEt_linear));
fprintf('  Mean |E|: %.6e V/m\n', mean(ModEt_linear));
fprintf('  Max E (dB): %.2f dB\n', max(ModEt_dB));
valid_dB = ModEt_dB(ModEt_dB > -200);
if ~isempty(valid_dB)
    fprintf('  Min E (dB): %.2f dB\n', min(valid_dB));
end
fprintf('\nPerformance Metrics (Computation Only - Excludes Plotting):\n');
fprintf('  Computation Time: %s\n', time_str);
fprintf('  Memory Used: %s\n', memory_used_str);
fprintf('\nAnalysis completed successfully!\n');
fprintf('=========================================================\n');
