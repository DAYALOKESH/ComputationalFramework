%% verify_efie_standalone.m
% Standalone verification script for EFIE FBM Model.
% This script runs the EFIE model exactly as provided in the reference code
% but loads the standard 'data/terrain_profile.txt' file.
% It outputs the 'Electric Field' plots for verification.

clear; close all; clc;

%% 1. Configuration
% Match the Unified Framework settings
f_MHz = 970;
tx_h = 52.0;
rx_h = 2.4;
terrain_file = fullfile('..', '..', 'data', 'terrain_profile.txt');

fprintf('EFIE Verification Run\n');
fprintf('Frequency: %.1f MHz\n', f_MHz);
fprintf('Terrain: %s\n', terrain_file);

%% 2. Load Data
if ~isfile(terrain_file)
    error('Terrain file not found: %s', terrain_file);
end
data = load(terrain_file);
X_full = data(:, 1);
Y_full = data(:, 2);

% --- Truncate to first 700m (User Request) ---
limit_dist = 700.0;
mask = X_full <= limit_dist;
if ~any(mask)
    warning('No points found within %.1fm. Using first 10 points.', limit_dist);
    mask(1:10) = true;
end

X_terrain = X_full(mask);
Y_terrain = Y_full(mask);

NumTerrainPts = length(X_terrain);
TerrainLength = X_terrain(end);

fprintf('Terrain loaded and truncated to: %.2f meters (%d points)\n', TerrainLength, NumTerrainPts);

%% 3. Constants
c = 299792458;
mu_0 = 4*pi*1e-7;
eps_0 = 8.854e-12;
f = f_MHz * 1e6;
lambda = c / f;
omega = 2 * pi * f;
beta_0 = omega * sqrt(mu_0 * eps_0);
Z_coeff = (beta_0^2) / (4.0 * omega * eps_0);

Xsource = X_terrain(1);
Ysource = Y_terrain(1) + tx_h;

%% 4. Discretization (Lambda/4)
discretization_factor = 4.0;
DeltaX = lambda / discretization_factor;
NoLinesubs = floor(TerrainLength / DeltaX);

fprintf('Discretization: %d segments (DeltaX = %.4f m)\n', NoLinesubs, DeltaX);

x_seg = (0:NoLinesubs-1)' * DeltaX;
y_seg = interp1(X_terrain, Y_terrain, x_seg, 'linear', 'extrap');

%% 5. Setup & Impedance
x_next = [x_seg(2:end); x_seg(end) + DeltaX];
y_next = [y_seg(2:end); y_seg(end)];
seg_length = sqrt((x_next - x_seg).^2 + (y_next - y_seg).^2);

R_source_seg = sqrt((Xsource - x_seg).^2 + (Ysource - y_seg).^2);

calc_H02 = @(Arg) besselj(0, Arg) - 1j * bessely(0, Arg);

fprintf('Computing Self Impedance...\n');
fun_real = @(x) besselj(0, beta_0 * x);
fun_imag = @(x) bessely(0, beta_0 * x);

Zself = zeros(NoLinesubs, 1);
for i = 1:NoLinesubs
    L_half = seg_length(i) / 2.0;
    val_real = quadgk(fun_real, 0, L_half);
    val_imag = quadgk(fun_imag, 0, L_half);
    Zself(i) = 2.0 * Z_coeff * (val_real - 1j * val_imag);
end

Ei = -Z_coeff * calc_H02(beta_0 * R_source_seg);

%% 6. Forward-Backward Solution
J = zeros(NoLinesubs, 1);
fprintf('Forward Sweep...\n');
J(1) = Ei(1) / Zself(1);
for p = 2:NoLinesubs
    q_idx = 1:(p-1);
    dist_pq = sqrt((x_seg(q_idx) - x_seg(p)).^2 + (y_seg(q_idx) - y_seg(p)).^2);
    Z_pq = Z_coeff * calc_H02(beta_0 * dist_pq);
    SUM = sum(seg_length(q_idx) .* Z_pq .* J(q_idx));
    J(p) = (Ei(p) - SUM) / Zself(p);
end

fprintf('Backward Sweep...\n');
for p = (NoLinesubs-1):-1:1
    q_idx = (p+1):NoLinesubs;
    dist_pq = sqrt((x_seg(q_idx) - x_seg(p)).^2 + (y_seg(q_idx) - y_seg(p)).^2);
    Z_pq = Z_coeff * calc_H02(beta_0 * dist_pq);
    SUM = sum(seg_length(q_idx) .* Z_pq .* J(q_idx));
    J(p) = J(p) - SUM / Zself(p);
end

%% 7. Calculate Field at Observation Points (Segments)
% The reference code plots E-field vs Distance *along the segments*.
% It computes Et at a height ObsHeight above the segments.

ObsHeight = rx_h; 
R_source_obs = sqrt((Xsource - x_seg).^2 + (Ysource - y_seg - ObsHeight).^2);
Ei_obs = -Z_coeff * calc_H02(beta_0 * R_source_obs);

Et = zeros(NoLinesubs, 1);

fprintf('Calculating Total Field at Observation Points...\n');
% Note: The reference code calculates Et at *every segment x-coordinate*.
% This is computationally heavy for plotting (N*N).
% To match reference exactly, we do it.

for idx = 1:NoLinesubs
    n_idx = 1:idx; % Sum only up to current? No, scattered field is from ALL segments.
    % Wait, reference code loop:
    % for idx = 1:NoLinesubs
    %    n_idx = 1:idx;
    %    ... sum(J(n_idx) ...
    % This implies it only sums contributions from 0 to x.
    % This is physically questionable for backscattering (reflections from ahead),
    % but if that's the reference logic, we follow it?
    % 
    % CHECK REFERENCE CODE PROVIDED:
    % "for idx = 1:NoLinesubs ... n_idx = 1:idx; ... sum(...)"
    % Yes, the reference code only sums scattered field from PREVIOUS segments.
    % This is effectively a "Forward Scattered" approximation for the field calculation loop,
    % even though J was solved with FBM.
    % We will replicate this EXACT behavior for verification.
    
    n_idx = 1:idx; 
    
    dist_n_idx = sqrt((x_seg(idx) - x_seg(n_idx)).^2 + ((y_seg(idx) + ObsHeight) - y_seg(n_idx)).^2);
    Z_n = Z_coeff * calc_H02(beta_0 * dist_n_idx);
    
    E_scattered = sum(J(n_idx) .* seg_length(n_idx) .* Z_n);
    Et(idx) = Ei_obs(idx) - E_scattered;
end

%% 8. Plotting (Match Reference Style)
ModEt_linear = abs(Et);

% Normalized dB: 20*log10( |E| / (1/sqrt(R)) )
% This removes the cylindrical spreading loss to show interference patterns.
ModEt_dB = 20.0 * log10(ModEt_linear ./ (1 ./ sqrt(R_source_obs)));

figure('Name', 'EFIE Verification: Electric Field');
subplot(2,1,1);
plot(x_seg, ModEt_dB, 'g-', 'LineWidth', 1.5);
xlabel('Distance (m)');
ylabel('Normalized E-Field (dB)');
title('Electric Field (Normalized by 1/sqrt(R))');
grid on;

subplot(2,1,2);
semilogy(x_seg, ModEt_linear, 'm-', 'LineWidth', 1.5);
xlabel('Distance (m)');
ylabel('|E| (V/m)');
title('Electric Field Magnitude (Linear)');
grid on;

% Save the figure
output_file = 'verify_efie_plot.png';
saveas(gcf, output_file);
fprintf('Verification Complete. Plot saved to %s\n', output_file);

%% 9. Path Loss Calculation & Plotting (Added Feature)
fprintf('Calculating Path Loss...\n');

% 3D Free Space Path Loss (FSPL) using direct distance (Source to Obs)
% Note: R_source_obs is 2D distance. For 3D FSPL approx, we assume z-axis is infinite line source?
% Actually, for comparability with point-source models (Hata, ITM), we usually
% take the 2D Diffraction Loss (Excess Loss) and add it to 3D FSPL.

% 1. Calculate Excess Loss (Diffraction Loss) from 2D simulation
% L_excess = -20 * log10( |Et| / |Ei| )  (positive dB = loss)
mag_Ei = abs(Ei_obs);
mag_Et = abs(Et);

% Prevent log(0)
mag_Et(mag_Et < 1e-20) = 1e-20;
mag_Ei(mag_Ei < 1e-20) = 1e-20;

L_excess = 20 * log10(mag_Ei ./ mag_Et);

% 2. Calculate 3D FSPL
% Distance in 3D: d_3d = sqrt((x - x_src)^2 + (y - y_src)^2 + (z_rx - z_tx)^2??)
% Here we treat the profile distance R_source_obs as the direct ray length.
fspl_3d = 20 * log10(4 * pi * R_source_obs / lambda);

% 3. Total Path Loss
PL = fspl_3d + L_excess;
PL(1) = 0; % First point

% Plot Path Loss
figure('Name', 'EFIE Verification: Path Loss');
plot(x_seg, PL, 'b-', 'LineWidth', 1.5);
xlabel('Distance (m)');
ylabel('Path Loss (dB)');
title('Total Path Loss (FSPL + EFIE Diffraction)');
grid on;

% Save Path Loss Plot
output_pl_file = 'verify_efie_pathloss.png';
saveas(gcf, output_pl_file);
fprintf('Path Loss Plot saved to %s\n', output_pl_file);
