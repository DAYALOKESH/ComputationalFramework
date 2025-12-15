function deygout_simulation()
    % DEYGOUT_SIMULATION Calculates Total Path Loss (FSPL + Deygout Diffraction)
    % for a moving receiver along a terrain profile.
    %
    % Configuration:
    %   Frequency: 970 MHz
    %   Tx Height: 52 m
    %   Rx Height: 2.4 m
    %   Max Edges: 3 (Principal + Left + Right)
    %   Input: X.txt (Distance [m], Height [m])

    % --- Configuration ---
    filename = 'X.txt';
    freq_mhz = 970.0;
    tx_h = 52.0; % Transmitter Height (m)
    rx_h = 2.4;  % Receiver Height (m)

    % --- Load Data ---
    if ~exist(filename, 'file')
        error('Error: File %s not found in the current directory.', filename);
    end
    data = load(filename);

    % Data format: Column 1 = Distance (m), Column 2 = Height (m)
    dist_m = data(:, 1);
    heights_m = data(:, 2);

    num_points = length(dist_m);

    % Pre-allocate results array
    total_path_loss = zeros(num_points, 1);

    % Wavelength Calculation
    c = 299792458;
    lambda = c / (freq_mhz * 1e6);

    fprintf('Starting Deygout Diffraction Simulation...\n');
    fprintf('  Points: %d\n', num_points);
    fprintf('  Frequency: %.1f MHz\n', freq_mhz);
    fprintf('  Tx Height: %.1f m, Rx Height: %.1f m\n', tx_h, rx_h);

    % --- Execution Loop with Metrics ---
    t_start = tic;

    % Loop through every point treating it as a receiver location
    % We start from index 2 because index 1 is the Transmitter location (Dist 0)
    for i = 2:num_points
        % Extract profile for the current link (Tx -> Current Point)
        curr_dists = dist_m(1:i);
        curr_heights = heights_m(1:i);

        % Calculate Deygout Diffraction Loss (Strict 3-edge limit)
        diff_loss_db = deygout_3_edges(curr_dists, curr_heights, tx_h, rx_h, lambda);

        % Calculate Free Space Path Loss (FSPL)
        d_km = dist_m(i) / 1000.0;
        if d_km > 0
            fspl_db = 32.44 + 20*log10(freq_mhz) + 20*log10(d_km);
        else
            fspl_db = 0;
        end

        % Total Loss
        total_path_loss(i) = fspl_db + diff_loss_db;
    end

    calc_time = toc(t_start);

    % --- Performance Reporting ---
    fprintf('\n--------------------------------------------------\n');
    fprintf('Performance Metrics (Core Computation):\n');
    fprintf('  Execution Time: %.6f seconds\n', calc_time);

    % Memory estimation (Approximate size of main data arrays)
    % Variables: dist_m, heights_m, total_path_loss (Doubles = 8 bytes)
    % Note: Intermediate variables in loop are transient.
    mem_usage_bytes = (length(dist_m)*8) * 2 + (length(total_path_loss)*8);
    fprintf('  Approx. Data Memory Usage: %.2f KB\n', mem_usage_bytes/1024);
    fprintf('--------------------------------------------------\n');

    % --- Plotting ---
    % Create a figure (invisible if running in non-GUI env, but standard commands used)
    h_fig = figure('Visible', 'off');

    % Subplot 1: Terrain Profile
    subplot(2, 1, 1);
    plot(dist_m, heights_m, 'k-', 'LineWidth', 1.5);
    hold on;
    % Draw Tx mast
    plot([dist_m(1), dist_m(1)], [heights_m(1), heights_m(1)+tx_h], 'r-', 'LineWidth', 2);
    % Draw example Rx mast at the end
    plot([dist_m(end), dist_m(end)], [heights_m(end), heights_m(end)+rx_h], 'b-', 'LineWidth', 2);
    grid on;
    title('Terrain Profile');
    xlabel('Distance (m)');
    ylabel('Height (m)');
    legend('Terrain', 'Tx Location', 'Rx Location (End)');

    % Subplot 2: Total Path Loss
    subplot(2, 1, 2);
    % Plot from index 2 to avoid the -Inf FSPL at distance 0
    plot(dist_m(2:end), total_path_loss(2:end), 'b-', 'LineWidth', 1.5);
    grid on;
    title('Total Pathloss vs Distance (FSPL + Deygout)');
    xlabel('Distance (m)');
    ylabel('Path Loss (dB)');

    % Save plot
    output_plot_file = 'simulation_results.png';
    print(h_fig, output_plot_file, '-dpng');
    fprintf('Results plotted and saved to %s\n', output_plot_file);
end

function L_total = deygout_3_edges(dists, heights, tx_h, rx_h, lambda)
    % Calculates diffraction loss using Deygout method (Max 3 edges)
    % Inputs:
    %   dists: vector of distances (m)
    %   heights: vector of terrain heights (m)
    %   tx_h: height of start antenna (relative to terrain)
    %   rx_h: height of end antenna (relative to terrain)
    %   lambda: wavelength (m)

    % 1. Principal Edge Identification
    [nu_p, idx_p] = find_max_nu(dists, heights, tx_h, rx_h, lambda);

    % If no principal edge or below threshold, return 0
    if isempty(nu_p) || nu_p <= -0.78
        L_total = 0;
        return;
    end

    % Principal Edge Loss
    J_p = diffraction_loss_db(nu_p);

    % 2. Left Edge (Tx to Principal)
    % Path: Start -> Principal Edge
    % The 'Receiver' at the Principal Edge is the obstacle tip itself.
    % Thus, effective receiver height relative to the obstacle height is 0.
    dists_left = dists(1:idx_p);
    heights_left = heights(1:idx_p);

    [nu_t, ~] = find_max_nu(dists_left, heights_left, tx_h, 0, lambda);

    J_t = 0;
    if ~isempty(nu_t) && nu_t > -0.78
        J_t = diffraction_loss_db(nu_t);
    end

    % 3. Right Edge (Principal to Rx)
    % Path: Principal Edge -> End
    % The 'Transmitter' at the Principal Edge is the obstacle tip.
    % Thus, effective transmitter height relative to the obstacle height is 0.
    dists_right = dists(idx_p:end);
    heights_right = heights(idx_p:end);

    [nu_r, ~] = find_max_nu(dists_right, heights_right, 0, rx_h, lambda);

    J_r = 0;
    if ~isempty(nu_r) && nu_r > -0.78
        J_r = diffraction_loss_db(nu_r);
    end

    % Apply Corrections
    if J_p <= 6
        T = J_p / 6.0;
    else
        T = 1.0;
    end

    total_dist_km = (dists(end) - dists(1)) / 1000.0;
    C = 8.0 + 0.04 * total_dist_km;

    % Final Deygout Loss
    L_total = J_p + T * (J_t + J_r + C);
end

function [max_nu, max_idx] = find_max_nu(dists, heights, tx_h, rx_h, lambda)
    % Finds the maximum Fresnel parameter nu along the path
    % Returns:
    %   max_nu: The maximum Fresnel parameter found
    %   max_idx: The absolute index in 'dists' where max_nu occurs

    num_p = length(dists);
    % Need at least 3 points (Start, Obstacle, End)
    if num_p < 3
        max_nu = [];
        max_idx = -1;
        return;
    end

    R_e = 8500000; % Effective Earth Radius (m) (Standard 4/3)

    d_total = dists(end) - dists(1);

    % Absolute heights of endpoints (Terrain + Antenna)
    h_tx_abs = heights(1) + tx_h;
    h_rx_abs = heights(end) + rx_h;

    max_nu = -inf;
    max_idx = -1;

    % Iterate over intermediate points only
    for i = 2:num_p-1
        d1 = dists(i) - dists(1);
        d2 = dists(end) - dists(i);

        % Line of Sight (LOS) height at this distance
        h_los = h_tx_abs + (h_rx_abs - h_tx_abs) * (d1 / d_total);

        % Earth bulge height
        h_bulge = (d1 * d2) / (2 * R_e);

        % Obstacle clearance parameter h
        % h = (Terrain + Bulge) - LOS
        h = heights(i) + h_bulge - h_los;

        % Fresnel Parameter nu
        % Avoid division by zero (should cover by i range, but for safety)
        if d1 > 0 && d2 > 0
            nu = h * sqrt(2 * (d1 + d2) / (lambda * d1 * d2));

            if nu > max_nu
                max_nu = nu;
                max_idx = i;
            end
        end
    end
end

function J = diffraction_loss_db(nu)
    % Calculates Knife-Edge Diffraction Loss (dB)
    % J(v) = 6.9 + 20log10(sqrt((v-0.1)^2 + 1) + v - 0.1)

    if nu <= -0.78
        J = 0;
    else
        val = sqrt((nu - 0.1)^2 + 1) + nu - 0.1;
        J = 6.9 + 20 * log10(val);
    end
end
