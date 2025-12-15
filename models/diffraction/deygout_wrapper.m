function PL = deygout_wrapper(Scenario, TerrainProfile)
%DEYGOUT_WRAPPER Calculates Path Loss using the Deygout (3-Edge) Diffraction method.
%
%   PL = deygout_wrapper(Scenario, TerrainProfile)
%
%   Input:
%       Scenario - Struct with Frequency_MHz, Tx_Height_m, Rx_Height_m.
%       TerrainProfile - Struct with Distance_m, Elevation_m.
%
%   Output:
%       PL - Vector of Total Path Loss (dB).

    %% 1. Unpack Inputs
    freq_mhz = Scenario.Frequency_MHz;
    tx_h = Scenario.Tx_Height_m;
    rx_h = Scenario.Rx_Height_m;
    
    dist_m = TerrainProfile.Distance_m(:);
    heights_m = TerrainProfile.Elevation_m(:);
    
    num_points = length(dist_m);
    
    % Wavelength
    c = 299792458;
    lambda = c / (freq_mhz * 1e6);
    
    %% 2. Initialization
    PL = zeros(num_points, 1);
    
    %% 3. Calculation Loop
    for i = 2:num_points
        % Extract profile for current link (Tx -> Current Point)
        curr_dists = dist_m(1:i);
        curr_heights = heights_m(1:i);
        
        % Calculate Deygout Loss
        diff_loss_db = deygout_3_edges(curr_dists, curr_heights, tx_h, rx_h, lambda);
        
        % Free Space Path Loss
        d_km = dist_m(i) / 1000.0;
        if d_km > 0
            fspl_db = 32.44 + 20*log10(freq_mhz) + 20*log10(d_km);
        else
            fspl_db = 0;
        end
        
        PL(i) = fspl_db + diff_loss_db;
    end

end

%% Local Functions

function L_total = deygout_3_edges(dists, heights, tx_h, rx_h, lambda)
    % Calculates diffraction loss using Deygout method (Max 3 edges)

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
    dists_left = dists(1:idx_p);
    heights_left = heights(1:idx_p);

    [nu_t, ~] = find_max_nu(dists_left, heights_left, tx_h, 0, lambda);

    J_t = 0;
    if ~isempty(nu_t) && nu_t > -0.78
        J_t = diffraction_loss_db(nu_t);
    end

    % 3. Right Edge (Principal to Rx)
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

    num_p = length(dists);
    if num_p < 3
        max_nu = [];
        max_idx = -1;
        return;
    end

    R_e = 8500000; % Effective Earth Radius (m) (Standard 4/3 approx)

    d_total = dists(end) - dists(1);

    % Absolute heights of endpoints (Terrain + Antenna)
    h_tx_abs = heights(1) + tx_h;
    h_rx_abs = heights(end) + rx_h;

    max_nu = -inf;
    max_idx = -1;

    for i = 2:num_p-1
        d1 = dists(i) - dists(1);
        d2 = dists(end) - dists(i);

        % LOS height
        h_los = h_tx_abs + (h_rx_abs - h_tx_abs) * (d1 / d_total);

        % Earth bulge
        h_bulge = (d1 * d2) / (2 * R_e);

        % Obstacle clearance parameter h
        h = heights(i) + h_bulge - h_los;

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
    if nu <= -0.78
        J = 0;
    else
        val = sqrt((nu - 0.1)^2 + 1) + nu - 0.1;
        J = 6.9 + 20 * log10(val);
    end
end
