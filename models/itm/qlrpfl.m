function prop = qlrpfl(pfl, prop_params)
% QLRPFL - Quick Longley-Rice Profile Analysis
%
% Inputs:
%   pfl: Profile structure (z, xi)
%   prop_params: Structure with h_g, N_s, etc.
%
% Output:
%   prop: Updated structure containing:
%         prop.d_L (Horizons)
%         prop.the (Horizon Angles)
%         prop.delta_h (Terrain Roughness)
%         prop.h_e (Effective Heights)
%         prop.a_eff (Effective Earth Radius)

    prop = struct();

    % 1. Physical Constants
    N_s = prop_params.N_s;
    % Calculate Effective Earth Radius
    % gamma_e = gamma_a * (1 - 0.04665 * exp(Ns/N1)) where N1=179.3 usually.
    % Simplified common ITM approx: k = 1 / (1 - 0.04665 * exp(Ns/179.3))
    % Standard ITM k-factor formula used in ITS code:
    k_factor = 1.0 / (1.0 - 0.04665 * exp(N_s / 179.3));
    a_earth = 6371000; % meters
    prop.a_eff = k_factor * a_earth;
    
    % 2. Horizons
    [prop.d_L, prop.the] = hzns(pfl, prop_params.h_g, prop.a_eff);
    
    % 3. Delta H
    prop.delta_h = dlthx(pfl);
    
    % 4. Effective Heights (h_e)
    % Logic: Fit smooth curve to "foreground"
    % ITM Standard: Foreground is 0 to d_L, but clamped to not exceed certain distance?
    % The 'Longley-Rice Model Implementation Guide' implies fitting from antenna to horizon.
    
    xi = pfl.xi;
    dist_total = (length(pfl.z) - 1) * xi;
    
    % --- Tx Effective Height ---
    d_L1 = prop.d_L(1);
    % Limit fit range (e.g., don't go beyond horizon)
    idx_L1 = floor(d_L1 / xi) + 1;
    % Clamp range (minimum 2 points, max full profile)
    idx_end_tx = max(2, min(length(pfl.z), idx_L1));
    
    z_eff_tx = zlsq1(pfl.z, xi, 1, idx_end_tx, 1);
    
    % h_e = h_g + z_ground - z_eff
    % (If flat ground, z_ground = z_eff, so h_e = h_g)
    % (If cliff, z_eff < z_ground, h_e increases)
    prop.h_e(1) = prop_params.h_g(1) + pfl.z(1) - z_eff_tx;
    
    % --- Rx Effective Height ---
    d_L2 = prop.d_L(2);
    % Distance from Rx back to its horizon
    idx_L2_from_end = floor(d_L2 / xi); 
    idx_start_rx = max(1, length(pfl.z) - idx_L2_from_end);
    idx_start_rx = min(length(pfl.z) - 1, idx_start_rx); % Ensure at least 1 point before end
    
    z_eff_rx = zlsq1(pfl.z, xi, idx_start_rx, length(pfl.z), length(pfl.z));
    
    prop.h_e(2) = prop_params.h_g(2) + pfl.z(end) - z_eff_rx;
    
    % Safety: h_e must not be negative or too small (physically impossible for propagation)
    % ITM often enforces a minimum h_e? 
    % We will clamp to a small positive value if calculation goes haywire.
    prop.h_e(1) = max(prop.h_e(1), 1.0);
    prop.h_e(2) = max(prop.h_e(2), 1.0);

end
