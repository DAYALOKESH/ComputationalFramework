function PL = bullington_wrapper(Scenario, TerrainProfile)
%BULLINGTON_WRAPPER Calculates Path Loss using the Bullington Knife-Edge Diffraction method.
%
%   PL = bullington_wrapper(Scenario, TerrainProfile)
%
%   Input:
%       Scenario - Struct with:
%           Frequency_MHz : Frequency in MHz
%           Tx_Height_m   : Transmitter height (AGL)
%           Rx_Height_m   : Receiver height (AGL)
%
%       TerrainProfile - Struct with:
%           Distance_m    : Vector of distances from Tx (meters)
%           Elevation_m   : Vector of terrain elevations (meters)
%
%   Output:
%       PL - Vector of Total Path Loss (dB) corresponding to TerrainProfile.Distance_m.
%            PL includes Free Space Path Loss + Diffraction Loss.

    %% 1. Unpack Inputs
    freq = Scenario.Frequency_MHz * 1e6; % Convert to Hz
    tx_h = Scenario.Tx_Height_m;
    rx_h = Scenario.Rx_Height_m;
    
    d_vec = TerrainProfile.Distance_m(:);
    z_vec = TerrainProfile.Elevation_m(:);
    
    num_points = length(d_vec);
    
    %% 2. Constants
    c = 3e8;
    lambda = c / freq;
    % Use standardized effective Earth radius for consistency across models
    [R_effective, ~] = get_effective_earth_radius(301);
    
    %% 3. Initialization
    PL = zeros(num_points, 1);
    
    % Tx Location (Index 1)
    d_tx = d_vec(1);
    z_tx = z_vec(1);
    h_tx_ant = z_tx + tx_h;
    
    PL(1) = 0; % Loss at Tx is 0 (or undefined)
    
    %% 4. Calculation Loop
    % Calculate PL for every point in the profile assuming it is the Receiver
    
    for i = 2:num_points
        % Rx Location (Current Point)
        d_rx = d_vec(i);
        z_rx = z_vec(i);
        h_rx_ant = z_rx + rx_h;
        
        D_total = d_rx - d_tx;
        
        % --- Free Space Path Loss ---
        fspl = 20 * log10(4 * pi * D_total / lambda);
        
        % --- Bullington Diffraction ---
        % Segment from Tx to current Rx
        x_segment = d_vec(1:i) - d_tx;      % Relative distance 0 to D_total
        z_segment = z_vec(1:i);
        
        % Earth Curvature Correction (Bulge)
        % h_eff = h + (d1 * d2) / (2 * Re)
        d1 = x_segment;
        d2 = D_total - x_segment;
        bulge = (d1 .* d2) / (2 * R_effective);
        h_eff = z_segment + bulge;
        
        L_diff = 0;
        
        if length(x_segment) > 2
            % Intermediate points only
            x_mid = x_segment(2:end-1);
            h_mid = h_eff(2:end-1);
            
            % 1. Slopes from Tx
            % m = (h - h_tx) / x
            slopes_tx = (h_mid - h_tx_ant) ./ x_mid;
            [m_tx, ~] = max(slopes_tx);
            
            % 2. Slopes from Rx (looking back)
            % m' = (h - h_rx) / (D - x)
            dist_from_rx = D_total - x_mid;
            slopes_rx = (h_mid - h_rx_ant) ./ dist_from_rx;
            [m_rx_prime, ~] = max(slopes_rx);
            
            % 3. Equivalent Obstacle Intersection
            % h_eq = h_tx + m_tx * x_eq
            % h_eq = h_rx + m_rx_prime * (D - x_eq)
            
            denominator = m_tx + m_rx_prime;
            
            if abs(denominator) > 1e-10
                x_eq = (h_rx_ant - h_tx_ant + m_rx_prime * D_total) / denominator;
                h_eq = h_tx_ant + m_tx * x_eq;
                
                % 4. LOS Height at x_eq
                % Line from Tx to Rx (straight in effective earth)
                slope_los = (h_rx_ant - h_tx_ant) / D_total;
                h_los_eq = h_tx_ant + slope_los * x_eq;
                
                % 5. Clearance height (h_obs)
                h_obs = h_eq - h_los_eq;
                
                % 6. Fresnel Parameter nu
                % nu = h_obs * sqrt(2(d1+d2)/(lambda*d1*d2))
                d1_eq = x_eq;
                d2_eq = D_total - x_eq;
                
                % Check validity
                if d1_eq > 0 && d2_eq > 0
                     nu = h_obs * sqrt(2 * (d1_eq + d2_eq) / (lambda * d1_eq * d2_eq));
                     
                     % 7. Diffraction Loss (dB) - ITU-R P.526
                     if nu > -0.78
                         J_nu = 6.9 + 20 * log10(sqrt((nu - 0.1)^2 + 1) + nu - 0.1);
                         L_diff = J_nu;
                     end
                end
            end
        end
        
        PL(i) = fspl + L_diff;
    end

end
