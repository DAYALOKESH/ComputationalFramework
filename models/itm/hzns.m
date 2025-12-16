function [dL, the] = hzns(pfl, h_g, a_eff)
% HZNS - Horizon Extraction Subroutine for Longley-Rice
%
% Inputs:
%   pfl:   Profile structure
%          pfl.xi - step size (meters)
%          pfl.z  - elevation array (ground heights, meters)
%   h_g:   [tx_height, rx_height] (structural heights, meters)
%   a_eff: Effective Earth Radius (meters)
%
% Outputs:
%   dL:    [dL1, dL2] - Horizon distances (meters)
%   the:   [the1, the2] - Horizon elevation angles (radians)
%
% Algorithm:
%   The horizon search identifies the point along the profile that subtends
%   the maximum elevation angle from each terminal, accounting for Earth
%   curvature through the effective radius.

    np = length(pfl.z);
    xi = pfl.xi;
    z = pfl.z;
    
    % Initialize outputs
    dL = zeros(1, 2);
    the = zeros(1, 2);
    
    % Handle edge case of very short profile (e.g. 2 points)
    if np < 2
        % If only 1 point, distance is 0.
        return;
    end
    
    % --- Transmitter Horizon (Looking Forward) ---
    z_tx_ground = z(1);
    h_tx_total = z_tx_ground + h_g(1);
    
    the_max = -Inf;
    d_L = (np - 1) * xi; % Default to end of profile if no obstruction found
    idx_horizon = np;     % Track index for clarity
    
    % Find point with maximum elevation angle
    for i = 2:np
        dist = (i - 1) * xi; 
        % Elevation angle accounting for Earth curvature
        theta = (z(i) - h_tx_total) / dist - dist / (2 * a_eff);
        
        if theta > the_max
            the_max = theta;
            idx_horizon = i;
            d_L = dist;
        end
    end
    
    dL(1) = d_L;
    the(1) = the_max;
    
    % --- Receiver Horizon (Looking Backward) ---
    z_rx_ground = z(end);
    h_rx_total = z_rx_ground + h_g(2);
    
    the_max = -Inf;
    d_L = (np - 1) * xi;
    idx_horizon = 1;      % Track index for clarity
    
    % Find point with maximum elevation angle (looking back from receiver)
    for i = (np-1):-1:1
        dist = (np - i) * xi; 
        % Elevation angle accounting for Earth curvature
        theta = (z(i) - h_rx_total) / dist - dist / (2 * a_eff);
        
        if theta > the_max
            the_max = theta;
            idx_horizon = i;
            d_L = dist;
        end
    end
    
    dL(2) = d_L;
    the(2) = the_max;

end