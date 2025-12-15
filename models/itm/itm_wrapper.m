function PL = itm_wrapper(Scenario, TerrainProfile)
%ITM_WRAPPER Calculates Path Loss using the Longley-Rice (ITM) Point-to-Point model.
%
%   PL = itm_wrapper(Scenario, TerrainProfile)
%
%   Input:
%       Scenario - Struct with:
%           Frequency_MHz : Frequency in MHz
%           Tx_Height_m   : Transmitter height (AGL)
%           Rx_Height_m   : Receiver height (AGL)
%           Environment   : (Optional) Climate code implicitly or explicit settings
%
%       TerrainProfile - Struct with:
%           Distance_m    : Vector of distances from Tx (meters)
%           Elevation_m   : Vector of terrain elevations (meters)
%
%   Output:
%       PL - Vector of Total Path Loss (dB).

    %% 1. Unpack Inputs
    freq_mhz = Scenario.Frequency_MHz;
    tx_h = Scenario.Tx_Height_m;
    rx_h = Scenario.Rx_Height_m;
    
    dist_m = TerrainProfile.Distance_m(:);
    elev_m = TerrainProfile.Elevation_m(:);
    
    num_points = length(dist_m);
    
    %% 2. Global Parameters & Configuration
    
    % Determine step size (assumed constant)
    xi = mean(diff(dist_m));
    if isnan(xi) || xi <= 0
        xi = 100; % Default fallback
    end
    
    % Calculate Global Delta-H (Terrain Roughness)
    % This ensures consistency across the path.
    full_pfl.z = elev_m;
    full_pfl.xi = xi;
    
    try
        % Try to use ITM library function if available
        global_delta_h = dlthx(full_pfl, dist_m(1), dist_m(end));
    catch
        % Fallback: Interdecile range (10% - 90%)
        global_delta_h = quantile(elev_m, 0.9) - quantile(elev_m, 0.1);
    end
    
    % ITM Options
    itm_opts = struct();
    itm_opts.pol = 0;           % Horizontal Polarization
    itm_opts.clim = 5;          % Continental Temperate
    itm_opts.conf = 0.5;        % 50% Confidence
    itm_opts.force_delta_h = global_delta_h; % Enforce global roughness
    
    %% 3. Calculation Loop
    PL = zeros(num_points, 1);
    PL(1) = 0; % Loss at Tx is 0/Undefined
    
    % Iterate through each point as a receiver
    for i = 2:num_points
        % Extract sub-profile from Tx to current Rx
        z_sub = elev_m(1:i);
        
        try
            [L_b, ~] = longley_rice_p2p(z_sub, xi, freq_mhz, tx_h, rx_h, itm_opts);
            PL(i) = L_b;
        catch
            PL(i) = NaN;
        end
    end

end
