function [loss_db, details] = longley_rice_p2p(profile_z, profile_xi, freq_mhz, h_tx, h_rx, options)
% LONGLEY_RICE_P2P - Point-to-Point Longley-Rice Model
%
% Inputs:
%   profile_z:  Terrain elevation array (meters)
%   profile_xi: Step size (meters)
%   freq_mhz:   Frequency (MHz)
%   h_tx, h_rx: Structural heights (meters)
%   options:    (Optional) Structure with fields:
%               .pol (0=Horiz, 1=Vert, Default 0)
%               .eps_r (Dielectric constant, Default 15)
%               .sigma (Conductivity, Default 0.005)
%               .N_s (Refractivity, Default 301)
%               .clim (Climate code 1-7, Default 5)
%               .conf (Confidence 0-1, Default 0.5)
%
% Outputs:
%   loss_db:    Basic Transmission Loss (dB)
%   details:    Structure with intermediate values (d_L, h_e, A_ref, etc.)

    % --- 1. Validation & Setup ---
    if nargin < 6, options = struct(); end
    
    % Defaults
    if ~isfield(options, 'pol'), options.pol = 0; end
    if ~isfield(options, 'eps_r'), options.eps_r = 15; end
    if ~isfield(options, 'sigma'), options.sigma = 0.005; end
    if ~isfield(options, 'N_s'), options.N_s = 301; end
    if ~isfield(options, 'clim'), options.clim = 5; end
    if ~isfield(options, 'conf'), options.conf = 0.5; end
    
    % Strict Input Ranges (Kimi's Requirement)
    if freq_mhz < 20 || freq_mhz > 20000
        error('Frequency must be between 20 and 20,000 MHz');
    end
    if options.eps_r < 1 || options.eps_r > 100
        error('Epsilon_r must be between 1 and 100');
    end
    if length(profile_z) < 10
        error('Profile must have at least 10 points');
    end
    
    % Build Internal Structs
    prop_params = options;
    prop_params.freq_mhz = freq_mhz;
    prop_params.h_g = [h_tx, h_rx];
    prop_params.clim_code = options.clim;
    
    pfl.z = profile_z;
    pfl.xi = profile_xi;
    
    % --- 2. Geometry (Preparatory Subroutines) ---
    prop = qlrpfl(pfl, prop_params);

    % Override Delta-H if provided to ensure consistency along a path
    if isfield(options, 'force_delta_h')
        prop.delta_h = options.force_delta_h;
    end
    
    % --- 3. Propagation Physics ---
    % Calculate for the full path distance
    dist_total = (length(pfl.z) - 1) * pfl.xi;
    
    [A_ref, mode] = lrprop(dist_total, prop, prop_params);
    
    % --- 4. Variability ---
    % Add distance information for climate variability lookup
    prop_params.dist_km = dist_total / 1000;
    [A_var, var_stats] = avar(A_ref, prop, prop_params);
    
    % --- 5. Final Calculation ---
    % Basic Transmission Loss L_b = L_bf (Free Space) + A_total
    lambda = 299.792458 / freq_mhz;
    L_bf = 20 * log10(4 * pi * dist_total / lambda);
    
    loss_db = L_bf + A_var;
    
    % Pack Details
    details.L_bf = L_bf;
    details.A_ref = A_ref;
    details.A_var = A_var;
    details.mode = mode;
    details.prop = prop;
    details.stats = var_stats;
    details.dist_total = dist_total;

end
