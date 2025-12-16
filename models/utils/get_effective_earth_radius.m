function [a_eff, k_factor] = get_effective_earth_radius(N_s)
%GET_EFFECTIVE_EARTH_RADIUS Calculate effective Earth radius
%
% Inputs:
%   N_s: Surface refractivity (N-units), default = 301 for standard atmosphere
%
% Outputs:
%   a_eff: Effective Earth radius (meters)
%   k_factor: Effective Earth radius factor (dimensionless)
%
% Standard atmosphere: N_s = 301, k = 4/3, a_eff = 8.495e6 m
%
% This function ensures geometric consistency across all propagation models
% for diffraction calculations, horizon determinations, and Earth curvature
% corrections.

    if nargin < 1 || isempty(N_s)
        N_s = 301; % Standard atmosphere
    end
    
    % Earth radius
    a_earth = 6371000; % meters
    
    % ITM/ITS standard formula for k-factor
    % k = 1 / (1 - 0.04665 * exp(N_s/179.3))
    k_factor = 1.0 / (1.0 - 0.04665 * exp(N_s / 179.3));
    
    % Effective radius
    a_eff = k_factor * a_earth;
    
    % For N_s = 301: k ≈ 1.3333 (4/3), a_eff ≈ 8.495e6 m
end
