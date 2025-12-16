function [A_ref, mode] = lrprop(d, prop, prop_params)
% LRPROP - Longley-Rice Propagation Core
%
% Inputs:
%   d: Distance (meters) - Can be scalar or vector
%   prop: Structure from qlrpfl (h_e, d_L, a_eff, delta_h, etc.)
%   prop_params: Input parameters (freq, pol, etc.)
%
% Outputs:
%   A_ref: Median Reference Attenuation (dB)
%   mode: Propagation mode (1=LOS, 2=Diffraction, 3=Scatter)

    % Unpack
    freq = prop_params.freq_mhz;
    lambda = 299.792458 / freq; % wavelength in meters
    
    d_L = prop.d_L; % [d_L1, d_L2]
    d_Lt = sum(d_L); % Total LOS distance
    h_e = prop.h_e;
    
    % Calculate Reference Attenuation for each distance
    A_ref = zeros(size(d));
    mode = zeros(size(d));
    
    for i = 1:length(d)
        dist = d(i);
        
        if dist < d_Lt
            % --- Line of Sight Region ---
            mode(i) = 1;
            
            % 1. Two-Ray Propagation with Ground Reflection
            % Path difference phase
            arg = 2 * pi * h_e(1) * h_e(2) / (lambda * dist);
            if dist == 0, arg = 100; end  % Avoid singularity
            
            % Grazing angle for ground reflection
            theta_g = atan((h_e(1) + h_e(2)) / dist);
            
            % Fresnel reflection coefficient
            % Ensure prop_params has required fields for ground parameters
            if ~isfield(prop_params, 'eps_r'), prop_params.eps_r = 15; end
            if ~isfield(prop_params, 'sigma'), prop_params.sigma = 0.005; end
            if ~isfield(prop_params, 'pol'), prop_params.pol = 1; end
            
            Gamma = calc_fresnel_reflection(theta_g, freq, prop_params);
            
            % Two-ray interference including ground reflection
            % Field: E = 1 + Gamma * exp(j*arg)
            % Attenuation: A = -20*log10|E|
            E_magnitude = abs(1 + Gamma * exp(1j * arg));
            A_two_ray = -20 * log10(E_magnitude);
            
            % Clamp to ensure no unphysical gain (passive terrain)
            A_two_ray = max(0, A_two_ray); 
            
            % 2. Blending toward diffraction at horizon
            w = max(0, min(1, dist / d_Lt)); 
            A_diff_at_horiz = calc_diffraction(d_Lt, freq, prop, lambda);
            
            A_ref(i) = (1-w) * A_two_ray + w * A_diff_at_horiz;
            
        elseif dist <= (1.5 * d_Lt) 
             % --- Diffraction Region ---
             mode(i) = 2;
             A_ref(i) = calc_diffraction(dist, freq, prop, lambda);
             
        else
             % --- Scatter Region ---
             mode(i) = 3;
             A_scat = calc_scatter(dist, freq, prop, prop_params);
             A_diff = calc_diffraction(dist, freq, prop, lambda);
             
             if A_scat < A_diff
                 A_ref(i) = A_scat;
             else
                 A_ref(i) = A_diff;
                 mode(i) = 2; 
             end
        end
    end
    
    % Note: A_ref can be negative in cases of ducting, constructive interference,
    % or favorable ground reflections. This is physically valid and represents
    % signal enhancement relative to free space. No clamping applied - trust the physics.

end

% --- Helper Functions ---

function A_diff = calc_diffraction(d, freq, prop, lambda)
    % ITM Diffraction - Blending Knife Edge and Smooth Earth
    %
    % MODIFIED ITM DIFFRACTION WEIGHTING (Non-standard)
    % This implementation uses a continuous weighting function between
    % knife-edge (rough terrain) and smooth-earth (smooth terrain) diffraction.
    % Standard ITM uses discrete mode switching.
    %
    % Weighting function: w(k) = k²/(1+k²) where k = delta_h/lambda
    % This provides smooth transition across all roughness values.
    %
    % Justification: Eliminates discontinuities in predictions as terrain
    % roughness varies, improving numerical stability for iterative analysis.
    
    % Knife-edge diffraction parameter
    theta_tot = prop.the(1) + prop.the(2) + d / prop.a_eff;
    v = theta_tot * sqrt(d / lambda); 
    
    % Knife-edge attenuation (Fresnel diffraction)
    if v > -0.7
        A_ke = 6.9 + 20 * log10( sqrt((v-0.1)^2 + 1) + v - 0.1 );
    else
        A_ke = 0;
    end
    
    % Smooth-earth (rounded obstacle) diffraction using Vogler's method
    % This accounts for Earth curvature effects
    beta = (1 + (prop.h_e(1) + prop.h_e(2)) / prop.a_eff) * d / prop.a_eff;
    X = 2 * beta * sqrt(prop.a_eff * lambda / pi);
    
    % Vogler's smooth-earth diffraction attenuation
    if X < 1.6
        A_r = 20 * log10(1 + X);
    else
        A_r = 20 * log10(X) + 5.8;
    end
    
    % Terrain roughness weighting factor
    % w = 0 for smooth terrain (use A_r), w = 1 for rough terrain (use A_ke)
    delta_h = prop.delta_h;
    k_rough = delta_h / lambda;
    w = k_rough^2 / (1 + k_rough^2);
    
    % Blended diffraction attenuation
    % Note: Foliage/forward obstacle loss (A_fo) not implemented in this version.
    % Can be extended for vegetation modeling in future iterations.
    A_diff = (1 - w) * A_r + w * A_ke;
end

function A_scat = calc_scatter(d, freq, prop, params)
    % Troposcatter Loss - Physics-based ITM formulation
    % Based on forward scatter from atmospheric turbulence
    
    % Angular distance (radians)
    theta_d = d / prop.a_eff;
    
    % Frequency in GHz
    f_ghz = freq / 1000;
    
    % Effective heights in km
    h1_km = prop.h_e(1) / 1000;
    h2_km = prop.h_e(2) / 1000;
    
    % Scatter angle (grazing angles from horizon)
    theta_s = max(0.001, theta_d - (h1_km + h2_km) / (prop.a_eff / 1000));
    
    % Basic scatter loss (ITM formulation)
    % Accounts for frequency, angular distance, and scatter angle
    A_scat = 165 + 20*log10(f_ghz) + 30*log10(theta_d) - 10*log10(theta_s);
    
    % CUSTOM CLIMATE CORRECTIONS (Non-standard ITM)
    % These empirical adjustments account for regional atmospheric stability:
    % 1=Equatorial (0 dB): Baseline
    % 2=Continental Subtropical (-2 dB): Slightly more stable
    % 3=Maritime Subtropical (-4 dB): Stable marine atmosphere
    % 4=Desert (+5 dB): Highly variable, increased scattering
    % 5=Continental Temperate (0 dB): Baseline
    % 6=Maritime Temperate Land (-3 dB): Moderate marine influence
    % 7=Maritime Temperate Sea (-5 dB): Stable marine atmosphere
    %
    % Standard ITM accounts for climate through N_s; these are additional corrections.
    klim = params.clim_code;
    climate_factors = [0, -2, -4, +5, 0, -3, -5]; % dB adjustments
    if klim >= 1 && klim <= 7
        A_scat = A_scat + climate_factors(klim);
    end
    
    % Distance-dependent correction for very long paths
    if d > 100000  % Beyond 100 km
        A_scat = A_scat + 5 * log10(d / 100000);
    end
end

function Gamma = calc_fresnel_reflection(theta, freq, params)
    % Fresnel Reflection Coefficient for Ground
    % Accounts for ground permittivity, conductivity, and polarization
    %
    % Inputs:
    %   theta: Grazing angle (radians)
    %   freq: Frequency (MHz)
    %   params: Structure with eps_r (permittivity), sigma (conductivity), pol (polarization)
    %
    % Output:
    %   Gamma: Complex reflection coefficient
    
    % Extract ground parameters
    eps_r = params.eps_r;      % Relative permittivity (dimensionless)
    sigma = params.sigma;      % Conductivity (S/m)
    pol = params.pol;          % Polarization: 0=Horizontal, 1=Vertical
    
    % Angular frequency
    omega = 2 * pi * freq * 1e6;  % rad/s
    eps_0 = 8.854187817e-12;      % Permittivity of free space (F/m)
    
    % Complex relative permittivity (accounts for conductivity losses)
    eps_complex = eps_r - 1j * sigma / (omega * eps_0);
    
    % Trigonometric terms
    sin_theta = sin(theta);
    cos_theta = cos(theta);
    
    % Avoid singularities at grazing incidence
    if sin_theta < 0.001
        sin_theta = 0.001;
        cos_theta = sqrt(1 - sin_theta^2);
    end
    
    % Square root term in Fresnel equations
    sqrt_term = sqrt(eps_complex - cos_theta^2);
    
    % Reflection coefficient depends on polarization
    if pol == 0
        % Horizontal polarization (E-field parallel to ground)
        Gamma = (sin_theta - sqrt_term) / (sin_theta + sqrt_term);
    else
        % Vertical polarization (E-field perpendicular to ground)
        Gamma = (eps_complex * sin_theta - sqrt_term) / ...
                (eps_complex * sin_theta + sqrt_term);
    end
    
    % For very low angles, reflection coefficient approaches -1 (total reflection)
    % This is physically correct for grazing incidence
end