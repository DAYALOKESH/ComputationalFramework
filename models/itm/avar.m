function [A_total, stats] = avar(A_ref, prop, prop_params)
% AVAR - ITM Variability Calculation
%
% Inputs:
%   A_ref: Median Reference Attenuation (dB)
%   prop: Profile properties (delta_h, etc.)
%   prop_params: climatic and confidence parameters
%
% Outputs:
%   A_total: Attenuation at specified confidence level (dB)
%   stats: Structure with deviations (Y_T, Y_L, Y_S)

    conf = prop_params.conf; % Fraction (e.g., 0.9)
    clim = prop_params.clim_code;
    
    % Standard Normal Deviate z_c
    % We want "Reliability" q. If q=0.9 (90% confident), we want the value
    % that is exceeded 90% of the time? Or we want the loss that is *not* exceeded?
    % Usually "Reliability" implies Service. High reliability -> Low Signal required?
    % ITM: A(q) = A_ref - V_med - Y(q)
    % "Higher reliability implies designing for higher attenuation".
    % So we want the upper tail of the Loss distribution.
    % z_c should be positive for q > 0.5.
    
    z_c = norminv(conf);
    
    % 1. Location Variability (Y_L)
    % Based on terrain irregularity parameter delta_h
    % Formula from ITM specification: sigma_L = 10 * (k*dh) / (k*dh + 13)
    % k is a frequency-dependent factor
    freq = prop_params.freq_mhz;
    
    % Frequency factor k (empirical from ITM measurements)
    % Based on ITM documentation: terrain effects increase with frequency
    % Use logarithmic scaling to capture frequency dependence
    k_f = 1 + log10(freq / 100);  % k=1 at 100MHz, k=2 at 1GHz, k=3.3 at 20GHz
    k_f = max(1, k_f);  % Ensure k >= 1
    
    delta_h = prop.delta_h;
    if delta_h == 0, delta_h = 10; end % Avoid zero, minimum terrain roughness
    
    % Location variability standard deviation
    sigma_L = 10 * (k_f * delta_h) / (k_f * delta_h + 13);
    
    Y_L = sigma_L * z_c;
    
    % 2. Time Variability (Y_T)
    % Climate-dependent based on NBS Technical Note 101 data
    % Varies with climate code (klim) and distance
    clim = prop_params.clim_code;
    
    % Time variability lookup table [climate][distance_range]
    % Distance ranges: <50km, 50-100km, 100-200km, >200km
    % Climate codes: 1=Equatorial, 2=Continental Subtropical, 3=Maritime Subtropical,
    %                4=Desert, 5=Continental Temperate, 6=Maritime Temperate Land,
    %                7=Maritime Temperate Sea
    sigma_T_table = [
        5, 6, 7, 8;    % 1: Equatorial
        4, 5, 6, 7;    % 2: Continental Subtropical
        4, 5, 6, 7;    % 3: Maritime Subtropical
        6, 7, 8, 9;    % 4: Desert (more variable)
        4, 5, 6, 7;    % 5: Continental Temperate
        3, 4, 5, 6;    % 6: Maritime Temperate, Over Land
        2, 3, 4, 5;    % 7: Maritime Temperate, Over Sea (most stable)
    ];
    
    % Determine distance range based on actual path distance if available
    % Distance ranges: <50km, 50-100km, 100-200km, >200km
    if isfield(prop_params, 'dist_km')
        d_km = prop_params.dist_km;
    else
        % Estimate from terrain irregularity (rough approximation)
        d_km = 75;  % Default to middle range
    end
    
    % Select appropriate sigma_T based on distance
    if d_km < 50
        col = 1;
    elseif d_km < 100
        col = 2;
    elseif d_km < 200
        col = 3;
    else
        col = 4;
    end
    
    if clim >= 1 && clim <= 7
        sigma_T = sigma_T_table(clim, col);
    else
        sigma_T = 5;  % Default fallback
    end
    
    Y_T = sigma_T * z_c;
    
    % 3. Situation Variability (Y_S)
    % PDF: sigma_S ~ 5-6 dB.
    sigma_S = 5;
    Y_S = sigma_S * z_c;
    
    % Total Variability Y
    % Combined: Y = sqrt(Y_T^2 + Y_L^2 + Y_S^2)
    % (Assuming independent variables)
    
    Y_total = sqrt(Y_T^2 + Y_L^2 + Y_S^2);
    
    % Note on Sign:
    % A_total = A_ref + Y_total (to get higher loss for high confidence)
    % If conf < 0.5, z_c is negative, Y_total calculated from z_c should handle sign?
    % But we squared them.
    % Correct logic: Y_total = sign(z_c) * sqrt(...)
    
    sign_z = sign(z_c);
    Y_total = sign_z * sqrt( (sigma_T*z_c)^2 + (sigma_L*z_c)^2 + (sigma_S*z_c)^2 );
    
    % Final Attenuation (Basic Transmission Loss)
    % L_b = L_fs + A_ref + Y
    % Wait, A_ref is "Reference Attenuation relative to Free Space".
    % So Total Attenuation relative to Free Space = A_ref + Y.
    
    A_total = A_ref + Y_total;
    
    stats.Y_T = Y_T;
    stats.Y_L = Y_L;
    stats.Y_S = Y_S;
    stats.Y_total = Y_total;

end
