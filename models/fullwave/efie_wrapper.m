function PL = efie_wrapper(Scenario, TerrainProfile)
%EFIE_WRAPPER Calculates Path Loss using EFIE Forward-Backward Method.
%   Matches logic of 'verify_efie_standalone.m' (Reference Implementation).
%
%   PL = efie_wrapper(Scenario, TerrainProfile)
%
%   Input:
%       Scenario - Struct with Frequency_MHz, Tx_Height_m, Rx_Height_m.
%       TerrainProfile - Struct with Distance_m, Elevation_m.
%
%   Output:
%       PL - Vector of Total Path Loss (dB).
%            PL = FSPL(3D) + DiffractionLoss(2D approximation).

    %% 1. Unpack Inputs
    f = Scenario.Frequency_MHz * 1e6;
    tx_h = Scenario.Tx_Height_m;
    rx_h = Scenario.Rx_Height_m;
    
    X_terrain = TerrainProfile.Distance_m(:);
    Y_terrain = TerrainProfile.Elevation_m(:);
    
    num_points = length(X_terrain);
    
    %% 2. Constants & Setup (Reference Implementation)
    c = 299792458; % Use precise c
    mu_0 = 4*pi*1e-7;
    eps_0 = 8.854e-12;
    lambda = c / f;
    omega = 2 * pi * f;
    beta_0 = omega * sqrt(mu_0 * eps_0);
    % eta_0 = sqrt(mu_0 / eps_0); % Not strictly needed for Ratio, but Z_coeff is derived
    
    % Impedance coefficient
    Z_coeff = (beta_0^2) / (4.0 * omega * eps_0);
    
    % Source Location
    Xsource = X_terrain(1);
    Ysource = Y_terrain(1) + tx_h;
    
    %% 3. Discretization (MoM Segments)
    % Use Lambda/4 for accuracy (High Fidelity)
    discretization_factor = 4.0;
    DeltaX = lambda / discretization_factor;
    
    TerrainLength = X_terrain(end);
    NoLinesubs = floor(TerrainLength / DeltaX);
    
    % Status Update
    fprintf('      [EFIE] Running High-Fidelity Simulation with N=%d unknowns...\n', NoLinesubs);
    
    % Pre-compute segment positions
    % We use the same logic as the reference: interpolate Y at uniform X steps
    x_seg = (0:NoLinesubs-1)' * DeltaX;
    
    % Linear interpolation for Y segment heights
    % (Reference code calculates index/prop manually, interp1 is equivalent and faster)
    y_seg = interp1(X_terrain, Y_terrain, x_seg, 'linear', 'extrap');
    
    %% 4. Pre-compute Geometry
    
    % Distance from source to each segment
    R_source_seg = sqrt((Xsource - x_seg).^2 + (Ysource - y_seg).^2);
    
    % Segment lengths (approximate as DeltaX or Euclidean distance)
    % Reference: seg_length(p) = sqrt(dx^2 + dy^2)
    % We need next point for diff.
    x_next = [x_seg(2:end); x_seg(end) + DeltaX];
    y_next = [y_seg(2:end); y_seg(end)]; % Flat extension assumption for last point
    seg_length = sqrt((x_next - x_seg).^2 + (y_next - y_seg).^2);

    %% 5. Impedance Matrix Elements
    
    % H0(2) Inline Function
    calc_H02 = @(Arg) besselj(0, Arg) - 1j * bessely(0, Arg);

    % Self Impedance (Numerical Integration via quadgk - matches Reference)
    fun_real = @(x) besselj(0, beta_0 * x);
    fun_imag = @(x) bessely(0, beta_0 * x);
    
    Zself = zeros(NoLinesubs, 1);
    % Optimization: If seg_length is roughly constant, Zself is roughly constant.
    % But for irregular terrain, it varies. We calculate all.
    for i = 1:NoLinesubs
        L_half = seg_length(i) / 2.0;
        val_real = quadgk(fun_real, 0, L_half);
        val_imag = quadgk(fun_imag, 0, L_half);
        Zself(i) = 2.0 * Z_coeff * (val_real - 1j * val_imag);
    end
    
    % Incident Field at Segments
    Ei = -Z_coeff * calc_H02(beta_0 * R_source_seg);
    
    %% 6. Forward-Backward Method (Solve for Current J)
    J = zeros(NoLinesubs, 1);
    
    % --- Forward Sweep ---
    J(1) = Ei(1) / Zself(1);
    for p = 2:NoLinesubs
        q_idx = 1:(p-1);
        dist_pq = sqrt((x_seg(q_idx) - x_seg(p)).^2 + (y_seg(q_idx) - y_seg(p)).^2);
        Z_pq = Z_coeff * calc_H02(beta_0 * dist_pq);
        SUM = sum(seg_length(q_idx) .* Z_pq .* J(q_idx));
        J(p) = (Ei(p) - SUM) / Zself(p);
    end
    
    % --- Backward Sweep ---
    for p = (NoLinesubs-1):-1:1
        q_idx = (p+1):NoLinesubs;
        dist_pq = sqrt((x_seg(q_idx) - x_seg(p)).^2 + (y_seg(q_idx) - y_seg(p)).^2);
        Z_pq = Z_coeff * calc_H02(beta_0 * dist_pq);
        SUM = sum(seg_length(q_idx) .* Z_pq .* J(q_idx));
        J(p) = J(p) - SUM / Zself(p);
    end
    
    %% 7. Calculate Field at Receiver Points (Original Profile)
    % Calculate Total Field (Et) at the requested Rx points
    
    Et_total = zeros(num_points, 1);
    Ei_total = zeros(num_points, 1);
    
    % Only calculate for points > source (index 2 onwards usually)
    for k = 1:num_points
        x_obs = X_terrain(k);
        y_obs = Y_terrain(k) + rx_h;
        
        dist_src_obs = sqrt((Xsource - x_obs)^2 + (Ysource - y_obs)^2);
        
        if dist_src_obs < 1e-3
            Ei_total(k) = NaN; % Source location
            Et_total(k) = NaN;
            continue;
        end
        
        % Incident Field (Line Source)
        Ei_total(k) = -Z_coeff * calc_H02(beta_0 * dist_src_obs);
        
        % Scattered Field
        dist_seg_obs = sqrt((x_seg - x_obs).^2 + (y_seg - y_obs).^2);
        Z_seg_obs = Z_coeff * calc_H02(beta_0 * dist_seg_obs);
        E_scattered = sum(J .* seg_length .* Z_seg_obs);
        
        Et_total(k) = Ei_total(k) - E_scattered;
    end
    
    %% 8. Convert to Path Loss
    % Strategy: 
    % 1. Calculate Diffraction Loss (Excess Loss) from the 2D EFIE result.
    %    L_diff = -20 * log10( |Et| / |Ei| )
    % 2. Add this to the standard 3D Free Space Path Loss (FSPL).
    %    PL_total = FSPL_3D + L_diff
    
    % 3D FSPL
    dist_3d = sqrt((X_terrain - Xsource).^2 + (Y_terrain + rx_h - (Ysource)).^2);
    fspl_3d = 20 * log10(4 * pi * dist_3d / lambda);
    
    % 2D Diffraction Loss
    mag_Ei = abs(Ei_total);
    mag_Et = abs(Et_total);
    
    % Prevent log(0)
    mag_Et(mag_Et < 1e-20) = 1e-20;
    mag_Ei(mag_Ei < 1e-20) = 1e-20;
    
    L_diff = 20 * log10(mag_Ei ./ mag_Et);
    
    % Combine
    PL = fspl_3d + L_diff;
    
    % Clean up first point
    PL(1) = 0;

end