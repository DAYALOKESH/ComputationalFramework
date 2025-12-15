function delta_h = dlthx(pfl)
% DLTHX - Terrain Irregularity Parameter Calculation
%
% Input:
%   pfl: Profile structure (pfl.z, pfl.xi)
%
% Output:
%   delta_h: Interdecile range of terrain residuals (meters)

    z = pfl.z;
    np = length(z);
    
    % FIX: Handle Short Paths Gracefully
    % If fewer than 10 points, we cannot reliably calculate statistics.
    % Instead of warning/error, return 0 (Smooth Earth assumption)
    % so simulation can proceed for near-field.
    if np < 10
        delta_h = 0;
        return;
    end
    
    % 1. Range Selection
    idx_start = floor(0.1 * np) + 1;
    idx_end = ceil(0.9 * np);
    
    if idx_start < 1, idx_start = 1; end
    if idx_end > np, idx_end = np; end
    
    z_sub = z(idx_start:idx_end);
    n_sub = length(z_sub);
    
    if n_sub < 2
        delta_h = 0;
        return;
    end
    
    x_sub = (0:(n_sub-1))' .* pfl.xi; 
    
    % 2. Trend Removal (Linear Fit)
    X = [ones(n_sub, 1), x_sub];
    coeffs = X \ z_sub(:);
    c = coeffs(1);
    m = coeffs(2);
    
    z_trend = m * x_sub + c;
    
    % 3. Residual Calculation
    residuals = z_sub(:) - z_trend;
    
    % 4. Interdecile Range Calculation
    sorted_res = sort(residuals);
    idx_10 = max(1, round(0.1 * n_sub));
    idx_90 = min(n_sub, round(0.9 * n_sub));
    
    v_10 = sorted_res(idx_10);
    v_90 = sorted_res(idx_90);
    
    delta_h = v_90 - v_10;
    delta_h = max(delta_h, 0);

end