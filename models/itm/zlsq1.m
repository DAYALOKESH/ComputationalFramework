function z_eff = zlsq1(z, xi, idx_start, idx_end, idx_target)
% ZLSQ1 - Least Squares Fit to Profile Segment
%
% Inputs:
%   z: Profile elevation array
%   xi: Step size
%   idx_start, idx_end: Indices defining the segment to fit (1-based)
%   idx_target: Index where we want to evaluate the fitted line (usually 1 or end)
%
% Output:
%   z_eff: Elevation of the fitted line at idx_target

    % Extract segment
    % Ensure indices are valid
    n = length(z);
    idx_start = max(1, min(n, idx_start));
    idx_end = max(1, min(n, idx_end));
    
    if idx_start == idx_end
        z_eff = z(idx_start);
        return;
    end

    z_seg = z(idx_start:idx_end);
    n_seg = length(z_seg);
    
    % X coordinates relative to the START of the whole profile (to align with target)
    x_seg = ((idx_start:idx_end)' - 1) * xi;
    
    % Fit Line L(x) = mx + c
    X = [ones(n_seg, 1), x_seg];
    coeffs = X \ z_seg(:);
    c = coeffs(1);
    m = coeffs(2);
    
    % Evaluate at target x
    x_target = (idx_target - 1) * xi;
    z_eff = m * x_target + c;

end
