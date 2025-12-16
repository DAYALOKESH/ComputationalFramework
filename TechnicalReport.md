# Technical Report: Unified Radio Propagation Computational Framework

## PART 1: TECHNICAL REPORT

---

## 1. System Architecture Overview

### 1.1 Framework Organization

The computational framework implements a hybrid architecture with distinct separation between terrain data acquisition and propagation modeling. The system comprises two primary components:

1. **Python-based Terrain Extraction Module** (`api/dem_profile_api.py`): Interfaces with Google Earth Engine API to retrieve elevation profiles along geodesic paths
2. **MATLAB Computational Engine** (`api/main_propagation_comparison.m`): Executes propagation algorithms and performs comparative analysis

Data exchange occurs through file-based communication using standardized text files containing distance-elevation pairs. This loose coupling enables independent development, testing, and version control of each component.

### 1.2 Directory Structure

```
ComputationalFramework/
├── api/
│   ├── dem_profile_api.py          # Python GEE interface
│   └── main_propagation_comparison.m # Master harness
├── models/
│   ├── itm/                        # Longley-Rice implementation
│   ├── empirical/                  # Hata/COST-231 models
│   ├── diffraction/                # Bullington and Deygout
│   ├── fullwave/                   # EFIE solver
│   └── utils/                      # Shared utilities
├── data/
│   └── terrain_profile.txt         # Terrain elevation data
├── output/
│   ├── comparison_plot.png         # Generated visualizations
│   ├── terrain_geometry.png
│   └── metrics.csv                 # Statistical results
└── setup_project.m                 # MATLAB path configuration
```

### 1.3 Execution Flow

The master harness (`main_propagation_comparison.m`) orchestrates execution:

1. Loads terrain profile from `data/terrain_profile.txt`
2. Validates input data and removes NaN values
3. Executes each propagation model through standardized wrapper functions
4. Generates comparison plots
5. Computes statistical metrics (RMSE, Bias) using EFIE as reference
6. Saves results to CSV file

---

## 2. Terrain Data Acquisition Pipeline

### 2.1 Python Implementation Details

The terrain extraction module (`dem_profile_api.py`) implements a class-based architecture using dataclasses for structured data representation.

#### 2.1.1 DEM Dataset Enumeration

The code defines available datasets through an enumeration (lines 25-31):

```python
class DEMDataset(Enum):
    SRTM_30M = "USGS/SRTMGL1_003"
    SRTM_90M = "CGIAR/SRTM90_V4"
    ALOS_30M = "JAXA/ALOS/AW3D30/V3_2"
    COPERNICUS_30M = "COPERNICUS/DEM/GLO30"
    NASADEM = "NASA/NASADEM_HGT/001"
```

The default dataset is `SRTM_30M` (USGS SRTMGL1_003).

#### 2.1.2 Great-Circle Path Interpolation

The `_interpolate_great_circle` method (lines 205-246) generates uniformly spaced coordinates along a geodesic path using spherical interpolation:

```python
# Angular distance calculation
d = 2 * math.asin(math.sqrt(
    math.sin((lat2 - lat1) / 2) ** 2 +
    math.cos(lat1) * math.cos(lat2) * math.sin((lon2 - lon1) / 2) ** 2
))

# Spherical linear interpolation
A = math.sin((1 - fraction) * d) / math.sin(d)
B = math.sin(fraction * d) / math.sin(d)

x = A * cos(lat1) * cos(lon1) + B * cos(lat2) * cos(lon2)
y = A * cos(lat1) * sin(lon1) + B * cos(lat2) * sin(lon2)
z = A * sin(lat1) + B * sin(lat2)
```

#### 2.1.3 Haversine Distance Calculation

Distance computation uses the haversine formula (lines 190-203):

```python
R = 6371000  # Earth radius in meters
a = (sin(delta_lat / 2) ** 2 + 
     cos(lat1_rad) * cos(lat2_rad) * sin(delta_lon / 2) ** 2)
c = 2 * atan2(sqrt(a), sqrt(1 - a))
distance = R * c
```

#### 2.1.4 Batch Elevation Retrieval

The `_fetch_elevations_batch` method (lines 259-294) retrieves elevation data in batches of 100 coordinates using Google Earth Engine's `sampleRegions` method:

```python
sampled = dem.select(band).sampleRegions(
    collection=points,
    scale=30,
    geometries=False
)
```

#### 2.1.5 Output Format

The module outputs distance-elevation pairs in ASCII format:
- Column 1: Cumulative distance in meters (starting at 0)
- Column 2: Terrain elevation in meters above WGS84 ellipsoid
- Space-delimited, one pair per line

---

## 3. Propagation Model Implementations

### 3.1 Electric Field Integral Equation (EFIE) Reference Model

#### 3.1.1 Physical Constants and Discretization

The EFIE implementation (`efie_wrapper.m`, lines 16-44) uses:

```matlab
c = 299792458;           % Speed of light (m/s)
mu_0 = 4*pi*1e-7;        % Permeability of free space
eps_0 = 8.854e-12;       % Permittivity of free space
lambda = c / f;          % Wavelength
omega = 2 * pi * f;      % Angular frequency
beta_0 = omega * sqrt(mu_0 * eps_0);  % Free-space wavenumber

% Impedance coefficient
Z_coeff = (beta_0^2) / (4.0 * omega * eps_0);

% Discretization: λ/4 sampling
discretization_factor = 4.0;
DeltaX = lambda / discretization_factor;
```

#### 3.1.2 Hankel Function Implementation

The zeroth-order Hankel function of the second kind (line 75):

```matlab
calc_H02 = @(Arg) besselj(0, Arg) - 1j * bessely(0, Arg);
```

This is constructed from Bessel functions: H₀⁽²⁾(x) = J₀(x) - jY₀(x)

#### 3.1.3 Self-Impedance Calculation

Self-impedance uses numerical integration via `quadgk` (lines 77-89):

```matlab
fun_real = @(x) besselj(0, beta_0 * x);
fun_imag = @(x) bessely(0, beta_0 * x);

for i = 1:NoLinesubs
    L_half = seg_length(i) / 2.0;
    val_real = quadgk(fun_real, 0, L_half);
    val_imag = quadgk(fun_imag, 0, L_half);
    Zself(i) = 2.0 * Z_coeff * (val_real - 1j * val_imag);
end
```

#### 3.1.4 Forward-Backward Iterative Solution

**Forward Sweep** (lines 98-105):
```matlab
J(1) = Ei(1) / Zself(1);
for p = 2:NoLinesubs
    q_idx = 1:(p-1);
    dist_pq = sqrt((x_seg(q_idx) - x_seg(p)).^2 + (y_seg(q_idx) - y_seg(p)).^2);
    Z_pq = Z_coeff * calc_H02(beta_0 * dist_pq);
    SUM = sum(seg_length(q_idx) .* Z_pq .* J(q_idx));
    J(p) = (Ei(p) - SUM) / Zself(p);
end
```

**Backward Sweep** (lines 108-114):
```matlab
for p = (NoLinesubs-1):-1:1
    q_idx = (p+1):NoLinesubs;
    dist_pq = sqrt((x_seg(q_idx) - x_seg(p)).^2 + (y_seg(q_idx) - y_seg(p)).^2);
    Z_pq = Z_coeff * calc_H02(beta_0 * dist_pq);
    SUM = sum(seg_length(q_idx) .* Z_pq .* J(q_idx));
    J(p) = J(p) - SUM / Zself(p);
end
```

#### 3.1.5 Path Loss Calculation Strategy

The code implements a 2D-to-3D conversion strategy (lines 147-182):

1. Calculate 2D distances from source
2. Normalize field magnitudes by cylindrical spreading (√ρ)
3. Compute excess loss as the ratio of normalized incident to total field
4. Combine with 3D free-space path loss

```matlab
% Normalize by cylindrical spreading
Ei_normalized = mag_Ei .* sqrt(dist_2d);
Et_normalized = mag_Et .* sqrt(dist_2d);

% Excess loss extraction
L_excess = 20 * log10(Ei_normalized ./ Et_normalized);

% 3D Free Space Path Loss
FSPL_3D = 20 * log10(4 * pi * dist_3d / lambda);

% Total path loss
PL = FSPL_3D + L_excess;
```

---

### 3.2 Longley-Rice Irregular Terrain Model (ITM)

#### 3.2.1 ITM Wrapper Structure

The `itm_wrapper.m` function (lines 20-75) iterates through terrain points, calculating path loss from transmitter to each subsequent point:

```matlab
% Global terrain roughness for consistency
global_delta_h = dlthx(full_pfl, dist_m(1), dist_m(end));

% ITM configuration
itm_opts = struct();
itm_opts.pol = 0;           % Horizontal Polarization
itm_opts.clim = 5;          % Continental Temperate
itm_opts.conf = 0.5;        % 50% Confidence
itm_opts.force_delta_h = global_delta_h;
```

#### 3.2.2 Profile Analysis (qlrpfl.m)

The `qlrpfl` function computes geometric parameters:

**Effective Earth Radius** (lines 21-26):
```matlab
k_factor = 1.0 / (1.0 - 0.04665 * exp(N_s / 179.3));
a_earth = 6371000; % meters
prop.a_eff = k_factor * a_earth;
```

For N_s = 301 N-units: k ≈ 4/3, a_eff ≈ 8.495×10⁶ m

**Horizon Extraction** (hzns.m):
```matlab
for i = 2:np
    dist = (i - 1) * xi;
    % Elevation angle accounting for Earth curvature
    theta = (z(i) - h_tx_total) / dist - dist / (2 * a_eff);
    if theta > the_max
        the_max = theta;
        d_L = dist;
    end
end
```

**Terrain Irregularity (dlthx.m)**:
```matlab
% Range selection (10% to 90%)
idx_start = floor(0.1 * np) + 1;
idx_end = ceil(0.9 * np);

% Linear trend removal
coeffs = X \ z_sub(:);
z_trend = m * x_sub + c;
residuals = z_sub(:) - z_trend;

% Interdecile range
delta_h = v_90 - v_10;
```

**Effective Heights via Least Squares (zlsq1.m)**:
```matlab
% Fit Line L(x) = mx + c
X = [ones(n_seg, 1), x_seg];
coeffs = X \ z_seg(:);
z_eff = m * x_target + c;
```

#### 3.2.3 Propagation Physics (lrprop.m)

**Mode Selection Logic**:
- `dist < d_Lt`: Line-of-Sight with two-ray propagation
- `d_Lt < dist <= 1.5*d_Lt`: Diffraction region
- `dist > 1.5*d_Lt`: Troposcatter with minimum of scatter/diffraction

**Two-Ray Propagation (LOS Region)**:
```matlab
arg = 2 * pi * h_e(1) * h_e(2) / (lambda * dist);
theta_g = atan((h_e(1) + h_e(2)) / dist);
Gamma = calc_fresnel_reflection(theta_g, freq, prop_params);
E_magnitude = abs(1 + Gamma * exp(1j * arg));
A_two_ray = -20 * log10(E_magnitude);
```

**Modified Diffraction Weighting (Non-Standard ITM)**:

The implementation uses a continuous weighting function (lines 91-137):
```matlab
% Knife-edge diffraction
if v > -0.7
    A_ke = 6.9 + 20 * log10(sqrt((v-0.1)^2 + 1) + v - 0.1);
else
    A_ke = 0;
end

% Smooth-earth diffraction (Vogler's method)
beta = (1 + (prop.h_e(1) + prop.h_e(2)) / prop.a_eff) * d / prop.a_eff;
X = 2 * beta * sqrt(prop.a_eff * lambda / pi);
if X < 1.6
    A_r = 20 * log10(1 + X);
else
    A_r = 20 * log10(X) + 5.8;
end

% Continuous blending (non-standard modification)
k_rough = delta_h / lambda;
w = k_rough^2 / (1 + k_rough^2);
A_diff = (1 - w) * A_r + w * A_ke;
```

**Troposcatter with Climate Corrections** (lines 140-182):
```matlab
theta_d = d / prop.a_eff;
theta_s = max(0.001, theta_d - (h1_km + h2_km) / (prop.a_eff / 1000));
A_scat = 165 + 20*log10(f_ghz) + 30*log10(theta_d) - 10*log10(theta_s);

% Custom climate corrections (non-standard)
climate_factors = [0, -2, -4, +5, 0, -3, -5]; % dB adjustments
A_scat = A_scat + climate_factors(klim);
```

**Fresnel Reflection Coefficient** (lines 184-232):
```matlab
eps_complex = eps_r - 1j * sigma / (omega * eps_0);
sqrt_term = sqrt(eps_complex - cos_theta^2);

if pol == 0  % Horizontal
    Gamma = (sin_theta - sqrt_term) / (sin_theta + sqrt_term);
else  % Vertical
    Gamma = (eps_complex * sin_theta - sqrt_term) / ...
            (eps_complex * sin_theta + sqrt_term);
end
```

The complex permittivity formulation accounts for both displacement current (real part, ε_r) and conduction current (imaginary part, σ/(ωε₀)) in the ground medium. This enables accurate modeling of ground reflection effects across the frequency range, where the reflection coefficient magnitude approaches unity at grazing incidence and decreases toward the pseudo-Brewster angle for vertical polarization. The implementation correctly handles both horizontal (TE) and vertical (TM) polarization cases.

#### 3.2.4 Variability Calculation (avar.m)

**Location Variability**:
```matlab
k_f = 1 + log10(freq / 100);  % Frequency factor
sigma_L = 10 * (k_f * delta_h) / (k_f * delta_h + 13);
Y_L = sigma_L * z_c;
```

**Time Variability** (Climate-dependent lookup table):
```matlab
sigma_T_table = [
    5, 6, 7, 8;    % 1: Equatorial
    4, 5, 6, 7;    % 2: Continental Subtropical
    4, 5, 6, 7;    % 3: Maritime Subtropical
    6, 7, 8, 9;    % 4: Desert
    4, 5, 6, 7;    % 5: Continental Temperate
    3, 4, 5, 6;    % 6: Maritime Temperate Land
    2, 3, 4, 5;    % 7: Maritime Temperate Sea
];
```

**Combined Variability**:
```matlab
Y_total = sign(z_c) * sqrt((sigma_T*z_c)^2 + (sigma_L*z_c)^2 + (sigma_S*z_c)^2);
A_total = A_ref + Y_total;
```

---

### 3.3 Bullington Knife-Edge Diffraction

#### 3.3.1 Algorithm Implementation (bullington_wrapper.m)

**Earth Curvature Correction** (lines 66-70):
```matlab
[R_effective, ~] = get_effective_earth_radius(301);
bulge = (d1 .* d2) / (2 * R_effective);
h_eff = z_segment + bulge;
```

**Maximum Slope Calculation** (lines 79-88):
```matlab
% Slopes from Tx
slopes_tx = (h_mid - h_tx_ant) ./ x_mid;
[m_tx, ~] = max(slopes_tx);

% Slopes from Rx (looking back)
dist_from_rx = D_total - x_mid;
slopes_rx = (h_mid - h_rx_ant) ./ dist_from_rx;
[m_rx_prime, ~] = max(slopes_rx);
```

**Equivalent Obstacle Intersection** (lines 95-103):
```matlab
denominator = m_tx + m_rx_prime;
x_eq = (h_rx_ant - h_tx_ant + m_rx_prime * D_total) / denominator;
h_eq = h_tx_ant + m_tx * x_eq;

slope_los = (h_rx_ant - h_tx_ant) / D_total;
h_los_eq = h_tx_ant + slope_los * x_eq;
h_obs = h_eq - h_los_eq;
```

**Fresnel Parameter and Diffraction Loss** (lines 109-121):
```matlab
nu = h_obs * sqrt(2 * (d1_eq + d2_eq) / (lambda * d1_eq * d2_eq));

% ITU-R P.526 diffraction loss
if nu > -0.78
    J_nu = 6.9 + 20 * log10(sqrt((nu - 0.1)^2 + 1) + nu - 0.1);
    L_diff = J_nu;
end

PL(i) = fspl + L_diff;
```

---

### 3.4 Deygout Multiple-Edge Diffraction

#### 3.4.1 Three-Edge Algorithm (deygout_wrapper.m)

**Principal Edge Identification** (lines 105-148):
```matlab
function [max_nu, max_idx] = find_max_nu(dists, heights, tx_h, rx_h, lambda)
    [R_e, ~] = get_effective_earth_radius(301);
    
    for i = 2:num_p-1
        d1 = dists(i) - dists(1);
        d2 = dists(end) - dists(i);
        
        % LOS height
        h_los = h_tx_abs + (h_rx_abs - h_tx_abs) * (d1 / d_total);
        
        % Earth bulge
        h_bulge = (d1 * d2) / (2 * R_e);
        
        % Obstacle clearance
        h = heights(i) + h_bulge - h_los;
        nu = h * sqrt(2 * (d1 + d2) / (lambda * d1 * d2));
        
        if nu > max_nu
            max_nu = nu;
            max_idx = i;
        end
    end
end
```

**Recursive Sub-Path Analysis** (lines 54-103):
```matlab
% Principal Edge
[nu_p, idx_p] = find_max_nu(dists, heights, tx_h, rx_h, lambda);
J_p = diffraction_loss_db(nu_p);

% Left Edge (Tx to Principal)
[nu_t, ~] = find_max_nu(dists(1:idx_p), heights(1:idx_p), tx_h, 0, lambda);
J_t = diffraction_loss_db(nu_t);

% Right Edge (Principal to Rx)
[nu_r, ~] = find_max_nu(dists(idx_p:end), heights(idx_p:end), 0, rx_h, lambda);
J_r = diffraction_loss_db(nu_r);
```

**Empirical Corrections** (lines 95-102):
```matlab
% Transition factor
if J_p <= 6
    T = J_p / 6.0;
else
    T = 1.0;
end

% Distance-dependent correction
C = 8.0 + 0.04 * total_dist_km;

% Total loss
L_total = J_p + T * (J_t + J_r + C);
```

---

### 3.5 Okumura-Hata and COST-231 Empirical Models

#### 3.5.1 Model Selection Logic (hata_model.m)

```matlab
use_cost231 = (f > 1500);
```

#### 3.5.2 Mobile Antenna Height Correction

**Urban Large City** (lines 60-69):
```matlab
if env == "urban_large"
    if f <= 200
        a_hr = 8.29 * (log10(1.54 * h_rx))^2 - 1.1;
    else
        a_hr = 3.2 * (log10(11.75 * h_rx))^2 - 4.97;
    end
else
    a_hr = (1.1 * log_f - 0.7) * h_rx - (1.56 * log_f - 0.8);
end
```

#### 3.5.3 Standard Okumura-Hata (150-1500 MHz)

**Urban Path Loss** (lines 113-122):
```matlab
A = 69.55; B = 26.16; C = 13.82; D = 44.9; E = 6.55;
L_urban = A + B * log_f - C * log_ht - a_hr + (D - E * log_ht) * log_d;
```

**Suburban Correction** (line 129):
```matlab
PL = L_urban - 2 * (log10(f / 28))^2 - 5.4;
```

**Rural Correction** (line 133):
```matlab
PL = L_urban - 4.78 * (log_f)^2 + 18.33 * log_f - 40.94;
```

#### 3.5.4 COST-231 Extension (1500-2000 MHz)

```matlab
A = 46.3; B = 33.9; C = 13.82; D = 44.9; E = 6.55;
C_m = 3;  % Metropolitan center correction
L_base = A + B * log_f - C * log_ht - a_hr + (D - E * log_ht) * log_d + C_m;
```

---

## 4. Comparative Analysis Framework

### 4.1 Statistical Metrics (main_propagation_comparison.m)

**RMSE Calculation** (line 196):
```matlab
rmse = @(est, ref) sqrt(nanmean((est - ref).^2));
```

**Bias Calculation** (line 197):
```matlab
bias = @(est, ref) nanmean(est - ref);
```

### 4.2 Reference Model Selection

The EFIE model serves as the reference standard when enabled:
```matlab
if Scenario.Run_EFIE
    Results.EFIE_dB = efie_wrapper(Scenario, TerrainProfile);
end
```

### 4.3 Output Generation

- Comparison plot: `output/comparison_plot.png`
- Terrain geometry: `output/terrain_geometry.png`
- Metrics CSV: `output/metrics.csv`

---

## 5. Simulation Parameters and Configuration

### 5.1 Default Parameters (main_propagation_comparison.m)

```matlab
Scenario.Frequency_MHz = 970;       % UHF frequency
Scenario.Tx_Height_m = 52.0;        % Transmitter height (AGL)
Scenario.Rx_Height_m = 2.4;         % Receiver height (AGL)
Scenario.Environment = 'suburban';  % Hata environment classification
```

### 5.2 ITM Ground Parameters (longley_rice_p2p.m)

```matlab
options.pol = 0;           % Horizontal Polarization
options.eps_r = 15;        % Dielectric constant
options.sigma = 0.005;     % Conductivity (S/m)
options.N_s = 301;         % Surface refractivity (N-units)
options.clim = 5;          % Continental Temperate climate
options.conf = 0.5;        % 50% Confidence level
```

### 5.3 Standardized Effective Earth Radius (get_effective_earth_radius.m)

```matlab
a_earth = 6371000;  % meters
k_factor = 1.0 / (1.0 - 0.04665 * exp(N_s / 179.3));
a_eff = k_factor * a_earth;
```

For N_s = 301: k ≈ 1.333 (4/3), a_eff ≈ 8.495×10⁶ m

---

## 6. Integration Mechanisms

### 6.1 File-Based Data Exchange

The Python terrain extraction module writes output in the format:
```
distance_m elevation_m
0.0 390.0
10.0 390.0
...
```

The MATLAB master harness reads this format:
```matlab
raw_data = readmatrix(Scenario.Terrain_File_Path);
TerrainProfile.Distance_m = raw_data(:, 1);
TerrainProfile.Elevation_m = raw_data(:, 2);
```

### 6.2 Wrapper Function Interface

All propagation models implement a standardized interface:
```matlab
PL = model_wrapper(Scenario, TerrainProfile)
```

Where:
- `Scenario`: Structure with frequency, antenna heights, environment
- `TerrainProfile`: Structure with distance and elevation vectors
- `PL`: Column vector of path loss values in dB

### 6.3 Path Configuration (setup_project.m)

```matlab
paths_to_add = {
    fullfile(root_dir, 'api'),
    fullfile(root_dir, 'models', 'itm'),
    fullfile(root_dir, 'models', 'empirical'),
    fullfile(root_dir, 'models', 'diffraction'),
    fullfile(root_dir, 'models', 'fullwave'),
    fullfile(root_dir, 'models', 'utils'),
    fullfile(root_dir, 'data'),
    fullfile(root_dir, 'output')
};
```

---

## 7. Technical Specifications

### 7.1 MATLAB Requirements

Based on functions used:
- `readmatrix` (R2019a+) - Required for terrain data loading
- `quadgk` (R2007b+) - Adaptive Gauss-Kronrod quadrature integration
- `besselj`, `bessely` (R2006a+) - Bessel functions of first and second kind
- `interp1` (R2006a+) - Linear interpolation
- `norminv` (Statistics Toolbox, R2006a+) - Inverse normal distribution

Minimum version: **MATLAB R2019a** with Statistics Toolbox (due to `readmatrix` dependency)

### 7.2 Python Requirements

```python
import ee                    # earthengine-api
import math
from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict, Any
from enum import Enum
import json
import os
```

Requires: **Python 3.7+** with `earthengine-api` package

### 7.3 Google Earth Engine Configuration

- Google Cloud Project with Earth Engine API enabled
- Authentication via `ee.Authenticate()` and `ee.Initialize(project=project_id)`
- Environment variable: `EE_PROJECT_ID` or `GOOGLE_CLOUD_PROJECT`

---

## 8. Model Limitations and Assumptions

### 8.1 EFIE Model

- Assumes 2D geometry (infinite extent perpendicular to propagation plane)
- Perfect electric conductor (PEC) ground boundary
- Converts 2D field solutions to 3D path loss through normalization
- Computational complexity: O(N²) with profile discretization
- Discretization: λ/4 sampling for accuracy

### 8.2 ITM Model

- Modified implementation with continuous diffraction weighting (non-standard)
- Allows negative reference attenuation (signal enhancement)
- Custom climate corrections for troposcatter (non-standard)
- Valid frequency range: 20 MHz to 20 GHz
- Valid path length: 1 to 2000 km

### 8.3 Bullington Model

- Reduces terrain to single equivalent knife-edge
- Most accurate for profiles with one dominant obstacle
- Does not account for atmospheric effects

### 8.4 Deygout Model

- Limited to three edges (principal + two secondary)
- Empirical correction factor: C = 8.0 + 0.04×D_km
- Accuracy degrades for profiles with many comparable obstacles

### 8.5 Hata/COST-231 Models

- Terrain-independent (ignores elevation profile)
- Statistical ensemble predictions (not site-specific)
- Valid ranges:
  - Okumura-Hata: 150-1500 MHz
  - COST-231: 1500-2000 MHz
- Environment classifications: urban, urban_large, suburban, rural

---

## 9. Code Quality and Reproducibility

### 9.1 Physical Constants

All physical constants are computed from fundamental values:
```matlab
c = 299792458;           % Exact speed of light
mu_0 = 4*pi*1e-7;        % Exact permeability
eps_0 = 8.854e-12;       % Approximate permittivity
```

### 9.2 Platform Independence

File paths use MATLAB's `fullfile` function for cross-platform compatibility:
```matlab
Scenario.Terrain_File_Path = fullfile('data', 'terrain_profile.txt');
```

### 9.3 Data Validation

Input sanitization removes NaN values:
```matlab
valid_mask = ~isnan(TerrainProfile.Elevation_m) & ~isnan(TerrainProfile.Distance_m);
TerrainProfile.Distance_m = TerrainProfile.Distance_m(valid_mask);
TerrainProfile.Elevation_m = TerrainProfile.Elevation_m(valid_mask);
```

---

## 10. Summary

This technical report documents the actual implementation of a unified radio propagation computational framework. The framework integrates:

1. **Python-based terrain extraction** using Google Earth Engine API with SRTM datasets
2. **Five propagation models** with standardized wrapper interfaces:
   - EFIE full-wave solver (reference)
   - Longley-Rice ITM (modified implementation)
   - Bullington knife-edge diffraction
   - Deygout multiple-edge diffraction
   - Okumura-Hata/COST-231 empirical models
3. **Comparative analysis** using RMSE and bias metrics

Key implementation characteristics:
- λ/4 discretization for EFIE accuracy
- Continuous diffraction weighting in ITM (non-standard)
- Standardized effective Earth radius across all models
- File-based integration between Python and MATLAB components

---

*Report generated based on code analysis of the ComputationalFramework repository.*
