# IV. IMPLEMENTATION

This chapter presents the detailed implementation of the computational framework developed for comparative analysis of radio propagation models over irregular terrain. The implementation comprises two distinct software modules: a terrain data extraction pipeline developed in Python interfacing with Google Earth Engine, and a MATLAB-based propagation modeling engine implementing five prediction algorithms. The following sections describe the software architecture, algorithmic implementations, and integration mechanisms.

## A. Software Development Environment

### 1) Development Platform and Tools

The computational framework was developed using a dual-language approach to leverage the strengths of both Python and MATLAB environments. Python 3.7 or later was selected for terrain data acquisition due to its superior integration with cloud-based geospatial services, while MATLAB R2019a or later was chosen for numerical computation due to its optimized matrix operations and signal processing toolboxes.

The development environment configuration includes:

- **Python Environment**: Python 3.7+ with the `earthengine-api` package for Google Earth Engine access, utilizing the `dataclasses` module for structured data representation and `enum` for type-safe dataset selection.

- **MATLAB Environment**: MATLAB R2019a or later with the Statistics and Machine Learning Toolbox (required for `norminv` function in variability calculations). Core numerical functions including `quadgk` (adaptive quadrature), `besselj`/`bessely` (Bessel functions), and `interp1` (interpolation) are utilized from the base MATLAB installation.

- **Version Control**: Git-based version control enables tracking of algorithm modifications and collaborative development.

### 2) Project Architecture

The framework employs a modular architecture facilitating independent development and testing of individual components. The directory structure organizes code by functional category:

```
ComputationalFramework/
├── api/                          # Interface layer
│   ├── dem_profile_api.py        # Terrain extraction module
│   └── main_propagation_comparison.m  # Master orchestration script
├── models/                       # Propagation model implementations
│   ├── itm/                      # Longley-Rice model (8 files)
│   ├── empirical/                # Hata/COST-231 models (4 files)
│   ├── diffraction/              # Bullington and Deygout (4 files)
│   ├── fullwave/                 # EFIE solver (5 files)
│   └── utils/                    # Shared utility functions
├── data/                         # Input terrain profiles
├── output/                       # Generated results and visualizations
└── setup_project.m               # Path configuration script
```

Each propagation model is encapsulated within a standardized wrapper function accepting common input structures, enabling uniform invocation from the master harness. This abstraction layer permits addition of new models without modification to the orchestration logic.

## B. Terrain Data Extraction Implementation

### 1) Google Earth Engine Integration

The terrain extraction module (`dem_profile_api.py`) implements a class-based architecture utilizing Python dataclasses for type-safe data representation. The `DEMProfileExtractor` class manages Earth Engine authentication, coordinate interpolation, and elevation sampling.

**Authentication and Initialization:**

```python
def _initialize_ee(self) -> None:
    """Initialize Google Earth Engine with project credentials."""
    if not self.project:
        self.project = input("Enter your Google Cloud Project ID: ").strip()
    
    try:
        ee.Initialize(project=self.project)
        DEMProfileExtractor._global_initialized = True
    except ee.EEException as e:
        if "authenticate" in str(e).lower():
            ee.Authenticate()
            ee.Initialize(project=self.project)
```

The implementation supports multiple Digital Elevation Model (DEM) datasets through an enumeration pattern:

| Dataset Identifier | Source | Resolution | Coverage |
|-------------------|--------|------------|----------|
| `USGS/SRTMGL1_003` | SRTM Global 1 arc-second | ~30 m | 60°N to 56°S |
| `CGIAR/SRTM90_V4` | SRTM 90m | ~90 m | 60°N to 56°S |
| `JAXA/ALOS/AW3D30/V3_2` | ALOS World 3D | ~30 m | Global |
| `COPERNICUS/DEM/GLO30` | Copernicus GLO-30 | ~30 m | Global |
| `NASA/NASADEM_HGT/001` | NASADEM | ~30 m | 60°N to 56°S |

### 2) Geodesic Path Computation

Profile extraction between transmitter and receiver coordinates employs great-circle path interpolation. The implementation computes intermediate coordinates through spherical trigonometry.

**Haversine Distance Calculation:**

For two geographic coordinates $(φ_1, λ_1)$ and $(φ_2, λ_2)$, the great-circle distance is computed as:

$$d = 2R \cdot \arcsin\left(\sqrt{\sin^2\left(\frac{φ_2 - φ_1}{2}\right) + \cos(φ_1)\cos(φ_2)\sin^2\left(\frac{λ_2 - λ_1}{2}\right)}\right)$$

where $R = 6{,}371{,}000$ meters represents Earth's mean radius. The Python implementation:

```python
@staticmethod
def haversine_distance(coord1: Coordinate, coord2: Coordinate) -> float:
    R = 6371000  # Earth radius in meters
    
    lat1_rad = math.radians(coord1.latitude)
    lat2_rad = math.radians(coord2.latitude)
    delta_lat = math.radians(coord2.latitude - coord1.latitude)
    delta_lon = math.radians(coord2.longitude - coord1.longitude)
    
    a = (math.sin(delta_lat / 2) ** 2 + 
         math.cos(lat1_rad) * math.cos(lat2_rad) * 
         math.sin(delta_lon / 2) ** 2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    
    return R * c
```

**Spherical Linear Interpolation (Slerp):**

To generate $N$ uniformly spaced points along the geodesic path, the implementation applies spherical linear interpolation. For interpolation parameter $f \in [0, 1]$:

$$\vec{P}(f) = \frac{\sin((1-f)θ)}{\sin θ}\vec{P}_1 + \frac{\sin(fθ)}{\sin θ}\vec{P}_2$$

where $θ$ represents the angular distance between endpoints and $\vec{P}_1$, $\vec{P}_2$ are the unit vectors corresponding to the start and end coordinates respectively. The Cartesian coordinates are computed as:

```python
A = math.sin((1 - fraction) * d) / math.sin(d)
B = math.sin(fraction * d) / math.sin(d)

x = A * math.cos(lat1) * math.cos(lon1) + B * math.cos(lat2) * math.cos(lon2)
y = A * math.cos(lat1) * math.sin(lon1) + B * math.cos(lat2) * math.sin(lon2)
z = A * math.sin(lat1) + B * math.sin(lat2)

lat = math.degrees(math.atan2(z, math.sqrt(x**2 + y**2)))
lon = math.degrees(math.atan2(y, x))
```

### 3) Batch Elevation Sampling

To optimize API throughput while maintaining query reliability, the implementation retrieves elevation data in batches of 100 coordinates per API call:

```python
def _fetch_elevations_batch(self, coordinates: List[Coordinate], 
                            batch_size: int = 100) -> List[Optional[float]]:
    dem = ee.Image(self.dataset.value)
    band = self._get_elevation_band()
    
    for batch_start in range(0, len(coordinates), batch_size):
        batch_coords = coordinates[batch_start:batch_start + batch_size]
        
        points = ee.FeatureCollection([
            ee.Feature(ee.Geometry.Point(coord.to_list()), {'index': i})
            for i, coord in enumerate(batch_coords)
        ])
        
        sampled = dem.select(band).sampleRegions(
            collection=points,
            scale=30,
            geometries=False
        )
        
        results = sampled.getInfo()
        # Process results...
```

### 4) Output Format Specification

The terrain extraction module generates ASCII text files with space-delimited columns:

- **Column 1**: Cumulative distance from transmitter (meters)
- **Column 2**: Terrain elevation above WGS84 ellipsoid (meters)

```
0.0 390.0
10.0 389.5
20.0 388.2
...
```

This format specification ensures direct compatibility with the MATLAB propagation solvers without requiring intermediate data transformation.

## C. EFIE Full-Wave Solver Implementation

### 1) Problem Formulation

The Electric Field Integral Equation (EFIE) implementation solves the two-dimensional scattering problem for a perfectly electrically conducting (PEC) surface illuminated by an infinite line source. The surface current distribution $J_s$ on the terrain surface satisfies:

$$\hat{n} \times \vec{E}_{inc}(\vec{r}) = -\hat{n} \times \vec{E}_{scat}(\vec{r}), \quad \vec{r} \in S$$

where $\vec{E}_{inc}$ is the incident field from the line source and $\vec{E}_{scat}$ is the scattered field due to induced surface currents.

### 2) Method of Moments Discretization

The terrain profile is discretized into $N$ linear segments with segment length determined by wavelength-dependent sampling:

```matlab
discretization_factor = 4.0;          % λ/4 sampling
DeltaX = lambda / discretization_factor;
NoLinesubs = floor(TerrainLength / DeltaX);
```

The λ/4 discretization satisfies Nyquist sampling requirements for accurate current representation, yielding segment count:

$$N = \left\lfloor \frac{L_{terrain}}{\lambda/4} \right\rfloor = \left\lfloor \frac{4L_{terrain}}{\lambda} \right\rfloor$$

For a 700-meter terrain profile at 970 MHz ($\lambda \approx 0.309$ m), this produces approximately 9,061 unknowns.

### 3) Impedance Matrix Construction

**Mutual Impedance:**

The mutual impedance between segments $p$ and $q$ is computed using the zeroth-order Hankel function of the second kind:

$$Z_{pq} = Z_{coeff} \cdot H_0^{(2)}(\beta_0 \rho_{pq})$$

where:
- $Z_{coeff} = \frac{\beta_0^2}{4\omega\epsilon_0}$ is the impedance coefficient
- $\beta_0 = \omega\sqrt{\mu_0\epsilon_0} = \frac{2\pi}{\lambda}$ is the free-space wavenumber
- $\rho_{pq} = \sqrt{(x_p - x_q)^2 + (y_p - y_q)^2}$ is the inter-segment distance
- $H_0^{(2)}(x) = J_0(x) - jY_0(x)$ is constructed from Bessel functions

The MATLAB implementation:

```matlab
% Physical constants
c = 299792458;
mu_0 = 4*pi*1e-7;
eps_0 = 8.854e-12;
lambda = c / f;
omega = 2 * pi * f;
beta_0 = omega * sqrt(mu_0 * eps_0);

% Impedance coefficient
Z_coeff = (beta_0^2) / (4.0 * omega * eps_0);

% Hankel function H0^(2)
calc_H02 = @(Arg) besselj(0, Arg) - 1j * bessely(0, Arg);

% Mutual impedance
Z_pq = Z_coeff * calc_H02(beta_0 * dist_pq);
```

**Self-Impedance:**

Self-impedance calculation requires numerical integration to handle the logarithmic singularity of $Y_0(x)$ as $x \to 0$:

$$Z_{self} = 2Z_{coeff} \int_0^{L/2} \left[J_0(\beta_0 x) - jY_0(\beta_0 x)\right] dx$$

The implementation employs MATLAB's adaptive Gauss-Kronrod quadrature:

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

### 4) Forward-Backward Iterative Solution

Rather than direct matrix inversion with $O(N^3)$ complexity, the implementation employs the Forward-Backward Method (FBM) achieving $O(N^2)$ complexity.

**Forward Sweep:**

The forward sweep computes preliminary current values marching from the source location:

$$J_{fwd}(p) = \frac{E_{inc}(p) - \sum_{q=1}^{p-1} J(q) Z_{pq} \ell_q}{Z_{self,p}}$$

```matlab
J(1) = Ei(1) / Zself(1);
for p = 2:NoLinesubs
    q_idx = 1:(p-1);
    dist_pq = sqrt((x_seg(q_idx) - x_seg(p)).^2 + ...
                   (y_seg(q_idx) - y_seg(p)).^2);
    Z_pq = Z_coeff * calc_H02(beta_0 * dist_pq);
    SUM = sum(seg_length(q_idx) .* Z_pq .* J(q_idx));
    J(p) = (Ei(p) - SUM) / Zself(p);
end
```

**Backward Sweep:**

The backward sweep refines the solution by incorporating contributions from segments ahead:

$$J(p) = J_{fwd}(p) - \frac{\sum_{q=p+1}^{N} J(q) Z_{pq} \ell_q}{Z_{self,p}}$$

```matlab
for p = (NoLinesubs-1):-1:1
    q_idx = (p+1):NoLinesubs;
    dist_pq = sqrt((x_seg(q_idx) - x_seg(p)).^2 + ...
                   (y_seg(q_idx) - y_seg(p)).^2);
    Z_pq = Z_coeff * calc_H02(beta_0 * dist_pq);
    SUM = sum(seg_length(q_idx) .* Z_pq .* J(q_idx));
    J(p) = J(p) - SUM / Zself(p);
end
```

### 5) Path Loss Extraction

The two-dimensional EFIE solution provides field magnitudes with cylindrical spreading characteristics ($1/\sqrt{\rho}$ dependence). To extract terrain-specific excess loss for combination with three-dimensional propagation:

**Step 1: Compute Scattered Field**

```matlab
for k = 1:num_points
    x_obs = X_terrain(k);
    y_obs = Y_terrain(k) + rx_h;
    
    dist_seg_obs = sqrt((x_seg - x_obs).^2 + (y_seg - y_obs).^2);
    Z_seg_obs = Z_coeff * calc_H02(beta_0 * dist_seg_obs);
    E_scattered = sum(J .* seg_length .* Z_seg_obs);
    
    Et_total(k) = Ei_total(k) - E_scattered;
end
```

**Step 2: Normalize by Cylindrical Spreading**

```matlab
% Remove cylindrical spreading to extract excess loss
Ei_normalized = mag_Ei .* sqrt(dist_2d);
Et_normalized = mag_Et .* sqrt(dist_2d);

% Excess loss (terrain interaction only)
L_excess = 20 * log10(Ei_normalized ./ Et_normalized);
```

**Step 3: Combine with 3D Free-Space Path Loss**

$$PL_{total} = 20\log_{10}\left(\frac{4\pi d}{\lambda}\right) + L_{excess}$$

```matlab
FSPL_3D = 20 * log10(4 * pi * dist_3d / lambda);
PL = FSPL_3D + L_excess;
```

## D. Longley-Rice ITM Implementation

### 1) Model Overview

The Longley-Rice Irregular Terrain Model (ITM) implementation follows the point-to-point prediction mode documented in NTIA Technical Report TR-ERL 79-ITS 67. The implementation includes several modifications from standard ITM variants as described below.

### 2) Profile Analysis (`qlrpfl.m`)

**Effective Earth Radius:**

The effective Earth radius accounts for atmospheric refraction:

$$a_{eff} = \frac{a_{earth}}{1 - 0.04665 \cdot \exp(N_s/179.3)}$$

where $a_{earth} = 6{,}371{,}000$ meters and $N_s$ represents surface refractivity (default: 301 N-units for continental temperate climate).

```matlab
k_factor = 1.0 / (1.0 - 0.04665 * exp(N_s / 179.3));
a_earth = 6371000;
prop.a_eff = k_factor * a_earth;
```

For $N_s = 301$: $k \approx 4/3$, yielding $a_{eff} \approx 8.495 \times 10^6$ meters.

**Horizon Search (`hzns.m`):**

Radio horizons are identified by finding the maximum grazing angle from each terminal:

```matlab
for i = 2:np
    dist = (i - 1) * xi;
    % Elevation angle with Earth curvature correction
    theta = (z(i) - h_tx_total) / dist - dist / (2 * a_eff);
    if theta > the_max
        the_max = theta;
        d_L = dist;
    end
end
```

**Terrain Irregularity (`dlthx.m`):**

The interdecile range parameter $\Delta h$ characterizes terrain roughness:

```matlab
% Select central 80% of profile
idx_start = floor(0.1 * np) + 1;
idx_end = ceil(0.9 * np);
z_sub = z(idx_start:idx_end);

% Linear trend removal
X = [ones(n_sub, 1), x_sub];
coeffs = X \ z_sub(:);
z_trend = coeffs(2) * x_sub + coeffs(1);
residuals = z_sub(:) - z_trend;

% Interdecile range
sorted_res = sort(residuals);
delta_h = sorted_res(round(0.9*n_sub)) - sorted_res(round(0.1*n_sub));
```

### 3) Propagation Mode Selection (`lrprop.m`)

The core prediction engine selects between three propagation modes based on path geometry:

| Distance Range | Mode | Primary Mechanism |
|---------------|------|-------------------|
| $d < d_{LOS}$ | Line-of-Sight | Two-ray ground reflection |
| $d_{LOS} < d \leq 1.5d_{LOS}$ | Diffraction | Knife-edge and smooth-earth |
| $d > 1.5d_{LOS}$ | Troposcatter | Forward scatter |

**Line-of-Sight Region:**

Two-ray propagation includes ground reflection with Fresnel coefficient computation:

```matlab
arg = 2 * pi * h_e(1) * h_e(2) / (lambda * dist);
theta_g = atan((h_e(1) + h_e(2)) / dist);
Gamma = calc_fresnel_reflection(theta_g, freq, prop_params);
E_magnitude = abs(1 + Gamma * exp(1j * arg));
A_two_ray = -20 * log10(E_magnitude);
```

The Fresnel reflection coefficient accounts for complex ground permittivity:

$$\varepsilon_{complex} = \varepsilon_r - j\frac{\sigma}{\omega\varepsilon_0}$$

```matlab
eps_complex = eps_r - 1j * sigma / (omega * eps_0);
sqrt_term = sqrt(eps_complex - cos_theta^2);

if pol == 0  % Horizontal polarization
    Gamma = (sin_theta - sqrt_term) / (sin_theta + sqrt_term);
else  % Vertical polarization
    Gamma = (eps_complex * sin_theta - sqrt_term) / ...
            (eps_complex * sin_theta + sqrt_term);
end
```

**Diffraction Region (Modified Implementation):**

The implementation employs a continuous terrain roughness weighting scheme rather than discrete mode switching:

$$w = \frac{k^2}{1 + k^2}, \quad k = \frac{\Delta h}{\lambda}$$

$$A_{diff} = (1 - w) \cdot A_{smooth} + w \cdot A_{knife-edge}$$

This modification eliminates prediction discontinuities when terrain roughness crosses threshold values:

```matlab
% Knife-edge diffraction
v = theta_tot * sqrt(d / lambda);
if v > -0.7
    A_ke = 6.9 + 20 * log10(sqrt((v-0.1)^2 + 1) + v - 0.1);
else
    A_ke = 0;
end

% Smooth-earth diffraction (Vogler's method)
beta = (1 + (h_e(1) + h_e(2)) / a_eff) * d / a_eff;
X = 2 * beta * sqrt(a_eff * lambda / pi);
if X < 1.6
    A_r = 20 * log10(1 + X);
else
    A_r = 20 * log10(X) + 5.8;
end

% Continuous blending
k_rough = delta_h / lambda;
w = k_rough^2 / (1 + k_rough^2);
A_diff = (1 - w) * A_r + w * A_ke;
```

### 4) Variability Calculation (`avar.m`)

Statistical variability accounts for location, time, and situation uncertainties:

$$Y_{total} = \text{sign}(z_c) \sqrt{(\sigma_T z_c)^2 + (\sigma_L z_c)^2 + (\sigma_S z_c)^2}$$

where $z_c$ is the standard normal deviate for the specified confidence level.

**Location Variability:**

```matlab
k_f = 1 + log10(freq / 100);  % Frequency factor
sigma_L = 10 * (k_f * delta_h) / (k_f * delta_h + 13);
Y_L = sigma_L * z_c;
```

**Time Variability (Climate-Dependent):**

```matlab
sigma_T_table = [
    5, 6, 7, 8;    % Equatorial
    4, 5, 6, 7;    % Continental Subtropical
    4, 5, 6, 7;    % Maritime Subtropical
    6, 7, 8, 9;    % Desert
    4, 5, 6, 7;    % Continental Temperate
    3, 4, 5, 6;    % Maritime Temperate Land
    2, 3, 4, 5;    % Maritime Temperate Sea
];
```

## E. Knife-Edge Diffraction Implementations

### 1) Bullington Method (`bullington_wrapper.m`)

The Bullington method reduces arbitrary terrain profiles to a single equivalent knife-edge obstacle through geometric construction.

**Maximum Slope Calculation:**

```matlab
% Forward slopes from transmitter
slopes_tx = (h_mid - h_tx_ant) ./ x_mid;
[m_tx, ~] = max(slopes_tx);

% Backward slopes from receiver
dist_from_rx = D_total - x_mid;
slopes_rx = (h_mid - h_rx_ant) ./ dist_from_rx;
[m_rx_prime, ~] = max(slopes_rx);
```

**Equivalent Obstacle Location:**

The intersection of maximum-slope tangent lines determines the equivalent obstacle:

```matlab
x_eq = (h_rx_ant - h_tx_ant + m_rx_prime * D_total) / (m_tx + m_rx_prime);
h_eq = h_tx_ant + m_tx * x_eq;
h_obs = h_eq - h_los_eq;  % Clearance above LOS
```

**Fresnel Parameter and Diffraction Loss:**

$$\nu = h_{obs} \sqrt{\frac{2(d_1 + d_2)}{\lambda d_1 d_2}}$$

$$J(\nu) = \begin{cases} 
0 & \nu \leq -0.78 \\
6.9 + 20\log_{10}\left(\sqrt{(\nu-0.1)^2 + 1} + \nu - 0.1\right) & \nu > -0.78
\end{cases}$$

### 2) Deygout Method (`deygout_wrapper.m`)

The Deygout method applies a recursive three-edge approach for multiple terrain obstacles.

**Principal Edge Identification:**

```matlab
function [max_nu, max_idx] = find_max_nu(dists, heights, tx_h, rx_h, lambda)
    for i = 2:num_p-1
        d1 = dists(i) - dists(1);
        d2 = dists(end) - dists(i);
        
        h_los = h_tx_abs + (h_rx_abs - h_tx_abs) * (d1 / d_total);
        h_bulge = (d1 * d2) / (2 * R_e);
        h = heights(i) + h_bulge - h_los;
        
        nu = h * sqrt(2 * (d1 + d2) / (lambda * d1 * d2));
        
        if nu > max_nu
            max_nu = nu;
            max_idx = i;
        end
    end
end
```

**Three-Edge Loss Calculation:**

$$L_{total} = J(\nu_p) + T \cdot [J(\nu_t) + J(\nu_r) + C]$$

where:
- $\nu_p$: Principal edge Fresnel parameter
- $\nu_t$, $\nu_r$: Secondary edge parameters (Tx-side and Rx-side)
- $T = \min(J(\nu_p)/6, 1)$: Transition factor
- $C = 8.0 + 0.04 \cdot D_{km}$: Distance-dependent correction

## F. Empirical Model Implementations

### 1) Okumura-Hata Model (`hata_model.m`)

For frequencies 150-1500 MHz, the implementation applies the standard Okumura-Hata formulation:

$$L_U = 69.55 + 26.16\log_{10}(f) - 13.82\log_{10}(h_B) - a(h_M) + [44.9 - 6.55\log_{10}(h_B)]\log_{10}(d)$$

**Mobile Antenna Height Correction:**

For small to medium-sized cities:
$$a(h_M) = [1.1\log_{10}(f) - 0.7]h_M - [1.56\log_{10}(f) - 0.8]$$

For large cities ($f \geq 300$ MHz):
$$a(h_M) = 3.2[\log_{10}(11.75 h_M)]^2 - 4.97$$

**Environment Corrections:**

```matlab
% Suburban
PL = L_urban - 2 * (log10(f / 28))^2 - 5.4;

% Rural (Open Area)
PL = L_urban - 4.78 * (log10(f))^2 + 18.33 * log10(f) - 40.94;
```

### 2) COST-231 Extension

For frequencies 1500-2000 MHz:

$$L = 46.3 + 33.9\log_{10}(f) - 13.82\log_{10}(h_B) - a(h_M) + [44.9 - 6.55\log_{10}(h_B)]\log_{10}(d) + C_m$$

where $C_m = 0$ dB for suburban/rural and $C_m = 3$ dB for metropolitan centers.

## G. Master Harness and Integration

### 1) Scenario Configuration

The master harness (`main_propagation_comparison.m`) configures simulation parameters:

```matlab
Scenario = struct();
Scenario.Frequency_MHz = 970;       % UHF frequency
Scenario.Tx_Height_m = 52.0;        % Transmitter height (AGL)
Scenario.Rx_Height_m = 2.4;         % Receiver height (AGL)
Scenario.Environment = 'suburban';  % Hata environment classification
```

### 2) Standardized Wrapper Interface

All propagation models implement a common interface:

```matlab
function PL = model_wrapper(Scenario, TerrainProfile)
% Input:
%   Scenario - Structure with Frequency_MHz, Tx_Height_m, Rx_Height_m
%   TerrainProfile - Structure with Distance_m, Elevation_m vectors
%
% Output:
%   PL - Column vector of path loss values (dB)
```

### 3) Statistical Metrics Computation

Model accuracy evaluation employs RMSE and bias metrics:

```matlab
rmse = @(est, ref) sqrt(nanmean((est - ref).^2));
bias = @(est, ref) nanmean(est - ref);

if Scenario.Run_EFIE
    rmse_itm = rmse(Results.ITM_dB, Results.EFIE_dB);
    bias_itm = bias(Results.ITM_dB, Results.EFIE_dB);
end
```

$$RMSE = \sqrt{\frac{1}{N}\sum_{i=1}^{N}(PL_{model,i} - PL_{EFIE,i})^2}$$

$$Bias = \frac{1}{N}\sum_{i=1}^{N}(PL_{model,i} - PL_{EFIE,i})$$

## H. Computational Considerations

### 1) Complexity Analysis

| Model | Time Complexity | Memory Complexity | Typical Runtime* |
|-------|-----------------|-------------------|------------------|
| EFIE | $O(N^2)$ | $O(N)$ | 30-120 seconds |
| ITM | $O(M \cdot N)$ | $O(N)$ | 1-5 seconds |
| Bullington | $O(M \cdot N)$ | $O(N)$ | < 1 second |
| Deygout | $O(M \cdot N^2)$ | $O(N)$ | 1-3 seconds |
| Hata | $O(M)$ | $O(M)$ | < 0.1 second |

*For 700-meter profile at 970 MHz; $N$ = EFIE segments (~9000), $M$ = profile points (~71)

### 2) Numerical Stability

The implementation includes safeguards for numerical stability:

```matlab
% Minimum distance threshold (1 mm) prevents singularity at source location
% This value is well below typical wavelengths (λ ≈ 0.3m at 970 MHz) and
% discretization scales, ensuring numerical stability without affecting results
dist_2d(dist_2d < 1e-3) = 1e-3;

% Field magnitude floor prevents log(0) in dB conversions
% Value chosen to be negligible compared to physically meaningful field magnitudes
mag_Et(mag_Et < 1e-20) = 1e-20;
mag_Ei(mag_Ei < 1e-20) = 1e-20;
```

### 3) Effective Earth Radius Standardization

To ensure geometric consistency across all models, a shared utility function provides the effective Earth radius:

```matlab
function [a_eff, k_factor] = get_effective_earth_radius(N_s)
    if nargin < 1, N_s = 301; end
    
    a_earth = 6371000;
    k_factor = 1.0 / (1.0 - 0.04665 * exp(N_s / 179.3));
    a_eff = k_factor * a_earth;
end
```

This standardization eliminates geometric discrepancies as a confounding variable in model comparisons.

## I. Summary

This chapter presented the detailed implementation of a computational framework for comparative analysis of radio propagation models. The key implementation contributions include:

1. **Dual-language architecture** separating terrain acquisition (Python/GEE) from propagation modeling (MATLAB), enabling independent development and testing.

2. **EFIE full-wave solver** with Forward-Backward acceleration achieving $O(N^2)$ complexity, providing high-fidelity reference solutions for model validation.

3. **Modified ITM implementation** with continuous diffraction weighting eliminating prediction discontinuities in standard ITM variants.

4. **Standardized model interfaces** enabling uniform invocation and consistent comparative analysis across five propagation prediction approaches.

5. **Geometric consistency** through shared effective Earth radius calculation across all terrain-aware models.

The implementation balances computational efficiency with physical accuracy, enabling practical execution of comparative studies while maintaining rigorous electromagnetic foundation in the EFIE reference model.
