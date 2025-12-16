# III. METHODOLOGY

## A. System Architecture

The computational framework developed for this research implements a hybrid architecture separating data acquisition from propagation modeling. The system consists of two primary modules: a Python-based terrain extraction interface accessing the Google Earth Engine API, and a MATLAB computational engine executing propagation algorithms and comparative analysis. This loose coupling through file-based data exchange enables independent development and testing of each component while maintaining reproducibility across different computing environments.

The terrain extraction module operates as a standalone utility that retrieves elevation profiles along geodesic paths between specified coordinate pairs. The propagation modeling engine accepts these profiles through a standardized text format and orchestrates the execution of multiple prediction models against identical input data. A master harness script (`main_propagation_comparison.m`) coordinates model execution, captures output from each algorithm, and performs statistical analysis using the full-wave EFIE solution as the reference standard.

## B. Terrain Data Acquisition Pipeline

### 1) Data Source Selection

Terrain elevation data was sourced from the NASA Shuttle Radar Topography Mission Global 1 arc-second dataset (SRTMGL1_003) accessed through the Google Earth Engine platform. This dataset provides approximately 30-meter horizontal resolution with vertical accuracy specifications of 16 meters absolute and 6 meters relative (90% confidence). The selection of SRTMGL1_003 balances spatial resolution requirements for sub-wavelength terrain features at UHF frequencies against API throughput constraints and computational tractability.

The extraction pipeline implements automatic fallback logic to alternative datasets (SRTM90, ALOS30, Copernicus GLO30) in cases where the primary dataset returns incomplete coverage. For this study, all test profiles successfully retrieved complete data from the SRTMGL1_003 source without requiring fallback mechanisms.

### 2) Geodesic Path Sampling

Profile extraction between transmitter and receiver coordinates employs great-circle path interpolation to generate uniformly spaced sampling points. The implementation calculates intermediate coordinates through spherical trigonometry, accounting for Earth's curvature over path lengths exceeding several kilometers. For a path defined by start coordinate $(φ_1, λ_1)$ and end coordinate $(φ_2, λ_2)$, the algorithm computes the great-circle distance $d$ and generates $N$ equally spaced points through interpolation parameter $f = i/(N-1)$ where $i ∈ [0, N-1]$.

The interpolation formula applies Vincenty's formulae for inverse and direct geodesic problems on the WGS84 ellipsoid, implemented through the haversine distance calculation and spherical linear interpolation. Batch processing retrieves elevation data in groups of 100 coordinates per API call to the Google Earth Engine `sampleRegions` method, optimizing throughput while maintaining query reliability.

### 3) Output Format Specification

The terrain extraction module writes profile data to ASCII text files with strict formatting: two space-delimited columns representing cumulative distance in meters (Column 1) and terrain elevation in meters above the WGS84 ellipsoid (Column 2). Distance values begin at zero for the transmitter location and increment uniformly to the total path length. This format specification ensures direct compatibility with the MATLAB propagation solvers without requiring data transformation or parsing logic.

## C. Propagation Model Implementations

### 1) Electric Field Integral Equation (EFIE) Reference Model

The full-wave reference implementation solves the two-dimensional Electric Field Integral Equation through the Method of Moments with Forward-Backward acceleration. The discretization strategy samples the terrain profile at spatial intervals of $\Delta x = λ/4$ to satisfy Nyquist sampling requirements for accurate current representation on the surface. This discretization yields $N = \lfloor L_{terrain}/\Delta x \rfloor$ linear segments where $L_{terrain}$ represents the total path length.

**Impedance Matrix Construction:** The mutual impedance between segments $p$ and $q$ separated by distance $\rho_{pq}$ is computed using the zeroth-order Hankel function of the second kind:

$$Z_{pq} = \frac{\beta_0^2}{4\omega\epsilon_0} H_0^{(2)}(\beta_0 \rho_{pq})$$

where $\beta_0 = 2\pi/\lambda$ represents the free-space wavenumber and $H_0^{(2)}$ denotes the Hankel function $H_0^{(2)}(x) = J_0(x) - jY_0(x)$ constructed from Bessel functions of the first and second kind.

Self-impedance calculation requires numerical integration to handle the logarithmic singularity at zero argument. The implementation employs adaptive quadrature integration:

$$Z_{self} = 2Z_{coeff} \int_0^{L/2} [J_0(\beta_0 x) - jY_0(\beta_0 x)] dx$$

where $L$ represents the segment length and $Z_{coeff} = \beta_0^2/(4\omega\epsilon_0)$. The MATLAB `quadgk` function with default tolerance settings performs this integration, automatically handling the endpoint singularity in the $Y_0$ term.

**Forward-Backward Iterative Solution:** The surface current distribution $J$ is solved through a two-sweep iterative process. The forward sweep computes preliminary current values marching from the transmitter location:

$$J_{fwd}(p) = \frac{E_{inc}(p) - \sum_{q=1}^{p-1} J(q) Z_{pq} \ell_q}{Z_{self,p}}$$

where $E_{inc}$ represents the incident field from the line source and $\ell_q$ denotes segment length. The backward sweep refines the solution marching from the receiver location:

$$J(p) = J_{fwd}(p) - \frac{\sum_{q=p+1}^{N} J(q) Z_{pq} \ell_q}{Z_{self,p}}$$

This iterative approach achieves $O(N^2)$ computational complexity compared to $O(N^3)$ for direct matrix inversion, enabling practical solution of profiles with several thousand segments.

**Path Loss Calculation:** The scattered electric field at each receiver position is computed by summing contributions from all surface current elements. The total path loss combines three-dimensional free-space path loss with the two-dimensional excess loss derived from the integral equation solution:

$$PL_{total} = 20\log_{10}\left(\frac{4\pi d}{λ}\right) + 20\log_{10}\left(\frac{|E_{inc,norm}|}{|E_{total,norm}|}\right)$$

where $d$ represents the direct distance between transmitter and receiver, and $|E_{total}| = |E_{inc} - E_{scattered}|$ accounts for the coherent superposition of incident and scattered fields.

The two-dimensional EFIE solution provides field magnitudes that inherently include cylindrical spreading characteristics (1/√ρ dependence where ρ is radial distance from the source). To extract the terrain-specific excess loss—representing diffraction, interference, and scattering effects—these field magnitudes are normalized by √ρ before computing the logarithmic ratio. This normalized excess loss is then combined with three-dimensional free-space path loss to approximate realistic propagation conditions. This approach separates geometric spreading (point source in 3D space) from terrain interaction effects (captured by 2D surface wave analysis), providing a physically consistent path loss estimate.

### 2) Longley-Rice Irregular Terrain Model (ITM)

The Longley-Rice implementation executes the standard point-to-point prediction mode following the algorithms documented in NTIA Technical Report TR-ERL 79-ITS 67. The model requires terrain profile data, terminal heights, frequency, and environmental parameters as input.

**Profile Analysis Subroutine (`qlrpfl`):** Initial profile processing extracts geometric parameters through the `qlrpfl` function. This subroutine computes the effective Earth radius accounting for atmospheric refractivity:

$$a_{eff} = \frac{a_{earth}}{1 - 0.04665 \exp(N_s/179.3)}$$

where $a_{earth} = 6.371 \times 10^6$ meters and $N_s$ represents surface refractivity (set to 301 N-units for standard atmospheric conditions). The horizon search algorithm (`hzns`) identifies radio horizons by finding the maximum grazing angle from each terminal, accounting for Earth curvature through the effective radius.

Terrain irregularity quantification employs the interdecile range parameter $\Delta h$ computed by the `dlthx` function. This routine fits a linear trend to the central 80% of the profile (10th to 90th percentile range), calculates elevation residuals relative to this trend, and extracts the difference between the 90th and 10th percentile residual values. The resulting $\Delta h$ parameter characterizes terrain roughness on a scale from smooth (0-30 meters) to mountainous (>100 meters).

Effective antenna heights $h_e$ are determined through least-squares fitting (`zlsq1`) of the foreground terrain. For each terminal, a linear approximation to the terrain slope is computed over the distance to its respective horizon, and the effective height represents the difference between structural antenna height and the fitted surface elevation at the terminal location.

**Propagation Mode Selection (`lrprop`):** The core prediction engine selects between line-of-sight, diffraction, and troposcatter modes based on path geometry. For distances $d < d_{LOS}$ where $d_{LOS} = d_{L1} + d_{L2}$ represents the sum of horizon distances, the model applies two-ray propagation including ground reflection. The Fresnel reflection coefficient is computed accounting for ground permittivity ($\epsilon_r = 15$), conductivity ($\sigma = 0.005$ S/m), and grazing angle. The resulting attenuation includes interference between direct and reflected rays with appropriate phase relationships.

For distances $d_{LOS} < d < 1.5 d_{LOS}$, the model enters diffraction mode. The implementation blends knife-edge diffraction loss (calculated through Fresnel parameter $\nu$ based on obstacle clearance) with smooth-Earth diffraction (computed via Vogler's method accounting for Earth curvature). The diffraction loss implementation employs a continuous terrain roughness weighting scheme rather than the discrete mode switching in some ITM variants. The weighting function w = (Δh/λ)²/[1 + (Δh/λ)²] provides smooth transition between smooth-earth diffraction (Vogler's method) for flat terrain and knife-edge diffraction for rough terrain. This modification eliminates prediction discontinuities that can occur in standard ITM when terrain roughness crosses threshold values, improving numerical stability for incremental path-by-path analysis along a profile. Additionally, the reference attenuation A_ref is not constrained to non-negative values, preserving physically valid signal enhancements that occur under favorable propagation conditions such as constructive two-ray interference or atmospheric ducting.

Beyond $d > 1.5 d_{LOS}$, tropospheric scatter becomes relevant. The scatter loss calculation accounts for angular distance, scatter angle geometry, frequency, and atmospheric climate conditions. The implementation applies empirical climate-dependent corrections to the basic troposcatter formula, ranging from -5 dB (stable maritime atmosphere) to +5 dB (desert conditions with high atmospheric variability). These corrections supplement the refractivity-based climate modeling in standard ITM. The final attenuation applies the minimum of scatter loss and extrapolated diffraction loss to ensure physical continuity across mode boundaries.

**Variability Calculation (`avar`):** Statistical variability accounts for location, time, and situation uncertainties through independent additive terms. Location variability $Y_L$ scales with terrain irregularity parameter $\Delta h$ and frequency-dependent factor $k_f$. Time variability $Y_T$ depends on climate code (Continental Temperate, climate code 5) and path distance. Situation variability $Y_S$ applies a standard deviation of 5 dB representing antenna pattern variations and local environmental effects. The total variability combines these components:

$$Y_{total} = \text{sign}(z_c) \sqrt{(σ_T z_c)^2 + (σ_L z_c)^2 + (σ_S z_c)^2}$$

where $z_c$ represents the standard normal deviate corresponding to the specified confidence level (50% confidence, $z_c = 0$).

### 3) Bullington Knife-Edge Diffraction

The Bullington method reduces arbitrary terrain profiles to a single equivalent knife-edge obstacle through geometric construction. The algorithm identifies the intersection point of two tangent lines: the steepest slope from the transmitter looking forward, and the steepest slope from the receiver looking backward.

For a profile discretized into $N$ points with coordinates $(x_i, z_i)$, the forward slope calculation determines:

$$m_{tx} = \max_{i=2}^{N-1} \left\{\frac{z_i + b(x_i) - (z_1 + h_{tx})}{x_i}\right\}$$

where $b(x_i) = x_i(D-x_i)/(2a_{eff})$ represents the Earth bulge correction, $h_{tx}$ denotes transmitter structural height, and $D$ represents total path length. Similarly, the backward slope from the receiver computes:

$$m_{rx} = \max_{i=2}^{N-1} \left\{\frac{z_i + b(x_i) - (z_N + h_{rx})}{D - x_i}\right\}$$

The equivalent obstacle position $x_{eq}$ and height $h_{eq}$ are determined by solving the linear system formed by the intersection of these two lines. The clearance parameter $h_{obs}$ represents the vertical distance between the equivalent obstacle and the direct line-of-sight path, accounting for Earth curvature.

The Fresnel diffraction parameter is computed as:

$$\nu = h_{obs} \sqrt{\frac{2(d_1 + d_2)}{λ d_1 d_2}}$$

where $d_1 = x_{eq}$ and $d_2 = D - x_{eq}$ represent distances from the equivalent edge to the transmitter and receiver respectively. The diffraction loss follows the ITU-R P.526 approximation:

$$J(\nu) = \begin{cases} 
0 & \nu \leq -0.78 \\
6.9 + 20\log_{10}(\sqrt{(\nu-0.1)^2 + 1} + \nu - 0.1) & \nu > -0.78
\end{cases}$$

Total path loss combines free-space path loss with the diffraction loss $J(\nu)$.

### 4) Deygout Multiple-Edge Diffraction

The Deygout implementation applies a recursive three-edge approach to handle multiple terrain obstacles. The algorithm first identifies the principal edge—the point with maximum Fresnel parameter $\nu_p$ across the entire path. This principal obstacle experiences the greatest diffraction effect and dominates the propagation behavior.

After computing the principal edge loss $J(\nu_p)$, the algorithm recursively analyzes two sub-paths: transmitter to principal edge, and principal edge to receiver. For each sub-path, the process repeats to identify secondary edges with Fresnel parameters $\nu_t$ and $\nu_r$ respectively. The recursion terminates at depth two (three total edges: principal plus two secondary), consistent with the model's original formulation.

The total diffraction loss combines these edge contributions with empirical correction factors:

$$L_{total} = J(\nu_p) + T \cdot [J(\nu_t) + J(\nu_r) + C]$$

where $T$ represents a transition factor:

$$T = \begin{cases}
J(\nu_p)/6 & J(\nu_p) < 6 \text{ dB} \\
1 & J(\nu_p) \geq 6 \text{ dB}
\end{cases}$$

and $C = 8.0 + 0.04 D_{km}$ provides distance-dependent correction with $D_{km}$ representing total path length in kilometers. This heuristic correction accounts for multiple scattering effects and ensures smooth transition behavior as additional edges become significant.

The effective Earth radius used in all diffraction calculations is computed using the ITM standard formula:

$$k = \frac{1}{1 - 0.04665 \cdot \exp(N_s/179.3)}$$

$$a_{eff} = k \cdot a_{earth}$$

For standard atmosphere ($N_s = 301$ N-units), this yields $k \approx 4/3$ and $a_{eff} \approx 8.495 \times 10^6$ meters.

### 5) Okumura-Hata and COST-231 Empirical Models

The empirical model implementation applies frequency-dependent selection logic to choose between the original Okumura-Hata formulation (valid 150-1500 MHz) and the COST-231 extension (valid 1500-2000 MHz). For the study frequency of 970 MHz, the system applies the standard Okumura-Hata equations.

The base path loss for urban environments is computed as:

$$L_U = 69.55 + 26.16\log_{10}(f) - 13.82\log_{10}(h_B) - a(h_M) + [44.9 - 6.55\log_{10}(h_B)]\log_{10}(d)$$

where $f$ represents frequency in MHz, $h_B$ represents base station height in meters, $d$ represents distance in kilometers, and $a(h_M)$ represents the mobile antenna height correction factor. For small to medium-sized cities:

$$a(h_M) = [1.1\log_{10}(f) - 0.7]h_M - [1.56\log_{10}(f) - 0.8]$$

Suburban and rural corrections subtract frequency-dependent terms from the urban baseline. The suburban correction applies:

$$L_{suburban} = L_U - 2[\log_{10}(f/28)]^2 - 5.4$$

The rural (open area) correction applies:

$$L_{rural} = L_U - 4.78[\log_{10}(f)]^2 + 18.33\log_{10}(f) - 40.94$$

The implementation treats the terrain profile as a distance vector and computes path loss as a function of distance from the transmitter, independent of specific terrain elevation variations. Environmental classification is specified through the `Scenario.Environment` parameter with options for urban, suburban, and rural conditions.

The Hata model's terrain-independent formulation presents a fundamental limitation for this comparative study. While EFIE, ITM, Bullington, and Deygout explicitly account for terrain profile variations, the Hata model applies a statistical average behavior derived from measurements in specific urban morphologies. Consequently, Hata predictions represent expected median loss over an ensemble of similar paths rather than site-specific predictions for the actual terrain profile. This distinction must be considered when interpreting comparative metrics: deviations between Hata and terrain-specific models reflect both model accuracy and the fundamental difference in prediction philosophy (statistical ensemble vs. deterministic path-specific).

## D. Simulation Parameters and Configuration

All simulations executed under identical parameter settings to ensure valid comparison across models. The carrier frequency was set to 970 MHz, corresponding to the lower portion of the UHF band commonly used for broadcasting and mobile communications. Transmitter antenna height was specified as 52.0 meters above ground level, representative of broadcast tower installations. Receiver antenna height was set to 2.4 meters above ground level, consistent with mobile terminal or fixed-receiver configurations.

For the Longley-Rice model, ground electrical parameters were set to dielectric constant $\epsilon_r = 15$ and conductivity $\sigma = 0.005$ S/m, representing average terrain conditions. Surface refractivity was specified as $N_s = 301$ N-units (continental temperate climate). Polarization was set to horizontal, and confidence level was maintained at 50% (median predictions).

The Hata/COST-231 model was configured for suburban environment classification, applying the appropriate correction factors to the urban baseline predictions. This classification was selected to match typical propagation conditions in the test scenarios without requiring site-specific calibration.

Geometric consistency across all propagation models was ensured through standardized calculation of the effective Earth radius. All models employ the ITM standard formula k = 1/(1 - 0.04665·exp(Ns/179.3)) with surface refractivity Ns = 301 N-units, yielding k ≈ 4/3 and effective radius aeff ≈ 8.495×10⁶ meters. This standardization ensures that diffraction calculations, horizon determinations, and Earth curvature corrections are consistent across the EFIE, ITM, Bullington, and Deygout implementations, eliminating geometric discrepancies as a confounding variable in model comparisons.

## E. Comparative Analysis Framework

### 1) Statistical Metrics

Model accuracy evaluation employs two primary statistical measures computed against the EFIE full-wave solution as ground truth. Root Mean Square Error quantifies the standard deviation of prediction errors:

$$RMSE = \sqrt{\frac{1}{N}\sum_{i=1}^{N}(PL_{model,i} - PL_{EFIE,i})^2}$$

Mean Bias Error measures systematic tendency toward over- or under-prediction:

$$Bias = \frac{1}{N}\sum_{i=1}^{N}(PL_{model,i} - PL_{EFIE,i})$$

where $PL_{model,i}$ represents the predicted path loss from the model under evaluation, $PL_{EFIE,i}$ represents the full-wave reference solution, and $N$ denotes the number of sample points along the profile. Positive bias indicates conservative (pessimistic) predictions with higher estimated path loss than the reference, while negative bias indicates optimistic predictions underestimating actual propagation loss.

The selection of EFIE as the reference standard warrants clarification regarding its limitations. The EFIE implementation employs a two-dimensional integral equation solver with perfect electric conductor (PEC) boundary conditions, providing high-fidelity modeling of diffraction and surface wave interactions. However, the conversion from 2D field solutions to 3D path loss involves approximations in separating geometric spreading from terrain-specific effects. Additionally, the PEC assumption overestimates ground interaction effects compared to realistic lossy soil. These limitations imply that EFIE should be considered a high-fidelity reference for relative model comparison rather than absolute ground truth. The RMSE and bias metrics therefore quantify heuristic model deviation from a consistent computational reference, not necessarily from measured propagation behavior.

### 2) Computational Implementation

The master harness script loads terrain profile data, validates input consistency (removing any NaN values or discontinuities), and sequentially executes each propagation model through standardized wrapper functions. Each wrapper accepts the common `Scenario` structure (containing frequency, antenna heights, and environment classification) and the `TerrainProfile` structure (containing distance and elevation vectors) as inputs, returning a path loss vector with consistent dimensionality.

Model execution order proceeds from computationally intensive to efficient: EFIE reference (if enabled), Longley-Rice ITM, empirical Hata model, Bullington diffraction, and Deygout diffraction. Timing measurements capture execution duration for each model using MATLAB's `tic`/`toc` functions, though these measurements are not included in the formal analysis due to hardware-dependent variability.

Results aggregation constructs a unified data structure containing all model outputs aligned by distance coordinate. The statistical analysis module computes RMSE and bias metrics for each heuristic model relative to the EFIE reference, conditional on EFIE execution being enabled. Output generation produces comparison plots showing all models on common axes, terrain geometry visualization with transmitter/receiver locations, and a CSV file containing metric summary data for external analysis.

### 3) Reproducibility Considerations

The implementation maintains strict separation between configuration parameters and algorithm code to enable reproducible execution across different scenarios. All physical constants (speed of light, wavelength calculations, impedance relations) are computed from fundamental values rather than hardcoded approximations. File paths employ platform-independent construction through MATLAB's `fullfile` function. The setup script (`setup_project.m`) automatically configures the MATLAB path to include all necessary subdirectories, eliminating manual path management and reducing configuration errors.

Version control through the project structure enables tracking of algorithm modifications and parameter adjustments. The modular architecture permits independent testing of individual propagation models through standalone runner scripts before integration into the comparative framework. This design facilitates debugging, validation against reference implementations, and incremental development of model improvements.

## F. Model Limitations and Assumptions

Each propagation model implementation embodies specific simplifying assumptions that constrain its domain of validity:

**EFIE Full-Wave Model**: Assumes two-dimensional geometry with infinite extent perpendicular to the propagation plane, perfect electric conductor ground, and combines 2D diffraction effects with 3D spreading approximations. Computational cost scales as O(N²) with profile discretization, limiting practical application to paths under several kilometers with λ/4 sampling.

**Longley-Rice ITM**: Employs median statistical predictions incorporating terrain irregularity as a bulk parameter rather than resolving individual terrain features. The implementation includes a modified continuous diffraction weighting scheme and preserves negative reference attenuation values representing signal enhancement, deviating slightly from some ITM variants. Valid for frequencies 20 MHz to 20 GHz and path lengths 1 to 2000 km.

**Bullington Diffraction**: Reduces arbitrary terrain to a single equivalent knife-edge, potentially underestimating loss when multiple comparable obstacles exist. Most accurate for profiles with one dominant ridge. Does not account for atmospheric effects or frequency-dependent ground interactions beyond the wavelength dependence in the Fresnel parameter.

**Deygout Diffraction**: Limits multiple-edge treatment to three obstacles (principal plus two secondary edges) with empirical correction factors (C = 8.0 + 0.04D_km) derived from limited measurement campaigns. Accuracy degrades for profiles with many comparable obstacles where higher-order scattering becomes significant.

**Okumura-Hata/COST-231**: Applies statistical median predictions derived from measurements in specific urban morphologies (Tokyo metropolitan area for original Hata formulation). Completely independent of actual terrain profile geometry, providing ensemble-average behavior rather than site-specific predictions. Valid only within specified parameter ranges (frequency, antenna heights, distances) and environmental classifications (urban, suburban, rural) that match the original measurement campaigns.

These limitations contextualize the comparative analysis: discrepancies between models reflect both prediction accuracy and fundamental differences in modeling philosophy (deterministic vs. statistical, terrain-specific vs. ensemble-average, full-wave vs. high-frequency approximation).
