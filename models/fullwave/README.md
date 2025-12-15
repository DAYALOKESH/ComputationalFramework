# EFIE_FBM

## Electric Field Integral Equation with Forward-Backward Method

This repository implements the Electric Field Integral Equation (EFIE) using the Forward-Backward Method (FBM) to calculate electromagnetic field scattering from an undulating terrain surface excited by an infinite line source.

### Overview

The implementation solves the EFIE on terrain data to compute:
- **Surface Current Distribution**: The induced current on the terrain surface
- **Electric Field**: The total electric field at observation points above the terrain

Based on Balanis EM Theory (Page 680) and the Forward-Backward iterative Method.

### Files

- `efie_fbm.m` - Main MATLAB/Octave script implementing the EFIE-FBM algorithm
- `X.txt` - Terrain profile data (x-coordinate, height)
- `efie.txt` - Reference C++ implementation

### Requirements

- MATLAB R2016b or later, OR
- GNU Octave 4.0 or later

### Usage

1. Ensure the terrain file `X.txt` is in the same directory as `efie_fbm.m`
2. Run the script in MATLAB or Octave:

```matlab
% In MATLAB
run('efie_fbm.m')

% In Octave (command line)
octave --no-gui efie_fbm.m
```

### Configuration

The discretization factor can be adjusted in `efie_fbm.m` to trade off between accuracy and computation time:

```matlab
% Recommended values:
%   4.0 = lambda/4 (highest accuracy, slowest)
%   2.0 = lambda/2 (good accuracy, balanced speed) [DEFAULT]
%   1.0 = lambda (moderate accuracy, fast)
%   0.5 = 2*lambda (lower accuracy, very fast - for testing)
discretization_factor = 2.0;
```

### Output

The script generates:

#### Data Files
- `J_current.dat` - Surface current magnitude vs distance
- `E_field.dat` - Electric field (dB) vs distance
- `E_field_linear.dat` - Electric field (linear) vs distance
- `performance_metrics.txt` - Computation time and memory usage

#### Plots
- `EFIE_FBM_Results.png` - Four-panel plot with terrain profile, surface current, and electric field
- `EFIE_FBM_Terrain_Overlay.png` - Surface current and electric field overlaid on terrain

### Performance Metrics

The script measures and reports:
- **Computation Time**: Time taken for all calculations (excludes plotting)
- **Memory Usage**: Memory consumed during computation

### Physical Parameters

Default settings (configurable in script):
- Operating Frequency: 970 MHz
- Source Location: (0, 442) meters
- Observation Height: 2.4 meters above terrain

### Algorithm

The Forward-Backward Method (FBM) is an iterative approach that:
1. **Forward Sweep**: Solves for surface currents from source to far end
2. **Backward Sweep**: Refines the solution from far end back to source
3. **Field Calculation**: Computes total electric field at observation points

This approach is more memory-efficient than direct matrix inversion methods for large problems.

### License

This project is for educational and research purposes.