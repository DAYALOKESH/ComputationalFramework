# Unified Radio Propagation Computational Framework

This framework consolidates multiple radio propagation models into a single, modular system. It allows for the comparison of different propagation models (Empirical, Diffraction-based, and Full-Wave) over a given terrain profile.

## Directory Structure

*   **`api/`**: Contains the Master Harness (`main_propagation_comparison.m`) and Data Ingestion tools (`dem_profile_api.py`).
*   **`models/`**: Contains the model implementations.
    *   `itm/`: Longley-Rice (ITM) model.
    *   `empirical/`: Hata and COST 231 models.
    *   `diffraction/`: Bullington and Deygout diffraction models.
    *   `fullwave/`: Electric Field Integral Equation (EFIE) model.
*   **`data/`**: Stores terrain profile data (`terrain_profile.txt`).
*   **`output/`**: Stores generated plots and results.

## Usage

### 1. Setup
Run the setup script in MATLAB to add all necessary folders to the path:
```matlab
setup_project
```

### 2. Generate Data (Optional)
If you need a new terrain profile, use the Python API:
```bash
python api/dem_profile_api.py --start-lat 37.7 --start-lon -122.4 --end-lat 37.8 --end-lon -122.3 --output data/terrain_profile.txt
```

### 3. Run Comparison
Execute the Master Harness in MATLAB:
```matlab
main_propagation_comparison
```
This script will:
1.  Load the terrain profile from `data/terrain_profile.txt`.
2.  Run enabled models (ITM, Hata, Bullington, Deygout, EFIE).
3.  Generate a comparison plot in `output/comparison_plot.png`.
4.  Display RMSE metrics relative to the ITM model.

## Models Included

1.  **Longley-Rice (ITM):** The reference standard for irregular terrain.
2.  **Okumura-Hata / COST 231:** Empirical models for urban/suburban environments.
3.  **Bullington:** Knife-edge diffraction model.
4.  **Deygout:** Dominant obstacle (3-edge) diffraction model.
5.  **EFIE (Full-Wave):** Physics-based forward-backward method (computationally intensive).