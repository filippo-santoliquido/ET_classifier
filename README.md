## Classifying Binary Black Holes from Population III Stars with the Einstein Telescope: A Machine-Learning Approach

This repository includes the code and resources to replicate all results and figures from [Santoliquido et al., 2024](https://ui.adsabs.harvard.edu/abs/2024A%26A...690A.362S/abstract).

### Setup Instructions

To use the notebooks, first download the data from this [Zenodo repository](https://zenodo.org/records/14013888). This dataset contains parameter estimations generated with [GWFish](https://github.com/janosch314/GWFish) of black hole mergers from Population III and Population I-II stars, based on models in [Iorio et al. 2023](https://ui.adsabs.harvard.edu/abs/2023MNRAS.524..426I/abstract) and [Santoliquido et al., 2023](https://ui.adsabs.harvard.edu/abs/2023MNRAS.524..307S/abstract).

### Notebooks Overview

1. `Figure1_relative_errors.ipynb`  
   - This notebook uses `GWFish` catalogs to recreate Figure 1 from [Santoliquido et al., 2024](https://ui.adsabs.harvard.edu/abs/2024A%26A...690A.362S/abstract).

2. `ResamplingWithPriors.ipynb`  
   - This notebook reads the `GWFish` catalogs and applies resampling based on prior ranges. Before running it, create the following folders in `clf_ET`:
     - `measured_values_v1`
     - `test_catalogs_v1`
     - `test_catalogs_v1/CornerPlotsToInspectSampling` (corner plots are generated only for the pessimistic 10-year scenario).
   - **Note**: This notebook takes time to execute and needs to be run for all configurations listed within it.

3. `XGBoost_ET.ipynb` 
   - The primary notebook of the project. Run this to reproduce all main results.

For any questions or further assistance, feel free to reach out to me at [filippo.santoliquido@gssi.it](mailto:filippo.santoliquido@gssi.it).
