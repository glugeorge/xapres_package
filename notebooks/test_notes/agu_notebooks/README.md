# Description of files in `agu_notebooks`
This file details all of the notebooks in the folder `agu_notebooks`. All notebooks and subdirectories were created by George Lu in 2023 to analyze ApRES data from the Greenland Lakes Project 2022-2023 field campaign. 

- `agu_figures`: Generates the figures used in George Lu's AGU presentation in 2023.
- `clipping_check`: Plots histograms of chirps to look for clipping.
- `correcting_clipped_data`: Investigates impact of applying half-chirp window on data that is clipped.
- `debug_filename`: Recreates bug where filename was cut off.
- `processing_steps`: Step-by-step description of how to obtain the complex profiles
- `profile_sensitivity`: Examines sensitivity of profiles to different attenuator settings and bandwidths
- `strain_rate_steps_new`: Step-by-step description of attaining strain rates using xarray's resample feature to sum up unit displacements to a regular grid
- `winter_data_check`: Looks at station data from winter 2022-2023, and investigates datagaps at station A14

The folder `strain_rate_notes` contains additional notebooks aiming to understand calculated strain rates.