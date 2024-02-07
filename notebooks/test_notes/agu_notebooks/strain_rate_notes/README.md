# Description of files in `strain_rate_notes`
This file details all of the notebooks in the folder `agu_notebooks`. All notebooks and subdirectories were created by George Lu in 2023 to analyze ApRES data from the Greenland Lakes Project 2022-2023 field campaign. 

- `ronne_test`: Testing the displacement calculation on data from Ronne ice shelf. 
- `seasonal_coherence`: Examines how coherence, and thus error, varies with depth and time. Compares to temperature. 
- `strain_rate_2burst`: Compares summing up displacements with directly calculating displacement between bursts further apart. 
- `strain_rate_drainage`: Initial notebook that looks at how the strain rates evlove during the lake drainage
- `vels_compare`: Comparing displacement calculations attained with Irena Vankova's MATLAB scripts with the xapres methods. 
- `winter_strain_rates`: An initial effort to determine strain rate estimates by fitting a line to calculated vertical velocities. 

The folder `winter_strain_rate_tests` contains additional notebooks that examine the sensitivity of calculated strain rates to various parameters. They are organized by station and by parameter that is varied (how many stacks, the coherence calculation window, how many displacements are summed), and they follow the same template, described by `template`. 