# Test Notebooks
This folder contains a series of notebooks used to test the scripts.
The following list will contain the notebook names, along with what they're testing:
- `error_deep_dive`: plots average error at select depths along with temperature and voltage for all data
- `interpolation`: explores impact of interpreting data onto regular grid for time
- `lake_drainage`: examines lake drainages across stations using new 'leapfrog' method
- `method_comparison`: compares resulting strain rates determined using summation method or 'leapfrog' method
- `strain_rate_method_temp`: applying leapfrog method to different lag windows, across all data
- `strain_rate_method`: development of the new 'leapfrog' method to calculate strain rate
- `toy_strain_rate`: a demonstration of how the two different methods leads to different results

Test notes from development by Jonny Kingslake
- `testing_xapres`: testing loading in .DAT files into Xarray and plotting some data
- `tests_for_github`: tests file selection and searching for implementation in Github tests
- `unattended_testing`: developing and testing loading and structuring ApRES data collected in attended mode. 
