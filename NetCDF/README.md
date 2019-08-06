Code snippets related to NetCDF files.

- `ncsparsify.cpp`: Remove frames from an AMBER NetCDF trajectory file.
- `nc_to_lammps_data.cpp`: Convert an AMBER NetCDF trajectory to a LAMMPS data file. Works for very large simulations.
- `ncfilter.py`: Extracts subsets (of atoms, specified by GROMACS-type .ndx files) 
  from AMBER NetCDF trajectories in parallel efficiently via MPI-IO-enabled NetCDF library.
  Requires `mpi4py` and `netCDF4`, as well as `GromacsWrapper` for .ndx file input and `MPIFileHandler.py` for MPI-IO logging.
  Tested on NEMO with 
  - `devel/python/3.6.5` (global module) 
  - `mpi4py/3.0.0-python-3.6.5-openmpi-3.1-gnu-7.3` (group-internal module)
  - `netcdf4-python/1.5.0` (group-internal module) 
- `ncjoin.py`: Merge multiple AMBER NetCDF trajectories into a single AMBER NetCDF trajectory.

See file headers for usage samples.
