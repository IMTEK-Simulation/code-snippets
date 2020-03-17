Code snippets related to NetCDF files.

- `forcesMinMax.nco`: NetCDF operator script to quickly extract 
  minimum and maximum forces with respective atom IDs from a
  LAMMPS NetCDF trajectory.
- `ncsparsify.cpp`: Remove frames from an AMBER NetCDF trajectory file.
- `nc_to_lammps_data.cpp`: Convert an AMBER NetCDF trajectory to a LAMMPS data
  file. Works for very large simulations.
- `ncfilter.py`: Extracts subsets (of atoms, specified by GROMACS-type .ndx
  files) from AMBER NetCDF trajectories in parallel efficiently via 
  MPI-IO-enabled NetCDF library. Requires `mpi4py` and `netCDF4`, as well as 
  `GromacsWrapper` for .ndx file input and `MPIFileHandler.py` 
  for MPI-IO logging. Tested on NEMO with
  - `devel/python/3.6.5` (global module)
  - `mpi4py/3.0.0-python-3.6.5-openmpi-3.1-gnu-7.3` (group-internal module)
  - `netcdf4-python/1.5.0` (group-internal module)
- `ncjoin.py`: Merge multiple AMBER NetCDF trajectories into a single AMBER
  NetCDF trajectory.
- `netcdf2data.py`: Same purpose as `nc_to_lammps_data.cpp` above, but
  `ovitos`-based.
- `selectByAtomType.sh`: Sample bash script using the following two NCO scipts.
- `selectByAtomType.nco`: Slices a sub-set of a LAMMPS NetCDF trajectory 
   by atom type efficiently via NetCDF operators.
- `sortByAtomID.nco`: Sorts along atom dimension by atom IDs
   using NetCDF operatos.
See file headers for usage samples.
