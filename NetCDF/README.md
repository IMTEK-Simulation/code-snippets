Code snippets related to NetCDF files.

- `forcesMinMax.nco`: NetCDF operator script to quickly extract 
  minimum and maximum forces with respective atom IDs from a
  LAMMPS NetCDF trajectory.
- `ncsparsify.cpp`: Remove frames from an AMBER NetCDF trajectory file.
- `nc_to_lammps_data.cpp`: Convert an AMBER NetCDF trajectory to a LAMMPS data
  file. Works for very large simulations.
- `selectByAtomType.sh`: Sample bash script using the following two NCO scipts.
- `selectByAtomType.nco`: Slices a sub-set of a LAMMPS NetCDF trajectory 
   by atom type efficiently via NetCDF operators.
- `sortByAtomID.nco`: Sorts along atom dimension by atom IDs
   using NetCDF operatos.
See file headers for usage samples.
