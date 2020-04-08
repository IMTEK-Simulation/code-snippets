Python modules related to LAMMPS data text files.

- `MergeLammpsDatafiles.py`: LAMMPS datafiles extracted from
  NetCDF lack topology information. Many tools such  as `ovitos` or VMD's
  TopoTools drop topology information when writing LAMMPS datafiles.
  `MergeLammpsDatafiles.py` and its command line wrapper `merge.py` utilize
  `pizza.py` to compare a possibly incomplete datafile against a reference and
  restore missing sections.
- `strip_comments.py`: strips all comments from LAMMPS data file (for processing with pizza.py)
- `to_hybrid.py`: Converts LAMMPS pair_coeff lines to 'pair_style hybrid' format.

See file headers for usage samples.
