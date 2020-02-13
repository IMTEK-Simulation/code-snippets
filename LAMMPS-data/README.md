Code snippets related to LAMMPS data text files.

- `merge.py`: command line interface for `MergeLammpsDatafiles.py`, see below.
- `MergeLammpsDatafiles.py`: LAMMPS datafiles extracted from
  NetCDF lack topology information. Many tools such  as `ovitos` or VMD's
  TopoTools drop topology information when writing LAMMPS datafiles.
  `MergeLammpsDatafiles.py` and its command line wrapper `merge.py` utilize
  `pizza.py` to compare a possibly incomplete datafile against a reference and
  restore missing sections.

See file headers for usage samples.
