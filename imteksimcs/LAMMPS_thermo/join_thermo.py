#!/usr/bin/env python
"""Concatenates two LAMMPS thermo.out-like white-space delimited tables"""
# 25 Feb 2019, jlh
import logging


def joinThermo(file1, file2, outfile=None, index_col=0, hashed_header=False,
  skiprows=0, delimiter='\s+'):
  """Joins (concatenates) two LAMMPS thermo output - like tables by comparing
     their indices, keeping rows from the table with a higher minimum index
     and discarding rows from the table with lower minimum index in the case
     of rows with overlapping indices.

    Parameters
    ----------
    file1, file2 : str,
        LAMMPS thermo.out file (thermo section from LAMMPS log file),
        white-space delimited with single-line column names header.
    outfile:       str, optional
        Output file name. Default: None
    index_col:     int, optional
        Position of column to use as index for table.
        Default: 0 (first column)
    """
  import pandas as pd
  logger = logging.getLogger('join_thermo.joinThermo')
  logger.info("Reading files '{:s}' and '{:s}'".format(file1,file2))

  if hashed_header:
    # Take original header from first file, ignore second (assumed same)
    part1, original_header = read_data_with_hashed_header(file1, skiprows,
                               index_col=index_col, delimiter=delimiter)
    part2, _               = read_data_with_hashed_header(file2, skiprows,
                               index_col=index_col, delimiter=delimiter)
  else:
    part1 = pd.read_csv(file1, sep=delimiter,
      index_col=index_col, skiprows=skiprows)
    part2 = pd.read_csv(file2, sep=delimiter,
      index_col=index_col, skiprows=skiprows)

  logger.info( "'{:s}' contains {:d} rows, '{:s}' {:d}.".format( file1,
    len(part1), file2, len(part2) ) )
  logger.info("Column '{:s}' designated as index.".format(part1.index.name))

  if min(part1.index) > min(part2.index):
    logger.info("{:s} predates {:s}.".format(file2,file1))
    tmp = part1
    part1 = part2
    part2 = tmp
  else:
    logger.info("{:s} predates {:s}.".format(file1,file2))

  first_overlapping_index = part2.index[0]
  if max(part1.index) < min(part2.index):
    logger.warn("No overlap between input files!")
  else:
    logger.info("First overlapping index is '{:s}': '{}'".format(
      part1.index.name, first_overlapping_index ) )

  # concatenate
  selection = part1.index < min(part2.index)
  joint = part1.loc[selection].append(part2)
  logger.info( "Joint table contains {:d} rows.".format( len(joint) ) )

  if outfile:
    logger.info( "Writing concatenated table to '{:s}'".format(outfile) )
    write_header_via_pandas = True
    with open(outfile, 'w') as f:
      if hashed_header and len(original_header) > 0:
        write_header_via_pandas = False # header assumed to be written now
        for header_line in original_header:
          f.write(header_line)
      joint.to_csv(f, sep=' ', header=write_header_via_pandas, index=True, float_format='%g')

  return joint


# for a LAMMPS thermo ave output, the file header may look like this:
#
#  # Time-averaged data for fix thermo_average
#  # TimeStep c_thermo_pe c_thermo_temp c_thermo_press [...]
#  1000 -1.86167e+06 293.853 -3.1683  [...]
def read_data_with_hashed_header(filename,skiprows=1,
  index_col=0, delimiter='\s+'):
  """Reads hashed line as header containing column names"""
  import pandas as pd
  logger = logging.getLogger('join_thermo.read_data_with_hashed_header')

  linecount = len(open(filename, 'r').readlines())
  logger.info("File '{:s}' has {:d} lines, {:d} of which are skipped.".format(
    filename, linecount, skiprows))

  header = pd.read_csv(filename,
    delimiter=delimiter, skiprows=skiprows,nrows=0)
  columns = header.columns[1:] # skip initial hash
  df = pd.read_csv( filename, delimiter=delimiter, header=None, comment='#',
    skiprows=skiprows, names=columns, index_col=index_col )

  header_linecount = linecount - len(df)
  original_header = []
  with open(filename, 'r') as f:
    for i in range(header_linecount):
      original_header.append( f.readline() )

  logger.info("Original header has {:d} lines.".format(header_linecount))

  return df, original_header
