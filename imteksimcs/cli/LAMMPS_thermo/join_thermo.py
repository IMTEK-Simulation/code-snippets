#!/usr/bin/env python
"""Concatenates two LAMMPS thermo.out-like white-space delimited tables"""
# 25 Feb 2019, jlh
import logging

from imteksimcs.LAMMPS_thermo.join_thermo import joinThermo


def main():
  import argparse
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument('--index', '-i',
    help="Zero-indexed index column, default: 0 (first column)",
    default=0)
  parser.add_argument('--hashed-header', '-hh', action='store_true',
    help="Indicate that header of table is preceded by has ('#').")
  parser.add_argument('--skiprows', '-s',
    help="Rows to skip before header, default: 0, but 1 for -hh",
    default=None)
  parser.add_argument('--verbose', '-v', action='store_true',
    help='Make this tool more verbose')
  parser.add_argument('--debug', action='store_true',
    help='Make this tool print debug info')
  parser.add_argument('first_input_file', help='''First input file (usually some
    white-space delimited LAMMPS thermo.out log file)''')
  parser.add_argument('second_input_file', help='Second input file')
  parser.add_argument('output_file', help='Name of concatenated output file',
    default='joint.out')
  args = parser.parse_args()

  if args.debug:
    loglevel = logging.DEBUG
  elif args.verbose:
    loglevel = logging.INFO
  else:
    loglevel = logging.WARNING

  logging.basicConfig(level = loglevel)
  logger = logging.getLogger('join_thermo.main')

  if args.hashed_header and args.skiprows is None:
    args.skiprows = 1
  elif args.hashed_header is None:
    args.skiprows = 0

  joinThermo( args.first_input_file, args.second_input_file, args.output_file,
    index_col=args.index, hashed_header=args.hashed_header, skiprows=args.skiprows )

  return


if __name__ == '__main__':
    main()
