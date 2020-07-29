#!/usr/bin/env python
#
# Copyright (C) 2019 IMTEK Simulation
# Author: Johannes Hoermann, johannes.hoermann@imtek.uni-freiburg.de
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
""" Extracts properties from LAMMPS data files and NetCDF trajectories

    8 Apr 2019, extracts properties from LAMMPS data files and NetCDF
    Johannes Hörmann, johannes.hoermann@imtek.uni-freiburg.de
    tested with Ovito 3.0.0-dev234 """

import sys, logging
import numpy as np

from imteksimcs.LAMMPS_data.lmp_extract_property import extract_rdf, extract_box_measures, extract_nothing
from imteksimcs.LAMMPS_data.lmp_extract_property import select_fcc_only

# property keyword - function mapping:
property_function = {
  'rdf':    extract_rdf,
  'box':    extract_box_measures,
  'none':   extract_nothing }

# modfier keywor - function mapping:
modifier_function = {
  'fcc':    select_fcc_only,
  'none':   None }

logger = logging.getLogger(__name__)
logfmt = "[%(levelname)s - %(filename)s:%(lineno)s - %(funcName)s() ] %(message)s (%(asctime)s)"
logging.basicConfig( format = logfmt )


def parse_property_function(s):
  logger.info("Looking up property '{:s}'.".format(s))
  try:
    return property_function[s]
  except:
    logger.exception("No extraction method exists for property '{:s}'!".format(s))
    return None


def parse_modifier_function(s):
  logger.info("Looking up modifier '{:s}'.".format(s))
  try:
    return modifier_function[s]
  except:
    logger.exception("No modifying method exists for modfifier '{:s}'!".format(s))
    return None


def parse_ranges(s):
    """ Parses ranges specified in a string into lists:
    Examples
    --------
    In [3]: s = '1, 2, 3, 4, 5, 6'
    In [4]: parseRanges(s)
    Out[4]: [1, 2, 3, 4, 5, 6]
    In [5]: parseRanges('3-6')
    Out[5]: [3, 4, 5]
    In [6]: parseRanges('3-16-3')
    Out[6]: [3, 6, 9, 12, 15]
    In [6]: parseRanges('1,20-22,3-16-3')
    Out[6]: [1, 3, 6, 9, 12, 15, 20, 21, 22]
    Adapted from
        https://stackoverflow.com/questions/4726168/parsing-command-line-input-for-numbers
    """
    result = set()
    for part in s.split(','):
        #x = part.split('-')
        if '-' in part:
            i = part.index('-')
            result.update( range(*map(int, part.split('-'))) )
        else:
            result.add( int(part) )
        #result.update(range(int(x[0]), int(x[-1]) + 1))
    return sorted(result)


def main():
  import argparse  

  parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter )
  parser.add_argument('--property', '-p',
    help="Properties to extract, one or several of {}".format(
      list(property_function.keys())),
    type      = parse_property_function,
    default   = "none",
    nargs     = '+',
    metavar   = "PROP",
    required  = True,
    choices   = list(property_function.values()))
  parser.add_argument('--modifier', '-m',
    help="Modifiers to apply before extraction, applied in order specified. Choose one or several from {}".format(
      list(modifier_function.keys())),
    type      = parse_modifier_function,
    default   = None,
    nargs     = '*',
    metavar   = "MOD",
    choices   = list(modifier_function.values()))
  parser.add_argument('--frames', '-f',
    help = "Frame selection, can be list or range of format '1,3,5-7'",
    type      = parse_ranges,
    default   = None)
  parser.add_argument('--style', '-s',
    help="LAMMPS atom style (default: full)",
    default   = "full")
  parser.add_argument('--trajectory', help='NetCDF trajectory file',
    metavar   = 'NETCDF', 
    default   = None)
  parser.add_argument('--verbose', '-v', action='store_true',
    help='Make this tool more verbose')
  parser.add_argument('--debug', action='store_true',
    help='Make this tool print debug info')
  parser.add_argument('topology_file', help='LAMMPS data file',
    metavar   = 'DATAFILE')
  parser.add_argument('outfile', 
    help='(List of) text files to store extracted properties',
    nargs     = '+', 
    metavar   = 'OUT')
  args = parser.parse_args()

  if args.debug:
    loglevel = logging.DEBUG
  elif args.verbose:
    loglevel = logging.INFO
  else:
    loglevel = logging.WARNING
  
  logger.setLevel(loglevel)

  if args.frames is not None:
    logger.info("Selected frames '{}'".format(args.frames))

  logger.info("Expecting atom style '{}'".format(args.style))

  for i, property in enumerate(args.property):
    if property is not None:
      try:
        outfile = args.outfile[i]
      except:
        logger.exception("No outfile specified for property function {:s}()!".format(
          property.__name__))
        return 1

      logger.info("Evaluating with {:s}()...".format(property.__name__))

      try:
        extractProperty( 
          topology_file   = args.topology_file,
          trajectory_file = args.trajectory,
          frames          = args.frames,
          property_func   = property,
          modifier_func   = args.modifier,
          outfile         = outfile,
          atom_style      = args.style)
      except:
        logger.exception("Unhandled exception while processing {:s}()!".format(
          property.__name__))          

  return 0


if __name__ == '__main__':
    main()
