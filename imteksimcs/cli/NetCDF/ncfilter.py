#! /usr/bin/env python

# ncfilter.py
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

# TODO:
# 2019/07/18: Hangs when more tasks than discrete entries in paralelization dim.
#             Porbably issue with MPI output.
'''
Filters a LAMMPS NetCDF trajectory by atom indices from GROMACS sytle .ndx file.

LAMMPS command "group2ndx" creates such a file from groups defined in LAMMPS if
LAMMPS has been compilesd with module "USER-COLVARS".
(refer to https://lammps.sandia.gov/doc/group2ndx.html)

Example: 

  mpirun -n 20 ncfilter.py --log ncfilter.log \\ 
      default.nc filtered.nc groups.ndx indenter

will select all atoms in the 'indenter' group as defined in 'groups.ndx'
by atom ids from 'default.nc' and write a copy of this subset to 'filtered.nc'.
While only error are logged to the terminal, specifying a log file will hold
messages on all verbosity levels.

Requires GromacsWrapper, tested with version 0.8.0 (April 30, 2019)
* https://gromacswrapper.readthedocs.io
* https://github.com/Becksteinlab/GromacsWrapper

Requires MPIFileHandler for logging,
* https://gist.github.com/chengdi123000/42ec8ed2cbef09ee050766c2f25498cb

Source: https://github.com/IMTEK-Simulation/code-snippets
'''
import logging
import argparse

from imteksimcs.NetCDF.ncfilter import ncfilter

def main():

  # in order to have both:
  # * preformatted help text and ...
  # * automatic display of defaults
  class ArgumentDefaultsAndRawDescriptionHelpFormatter(
      argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass

  parser = argparse.ArgumentParser(description=__doc__,
      formatter_class = ArgumentDefaultsAndRawDescriptionHelpFormatter)

  parser.add_argument(
    'infile', default='default.nc', help='NetCDF trajectory')
  parser.add_argument(
    'outfile', default='filtered.nc', help='NetCDF trajectory')
  parser.add_argument(
    'ndx', default='group.ndx', help='GROMACS style .ndx file')
  parser.add_argument(
    '-n', '--ndx', default='group.ndx', dest='ndx',
    help='GROMACS style .ndx file')
  parser.add_argument(
    'group', default='indenter', help='Group selection as named in .ndx')
  parser.add_argument(
    '-g','--group', default='indenter', dest='group',
    help='Group selection as named in .ndx')
  parser.add_argument(
    '-sv','--selection-variable', default='id', dest='selvarname',
    help='Name of NetCDF variable to compare against indices in .ndx file')
  parser.add_argument(
    '-sd','--selection-dimension', default='atom', dest='seldimname',
    help='Name of NetCDF dimension to reduce')
  parser.add_argument(
    '-pd','--parallelization-dimension', default='frame', dest='pardimname',
    help='Name of NetCDF dimension for chunk-wise parallel processing.')
  parser.add_argument(
    '-f','--fmt', '--format', default='NETCDF4', dest='format',
    help=' '.join(('NetCDF format. Refer to',
    'http://unidata.github.io/netcdf4-python/netCDF4/index.html#section1'))
    )
  parser.add_argument(
    '-i','--independent', action='store_false', dest='collective',
    help=' '.join(('Parallel IO mode. Independent mode is used',
      'and unlimited are written as finite dimensions if specified.',
      'Otherwise, collective mode is used by default. Please refert to',
      'http://unidata.github.io/netcdf4-python/netCDF4/index.html#section13'))
    )
  parser.add_argument(
    '-l','--log', default=None, const='log.out', dest='logfile', nargs='?',
    help=' '.join(("Write output to 'log.out' or any other log file whose name",
      "is specified optionally after this flag, instead of stream output to",
      "the terminal"))
    )
  parser.add_argument('--verbose', '-v', action='count', dest='verbose',
      default=0, help='Make this tool more verbose')
  parser.add_argument('--debug','-d', action='store_true',
      help='Make this tool print debug info')
  args = parser.parse_args()

  loglevel = logging.ERROR
  
  if args.verbose > 0:
    loglevel = logging.WARN
  if args.verbose > 1:
    loglevel = logging.INFO
  if args.debug or (args.verbose > 2):
    loglevel = logging.DEBUG

  ncfilter(
    infile      = args.infile,
    outfile      = args.outfile,
    ndx_file    = args.ndx,
    g_sel       = args.group,
    selvarname  = args.selvarname,
    seldimname  = args.seldimname,
    pardimname  = args.pardimname,
    format      = args.format,
    collective  = args.collective,
    loglvl      = loglevel,
    logout      = args.logfile
  )


if __name__ == '__main__':
  main()
