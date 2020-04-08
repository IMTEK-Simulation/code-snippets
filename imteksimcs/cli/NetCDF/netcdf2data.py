#!/usr/bin/env ovitos
#
# netcdf2data.py
#
# Copyright (C) 2018, 2019 IMTEK Simulation
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

""" Extracts selected frames from NetCDF to LAMMPS data files """
# 22 Oct 2018, ovitos script for extracting NetCDF frames as LAMMPS data files

import argparse
import logging

from imteksimcs.NetCDF.netcdf2data import extractFrames
from imteksimcs.util.parse_ranges import parseRanges

# adapted from
# https://ovito.org/manual/python/modules/ovito_modifiers.html#ovito.modifiers.LoadTrajectoryModifier
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--frames', '-f',
        help = "Frames to extract, can be list or range of format '1,3,5-7'",
        type = parseRanges,
        default = None)
    parser.add_argument('--style', '-s',
        help="LAMMPS atom style (default: full)",
        default="full")
    parser.add_argument('--verbose', '-v', action='store_true',
        help='Make this tool more verbose')
    parser.add_argument('--debug', action='store_true',
        help='Make this tool print debug info')
    parser.add_argument('topology_file', help='LAMMPS data file')
    parser.add_argument('trajectory_file', help='NetCDF trajectory file')
    args = parser.parse_args()

    if args.debug:
        loglevel = logging.DEBUG
    elif args.verbose:
        loglevel = logging.INFO
    else:
        loglevel = logging.WARNING

    logging.basicConfig(level = loglevel)
    logger = logging.getLogger('netcdf2data.main')

    if args.frames is not None:
        logger.info("Extracting selected frames '{}'".format(args.frames))

    logger.info("Expecting atom style '{}'".format(args.style))

    extractFrames(
        topology_file   = args.topology_file,
        trajectory_file = args.trajectory_file,
        frames          = args.frames,
        atom_style      = args.style)

    return


if __name__ == '__main__':
  main()
