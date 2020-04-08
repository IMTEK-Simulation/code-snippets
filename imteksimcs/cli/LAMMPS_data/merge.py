#!/bin/sh
#
# mergy.py
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
#
''':'
exec pizza.py -f "$(basename $0)" "$@"
'''
# Above assures that script is executed by pizza.py.
#
# ATTENTION: Script must be on PIZZA_SCRIPTS path defined in pizza.py's
# src/DEFAULTS.py file.
#
# Wrapper for MergeLammpsDataFiles.merge_lammps_datafiles:
#
# Takes a full, reliable system description from reffile and a newly generated
# system from datafile. Checks for sections that are missing in datafile, but
# present in reffile. These sections are either copied as they are, or specific
# mappings are applied if available:
#
# Possible applications:
#
#  TopoTools v1.7-generated LAMMPS datafiles do not contain any Coeffs sections.
# However, they contain commented index <--> name mappings for all Coeffs
# sections. Each missing section is checked for such mappings and they are
# applied if available.
#
# Next to the outfile, reffile and datafile are stripped off their comments
# and written to reffile_stripped and datafile_stripped, since pizza.py
# does not process commented files.
#
#  Ovito drops topology information, i.e. angles and dihedrals, when writing
# LAMMPS data files. Use MergeLammpsDataFiles to restore.
#
# execute with
#   pizza.py -f merge.py datafile reffile outfile
#
import sys
import argparse
import traceback
from imteksimcs.LAMMPS_data.MergeLammpsDataFiles import merge_lammps_datafiles

# recommendation from https://pizza.sandia.gov/doc/Section_basics.html#3_2

if not globals().has_key("argv"): argv = sys.argv
print("Called with '{}'".format(argv))

parser = argparse.ArgumentParser(
    description='Merges LAMMPS data file with missing entries from reference.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('datafile', nargs='?', metavar='datafile.lammps',
    default='377_SDS_on_AU_111_51x30x2_monolayer_with_counterion_100Ang_stepped.lammps',
    help="LAMMPS data file to process.")
parser.add_argument('reffile', nargs='?', metavar='reffile.lammps',
    default='377_SDS_on_AU_111_51x30x2_monolayer_with_counterion_psfgen.data',
    help="Reference data file containing complete system information.")
parser.add_argument('outfile', nargs='?', metavar='outfile.lammps',
    default='377_SDS_on_AU_111_51x30x2_monolayer_with_counterion_100Ang_stepped_parametrized.lammps',
    help="Merged output data file.")
parser.add_argument('-e', '--except', metavar='SECTION', action='append',
    default=[], dest='exceptions',
    help="Skips sections. Can be passed multiple times.")
args = parser.parse_args(argv[1:])

print("Using datafile: {}, reffile: {}; outfile: {}".format( args.datafile, args.reffile, args.outfile ) )
print("Skipping sections: {}".format( args.exceptions ) )

try:
    merge_lammps_datafiles( args.datafile, args.reffile, args.outfile, args.exceptions )
except Exception:
    traceback.print_exc()
    exit(1)

exit()
