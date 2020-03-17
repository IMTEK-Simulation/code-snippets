#!/usr/bin/env ovitos
#
# ncfilter.py
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


import logging

# adapted from
# https://ovito.org/manual/python/modules/ovito_modifiers.html#ovito.modifiers.LoadTrajectoryModifier
def main():
    import argparse
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

def parseRanges(s):
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

    # if '-' in s:
    #     i = s.index('-')
    #     return range(*map(int, s.split('-')))
    # return map(int, s.split(','))


def extractFrames(
    topology_file   = 'datafile.lammps', # lammps data file
    trajectory_file = 'trajectory.nc', # netcdf trajectory
    frames          = None, # list of frame numbers
    outfile_name    = None, # list of outfiles, same length as frames
    outfile_pattern = 'frame_{:08d}.lammps',
    atom_style      = 'full'):
    """Extracts selected or all frames from a NetCDF trajectory.

    Parameters
    ----------
    topology_file : str, optional
        LAMMPS data file, per default atom style 'full'
    trajectory_file: str, optional
        Trajectory file (atom positions) in NetCDF format
    frames: :obj:`list` of :obj:`int`, optional
        Arbitrary selection of frames. Per default, all frames extracted.
    outfile_name :obj:`dict` of :obj:`int` : :obj:`str`, optional
        Dictionary assigning output file name for every selected frame.
    outfile_pattern: str, optional
        If `outfile_name` not given, output file names are constructed from
        this pattern and frame number. Must conatin some {:d} integer
        placeholder, i.e. 'frame_{:08d}.lammps'
    atom_style: str, optional
        LAMMPS atom sytle as understood by ovitos
    """

    from ovito.io import import_file, export_file
    from ovito.modifiers import LoadTrajectoryModifier

    logger = logging.getLogger('netcdf2data.extractFrames')
    logger.info("Reading topology from '{:s}'".format(topology_file))
    # Load static topology data from a LAMMPS data file.
    node = import_file(topology_file, atom_style = atom_style)

    logger.info("Reading trajectory from '{:s}'".format(trajectory_file))
    # Load atom trajectories from separate LAMMPS dump file.
    traj_mod = LoadTrajectoryModifier()
    traj_mod.source.load(trajectory_file, multiple_frames = True)

    # Insert modifier into modification pipeline.
    node.modifiers.append(traj_mod)

    if frames is None: frames = range(traj_mod.source.num_frames)
    if outfile_name is None:
        outfile_name = {
            frame: outfile_pattern.format(frame) for frame in frames }

    for frame in frames:
        try:
            logger.info("Storing frame {:08d} to '{:s}'".format(
                frame, outfile_name[frame]))
        except IndexError:
            logger.error(
                "No outfile name specified for frame {:d}!".format(frame))
            pass # out of range

        export_file(
            node, outfile_name[frame],
            'lammps_data', atom_style = atom_style, frame = frame)

    return

if __name__ == '__main__':
    main()
