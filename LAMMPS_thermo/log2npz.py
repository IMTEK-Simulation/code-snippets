#!/usr/bin/env python3

import re
import sys
import argparse

from numpy import genfromtxt, savez_compressed


__author__ = "Lucas Frerot"
__email__ = "lucas.frerot@imtek.uni-freiburg.de"


def get_line_nums(fp, regexp):
    """Return line numbers matching regexp."""
    return [num for num, line in enumerate(fp) if regexp.match(line)]


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Read LAMMPS log from STDIN and write numerical data to "
            "STDOUT in compressed Numpy format (binary, one frame per array)"
        ),
        epilog=(
            "Example usage:\n"
            "\n"
            "    ./log2npz.py < log.lammps > lammps_thermo_logz.npz\n"
            "\n"
            "How to use generated file:\n"
            "\n"
            "    # Loops over all the runs\n"
            "    for thermo_data in np.load('lammps_thermo_logs.npz'):\n"
            "        thermo_data['Temp'], thermo_data['PotEng'], ..."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.parse_args()

    # Load whole log in memory
    lines = sys.stdin.readlines()

    # Compute start and end of individual frames
    starts = get_line_nums(lines, re.compile("^Per MPI rank"))
    ends = get_line_nums(lines, re.compile("^Loop time"))

    # Sanity check on the number of frames, to avoid incomplete frames
    assert len(starts) == len(ends),\
        "Problem at start/end of frames"

    # Construct frame list
    frame_lines = map(lambda start, end: lines[start+1:end], starts, ends)

    # Generate composed numpy array from data (dtype=None infers each col type)
    frames = (genfromtxt(frame,
                         names=True,
                         dtype=None,
                         encoding='utf-8',
                         usemask=True) for frame in frame_lines)

    # Write to stdout in compressed Numpy format (binary)
    savez_compressed(sys.stdout.buffer, *frames)


if __name__ == "__main__":
    main()
