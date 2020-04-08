#!/usr/bin/env python
"""
Strips all comments from LAMMPS data file to be processed by pizza.py
"""
from imteksimcs.LAMMPS_data.strip_comments import strip_comments

def main():
    import argparse

    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('infile',
        help='LAMMPS input data file')
    parser.add_argument('outfile',
        help='LAMMPS output data file readable by pizza.py')
    args = parser.parse_args()

    strip_comments(args.infile, args.outfile)

if __name__ == '__main__':
    main()

