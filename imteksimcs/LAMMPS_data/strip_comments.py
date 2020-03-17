#!/usr/bin/env python
"""
Strips all comments from LAMMPS data file to be processed by pizza.py
"""

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

def strip_comments(infile, outfile):
    """Removes all trailing comments from a LAMMPS data file.
       Necessary to make them pizza.py-processible"""
    import re
    regex = re.compile(r"\s*#.*$")
    with open(infile) as i:
        with open(outfile, 'w') as o:
            for line in i:
                line = regex.sub('',line)
                o.write(line)

if __name__ == '__main__':
    main()

