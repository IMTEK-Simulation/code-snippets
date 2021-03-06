#!/usr/bin/env python
"""
Converts LAMMPS pair_coeff lines to 'pair_style hybrid' format.
"""
# 29 Mar 2019, converts LAMMPS force field pair_coeff entries to hybrid style

from imteksimcs.LAMMPS_data.to_hybrid import convert


def main():
    import argparse

    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--pair-style',
        help="LAMMPS pair style to attach to pair_coeff entries.",
        default='lj/charmmfsw/coul/long')
    parser.add_argument('--regular-expression',
        help="""Regular expression to match pair_coeff lines.
          Insertion of 'pair-style' between first and second group.""",
        default='(^pair_coeff\s+\d+\s+\d+)\s+([^\s].*$)')
    parser.add_argument('infile',
        help='LAMMPS input file holding force field parameter',
        default='coeff.input', nargs='?')
    parser.add_argument('outfile',
        help='LAMMPS input file to store modified force field parameter entries',
        default='coeff_hybrid.input', nargs='?')
    args = parser.parse_args()

    convert(args.infile, args.outfile, args.pair_style, args.regular_expression)


if __name__ == '__main__':
    main()
