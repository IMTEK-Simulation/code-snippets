#!/usr/bin/env python
#
# build_smamp_pdb.py
#
# Copyright (C) 2020 IMTEK Simulation
# Authors: Sarah Elsayed, Johannes Hoermann, Richard Leute, Jan Mees
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
"""Build SMAMP molecule with terA -- (N-2) OXO0 -- terB.

Examples
--------

    build_smamp_pdb.py 10

builds 'SMAMP10.pdb' from 'terB.pdb', 8 'OXO.pdb' and 'terA.pdb'.

    build_smamp_pdb.py --up-offset 0.1 0.1 0.1 --down-offset -0.1 -0.1 -0.1 \
        --repeat-unit another_OXO.pdb -- 5 another_SMAMP.pdb

builds 'another_SMAMP.pdb' from 'terB.pdb', 3 'another_OXO.pdb' and 'terA.pdb'
with custom offsets at the residue junctions."""
import argparse
import mbuild as mb

import warnings
warnings.filterwarnings('ignore')


DEFAULTS = {
    'up_anchor':   'CD6',
    'down_anchor': 'CD1',
    'up_offset':   [0, 0.154/2, 0],
    'down_offset': [0, -0.154/2, 0],
    'terB_infile': 'terB.pdb',
    'terA_infile': 'terA.pdb',
    'repeat_unit_infile': 'OXO0.pdb',
    'outfile': 'SMAMP{:02d}.pdb',
}


# in analogy to http://mosdef.org/mbuild/tutorials/tutorial_polymers.html
class OXO0(mb.Compound):
    def __init__(self,
            infile=DEFAULTS['repeat_unit_infile'],
            up_offset=DEFAULTS['up_offset'],
            down_offset=DEFAULTS['down_offset'],
            up_anchor=DEFAULTS['up_anchor'],
            down_anchor=DEFAULTS['down_anchor']):
        super(OXO0, self).__init__()

        mb.load(infile, compound=self)

        atom_names = [a.name for a in self]
        nterminal = atom_names.index(down_anchor)
        cterminal = atom_names.index(up_anchor)

        self.add(mb.Port(anchor=self[nterminal]), label='down')
        mb.translate(self['down'], up_offset)

        self.add(mb.Port(anchor=self[cterminal]), label='up')
        mb.translate(self['up'], down_offset)


class terB(mb.Compound):
    def __init__(self,
            infile=DEFAULTS['terB_infile'],
            down_offset=DEFAULTS['down_offset'],
            down_anchor=DEFAULTS['down_anchor']):
        super(terB, self).__init__()
        mb.load(infile, compound=self)

        atom_names = [a.name for a in self]
        nterminal = atom_names.index(down_anchor)

        self.add(mb.Port(anchor=self[nterminal]), label='down')
        mb.translate(self['down'], down_offset)


class terA(mb.Compound):
    def __init__(self,
            infile=DEFAULTS['terA_infile'],
            up_offset=DEFAULTS['up_offset'],
            up_anchor=DEFAULTS['up_anchor']):
        super(terA, self).__init__()

        mb.load(infile, compound=self)

        atom_names = [a.name for a in self]
        cterminal = atom_names.index(up_anchor)

        self.add(mb.Port(anchor=self[cterminal]), label='up')
        mb.translate(self['up'], up_offset)


class SMAMP(mb.Compound):
    def __init__(self, n,
            terB_infile=DEFAULTS['terB_infile'],
            terA_infile=DEFAULTS['terA_infile'],
            repeat_unit_infile=DEFAULTS['repeat_unit_infile'],
            up_offset=DEFAULTS['up_offset'],
            down_offset=DEFAULTS['down_offset'],
            up_anchor=DEFAULTS['up_anchor'],
            down_anchor=DEFAULTS['down_anchor']):

        super(SMAMP, self).__init__()

        last_monomer = terB(
            infile=terB_infile,
            down_offset=down_offset,
            down_anchor=down_anchor)
        self.add(last_monomer)

        for i in range(n-2):  # terminals each include one repeat unit
            current_monomer = OXO0(
                infile=repeat_unit_infile,
                down_offset=down_offset,
                down_anchor=down_anchor,
                up_offset=up_offset,
                up_anchor=up_anchor)

            mb.force_overlap(move_this=current_monomer,
                             from_positions=current_monomer['up'],
                             to_positions=last_monomer['down'])
            self.add(current_monomer)
            last_monomer = current_monomer

        current_monomer = terA(
            infile=terA_infile,
            up_offset=up_offset,
            up_anchor=up_anchor)

        mb.force_overlap(move_this=current_monomer,
                         from_positions=current_monomer['up'],
                         to_positions=last_monomer['down'])
        self.add(current_monomer)


def main():
    """Build SMAMP molecule with terB -- (N-2) OXO0 -- terA.

    Examples:

        build_smamp_pdb.py 10

    builds 'SMAMP10.pdb' from 'terB.pdb', 8 'OXO.pdb' and 'terA.pdb'.

        build_smamp_pdb.py --up-offset 0.1 0.1 0.1 --down-offset -0.1 -0.1 -0.1 \
            --repeat-unit another_OXO.pdb -- 5 another_SMAMP.pdb

    builds 'another_SMAMP.pdb' from 'terB.pdb', 3 'another_OXO.pdb' and 'terA.pdb'
    with custom offsets at the residue junctions."""

    # in order to have both:
    # * preformatted help text and ...
    # * automatic display of defaults
    class ArgumentDefaultsAndRawDescriptionHelpFormatter(
            argparse.ArgumentDefaultsHelpFormatter,
            argparse.RawDescriptionHelpFormatter):
        pass

    class StoreAsNumpyArray(argparse._StoreAction):
        def __call__(self, parser, namespace, values, option_string=None):
            values = np.array(values, ndmin=1)
            return super().__call__(
                parser, namespace, values, option_string)

    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=ArgumentDefaultsAndRawDescriptionHelpFormatter)

    parser.add_argument('chain_length', metavar='N',
                        default=3, type=int,
                        help='desired chain length.')

    parser.add_argument('outfile', metavar='OUT', type=str, nargs='?',
                        default=DEFAULTS['outfile'],
                        help='.pdb format output file.')

    parser.add_argument('--terA', metavar='PDB', type=str,
                        default=DEFAULTS['terA_infile'], dest='terA_infile',
                        help='terA residue .pdb input file.')

    parser.add_argument('--terB', metavar='PDB', type=str,
                        default=DEFAULTS['terB_infile'], dest='terB_infile',
                        help='terB residue .pdb input file.')

    parser.add_argument('--repeat-unit', metavar='PDB',
                        default=DEFAULTS['repeat_unit_infile'], dest='repeat_unit_infile',
                        help='repeat unit residue .pdb input file.')

    parser.add_argument('--up-offset', default=DEFAULTS['up_offset'], 
                        type=float, nargs=3,
                        metavar=('X','Y','Z'), required=False, dest="up_offset",
                        help='Offset at "up" anchor (c-terminal).')

    parser.add_argument('--down-offset', default=DEFAULTS['down_offset'],
                        type=float, nargs=3,
                        metavar=('X','Y','Z'), required=False, dest="down_offset",
                        help='Offset at "down" anchor (n-terminal).')

    parser.add_argument('--up-anchor', default=DEFAULTS['up_anchor'], type=str,
                        metavar='ATM', required=False, dest="up_anchor",
                        help='"up" anchor (c-terminal).')

    parser.add_argument('--down-anchor', default=DEFAULTS['down_anchor'], type=str,
                        metavar='ATM', required=False, dest="down_anchor",
                        help='"down" anchor (n-terminal).')

    args = parser.parse_args()

    polymer = SMAMP(args.chain_length,
            terA_infile=args.terA_infile,
            terB_infile=args.terB_infile,
            repeat_unit_infile=args.repeat_unit_infile,
            up_offset=args.up_offset,
            down_offset=args.down_offset,
            up_anchor=args.up_anchor,
            down_anchor=args.down_anchor)

    # for assigning residue names according to class names as above
    pmd_smamp = polymer.to_parmed(infer_residues=True)
    # for writing 4-letter instead of 3-letter residue names, charmm=True
    pmd_smamp.write_pdb(
        args.outfile.format(args.chain_length), standard_resnames=False, charmm=True)


if __name__ == '__main__':
    main()
