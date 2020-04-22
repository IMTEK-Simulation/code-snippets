#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK
#
# gmx_tools.py
#
# Copyright (C) 2020 IMTEK Simulation
# Author: Johannes Hoermann, johannes.hoermann@imtek.uni-freiburg.de
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""Manipulates GROMACS-related file formats.

Examples
--------

    gmx_tools --verbose make pull_groups \\
        --topology-file default.top --coordinates-file default.gro \\
        --atom-name C12 --residue-name SDS --reference-group-name Substrate -- \\
        in.ndx out.ndx pull.mdp.template pull.mdp
"""

import argparse
import os
import sys  # for stdout and stderr
import tempfile
import logging

from imteksimcs.GROMACS import extend_ndx_by_inversion_group
from imteksimcs.GROMACS import extend_ndx_by_per_atom_groups
from imteksimcs.GROMACS import add_pull_groups_to_mdp


# in order to have both:
# * preformatted help text and ...
# * automatic display of defaults
class ArgumentDefaultsAndRawDescriptionHelpFormatter(
        argparse.ArgumentDefaultsHelpFormatter,
        argparse.RawDescriptionHelpFormatter):
    """Allows for both preformatted help and automatic defaults display."""
    pass


# action blocks

def extend_ndx_by_inversion_group_action(args):
    """Extend existing .ndx index file by negated residue selection.

    Returns
    -------
        str: group name
    """
    logger = logging.getLogger(__name__)
    indices = extend_ndx_by_inversion_group.get_non_residue_indices(
        gro=args.coordinates_file, top=args.topology_file,
        residue_name=args.residue_name )
    logger.debug("Indices: [{:s}]".format(','.join([str(i) for i in indices])))

    # Handing file objects instead of file names to ParmEd does not work well.
    # Thus, we allow for reading from and writing to stdin and stdout by
    # rerouting through named temporary files:
    if not args.ndxin:
        logger.info('Reading from stdin.')
        tmpin = tempfile.NamedTemporaryFile(mode='w+', suffix='.ndx')
        tmpin.write(sys.stdin.read())
        infile = tmpin.name
    else:
        infile = args.ndxin
        logger.info("Reading from '{:s}'.".format(infile))

    if not args.ndxout:
        tmpout = tempfile.NamedTemporaryFile(mode='w+', suffix='.ndx')
        outfile = tmpout.name
        logger.info('Writing to stdout.')
    else:
        outfile = args.ndxout
        logger.info("Writing to '{:s}'.".format(outfile))

    group_name = extend_ndx_by_inversion_group.extend_ndx_by_group(
        indices, ndxin=infile, ndxout=outfile, group_name=args.group_name)

    if not args.ndxout: # no explicit output file, write results to stdout
        sys.stdout.write(tmpout.read())

    logger.debug("Index group: [{:s}]".format(','.join(group_name)))
    logger.info("Done.")

    return group_name
    # Temporary files are released automatically

def extend_ndx_by_per_atom_groups_action(args):
    """Extend existing .ndx index file by per-atom groups.

    Returns
    -------
        list of str: group names
    """
    logger = logging.getLogger(__name__)
    indices = extend_ndx_by_per_atom_groups.get_atom_indices(
        gro=args.coordinates_file, top=args.topology_file,
        atom_name=args.atom_name, residue_name=args.residue_name )
    logger.debug("Indices: [{:s}]".format(','.join([ str(i) for i in indices])))

    # Handing file objects instead of file names to ParmEd does not work well.
    # Thus, we allow for reading from and writing to stdin and stdout by
    # rerouting through named temporary files:
    if not args.ndxin:
        logger.info('Reading from stdin.')
        tmpin = tempfile.NamedTemporaryFile(mode='w+',suffix='.ndx')
        tmpin.write(sys.stdin.read())
        infile = tmpin.name
    else:
        infile = args.ndxin
        logger.info("Reading from '{:s}'.".format(infile))

    if not args.ndxout:
        tmpout = tempfile.NamedTemporaryFile(mode='w+',suffix='.ndx')
        outfile = tmpout.name
        logger.info('Writing to stdout.')
    else:
        outfile = args.ndxout
        logger.info("Writing to '{:s}'.".format(outfile))

    group_names = extend_ndx_by_per_atom_groups.extend_ndx_by_per_atom_groups(
        indices, ndxin=infile, ndxout=outfile, group_prefix=args.group_prefix )

    if not args.ndxout: # no explicit output file, write results to stdout
        sys.stdout.write(tmpout.read())

    logger.debug("Index groups: [{:s}]".format(','.join(group_names)))
    logger.info("Done.")

    return group_names
    # Temporary files are released automatically


def add_pull_groups_to_mdp_action(args):
    """Adds pull groups to existing .mdp parameter file."""
    logger = logging.getLogger(__name__)

    if not args.mdpin:
        logger.info('Reading from stdin.')
        tmpin = tempfile.NamedTemporaryFile(mode='w+',suffix='.mdp')
        tmpin.write(sys.stdin.read())
        infile = tmpin.name
    else:
        infile = args.mdpin
        logger.info("Reading from '{:s}'.".format(infile))

    if not args.mdpout:
        tmpout = tempfile.NamedTemporaryFile(mode='w+',suffix='.mdp')
        outfile = tmpout.name
        logger.info('Writing to stdout.')
    else:
        outfile = args.mdpout
        logger.info("Writing to '{:s}'.".format(outfile))

    add_pull_groups_to_mdp.add_pull_groups_to_mdp(
        pull_group_names=args.pull_group_names,
        reference_group_name=args.reference_group_name,
        nsteps=args.nsteps,
        k=args.k,
        rate=args.rate,
        mdpin=infile,
        mdpout=outfile)

    if not args.mdpout: # no explicit output file, write results to stdout
        sys.stdout.write(tmpout.read())
    logger.info("Done.")


def make_per_atom_pull_groups_action(args):
    """Add pull groups to existing .ndx and .mdp files."""
    logger = logging.getLogger(__name__)
    args.pull_group_names = extend_ndx_by_per_atom_groups_action(args)
    logger.info("Identified {:d} groups".format(len(args.pull_group_names)))
    add_pull_groups_to_mdp_action(args)

# positional arguments blocks

def positional_arguments_ndx(
        parser, infile='ndxin', outfile='ndxout',
        inmeta='NDXIN', outmeta='NDXOUT', optional=True):
    """Add .ndx file positional arguments to parser"""
    nargs = '?' if optional else None
    parser.add_argument(infile, metavar=inmeta, default=None, nargs=nargs,
                        help='GROMACS .ndx index file. If not specified, read from stdin.')
    parser.add_argument(outfile, metavar=outmeta, default=None, nargs=nargs,
                        help='GROMACS .ndx index file. If not specified, read from stdout.')



def positional_arguments_mdp(
        parser, infile='mdpin', outfile='mdpout',
        inmeta='MDPIN', outmeta='MDPOUT', optional=True):
    nargs = '?' if optional else None
    """Add .mdp file positional arguments to parser"""
    parser.add_argument(infile, metavar=inmeta, default=None, nargs=nargs,
                        help='GROMACS .mdp file. If not specified, read from stdin.')
    parser.add_argument(outfile, metavar=outmeta, default=None, nargs=nargs,
                        help='GROMACS .mdp index file. If not specified, read from stdout.')


# optional arguments blocks
def optional_arguments_extend_ndx_by_inversion_group(parser):
    """Add extend_ndx_by_inversion_group optional arguments to parser"""
    ### options
    parser.add_argument('--topology-file','-t',
                        default=extend_ndx_by_inversion_group.defaults['top'], type=str,
                        metavar='TOP', required=False, dest="topology_file",
                        help='GROMACS .top topology file')

    parser.add_argument('--coordinates-file','-c',
                        default=extend_ndx_by_inversion_group.defaults['gro'], type=str,
                        metavar='GRO', required=False, dest="coordinates_file",
                        help='GROMACS .gro coordinates file')

    parser.add_argument('--residue-name','-r',
                        default=extend_ndx_by_inversion_group.defaults['residue_name'], type=str,
                        metavar='RES', required=False, dest="residue_name",
                        help='name of residue to exclude / negate')

    parser.add_argument('--group-name','-g',
                        default=extend_ndx_by_inversion_group.defaults['group_name'], type=str,
                        metavar='GROUP_NAME', required=False, dest="group_name",
                        help='index group name')


def optional_arguments_extend_ndx_by_per_atom_groups(parser):
    """Add extend_ndx_by_per_atom_groups optional arguments to parser"""
    ### options
    parser.add_argument('--topology-file','-t',
                        default=extend_ndx_by_per_atom_groups.defaults['top'], type=str,
                        metavar='TOP', required=False, dest="topology_file",
                        help='GROMACS .top topology file')

    parser.add_argument('--coordinates-file','-c',
                        default=extend_ndx_by_per_atom_groups.defaults['gro'], type=str,
                        metavar='GRO', required=False, dest="coordinates_file",
                        help='GROMACS .gro coordinates file')

    parser.add_argument('--atom-name','-a',
                        default=extend_ndx_by_per_atom_groups.defaults['atom_name'], type=str,
                        metavar='ATM', required=False, dest="atom_name",
                        help='name of atom within residue')

    parser.add_argument('--residue-name','-r',
                        default=extend_ndx_by_per_atom_groups.defaults['residue_name'], type=str,
                        metavar='RES', required=False, dest="residue_name",
                        help='name of residue containing atom of interest')

    parser.add_argument('--group-prefix','-g',
                        default=extend_ndx_by_per_atom_groups.defaults['group_prefix'], type=str,
                        metavar='GROUP_PREFIX', required=False, dest="group_prefix",
                        help='name of atom within residue')


def optional_arguments_add_pull_groups_to_mdp(parser):
    """Add add_pull_groups_to_mdp optional arguments to parser"""
    parser.add_argument('--reference-group-name',
                        default=add_pull_groups_to_mdp.defaults['reference_group_name'],
                        type=str,
                        metavar='GRP', required=False,
                        dest="reference_group_name",
                        help='Reference group name.')
    parser.add_argument('--pull-group-names',
                        default=[],
                        type=str,
                        nargs='+',
                        metavar='GRP', required=False,
                        dest="pull_group_names",
                        help='Pull group names.')
    parser.add_argument('--nsteps','-n',
                        default=add_pull_groups_to_mdp.defaults['k'],
                        type=float,
                        metavar='N', required=False, dest="nsteps",
                        help='Overrides number of steps in .mdp file.')
    parser.add_argument('-k',
                        default=add_pull_groups_to_mdp.defaults['k'],
                        type=float,
                        metavar='K', required=False, dest="k",
                        help='Harmonic pulling spring constant [kJ mol^-1 nm^-2].')
    parser.add_argument('--rate',
                        default=add_pull_groups_to_mdp.defaults['rate'],
                        type=float,
                        metavar='R', required=False, dest="rate",
                        help='Harmonic pulling rate [nm ps^-1].')


# subparser blocks
def subparser_extend_ndx_by_inversion_group(subparsers,
        command='extend_by_inversion_group', aliases=['invert']):
    parser = subparsers.add_parser(
        command, aliases=aliases,
        help='Extend exisiting .ndx index files by negating a residue selection.',
        formatter_class=ArgumentDefaultsAndRawDescriptionHelpFormatter)
    positional_arguments_ndx(parser)
    optional_arguments_extend_ndx_by_inversion_group(parser)
    parser.set_defaults(func=extend_ndx_by_inversion_group_action)


def subparser_extend_ndx_by_per_atom_groups(subparsers,
        command='extend_by_per_atom_groups', aliases=['atom_groups', 'groups']):
    parser = subparsers.add_parser(
        command, aliases=aliases,
        help='Extend exisiting .ndx index files by per-atom groups.',
        formatter_class=ArgumentDefaultsAndRawDescriptionHelpFormatter)
    positional_arguments_ndx(parser)
    optional_arguments_extend_ndx_by_per_atom_groups(parser)
    parser.set_defaults(func=extend_ndx_by_per_atom_groups_action)


def subparser_add_pull_groups_to_mdp(subparsers,
        command="add_pull_groups", aliases=['pull_groups']):
    """Set up add_pull_groups_to_mdp parser"""
    parser = subparsers.add_parser(
        command, aliases=aliases,
        help='Extend .mdp parameter file by pulling groups.',
        formatter_class=ArgumentDefaultsAndRawDescriptionHelpFormatter)
    positional_arguments_mdp(parser)
    optional_arguments_add_pull_groups_to_mdp(parser)
    parser.set_defaults(func=add_pull_groups_to_mdp_action)


def subparser_make_per_atom_pull_groups(subparsers,
        command="make_per_atom_pull_groups", aliases=['pull_groups']):
    """Set up make_per_atom_pull_groups parser"""
    parser = subparsers.add_parser(
        command, aliases=aliases,
        help='Add per-atom pull groups to .ndx and .mdp.',
        formatter_class=ArgumentDefaultsAndRawDescriptionHelpFormatter)
    positional_arguments_ndx(parser, optional=False)
    positional_arguments_mdp(parser, optional=False)
    optional_arguments_extend_ndx_by_per_atom_groups(parser)
    optional_arguments_add_pull_groups_to_mdp(parser)
    parser.set_defaults(func=make_per_atom_pull_groups_action)


def main():
    """Manipulate GROMACS-related file formats from command line."""

    # root parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=ArgumentDefaultsAndRawDescriptionHelpFormatter)

    # global logging options
    parser.add_argument('--debug', default=False, required=False,
                        action='store_true', dest="debug", help='debug')
    parser.add_argument('--verbose', default=False, required=False,
                        action='store_true', dest="verbose", help='verbose')

    parser.add_argument('--log', required=False, nargs='?', dest="log",
                        default=None, const='out.log', metavar='LOG',
                        help='Write out.log, optionally specify log name')

    # 1st level subparsers
    subparsers = parser.add_subparsers(
        help='command', dest='command', title='commands')
    ndx_parser = subparsers.add_parser(
        'ndx',
        help='Manipulate .ndx index file.',
        formatter_class=ArgumentDefaultsAndRawDescriptionHelpFormatter)
    mdp_parser = subparsers.add_parser(
        'mdp',
        help='Manipulate .mdp parameter file.',
        formatter_class=ArgumentDefaultsAndRawDescriptionHelpFormatter)
    make_parser = subparsers.add_parser(
        'make',
        help='Carry out a whole set of manipulations on several files.',
        formatter_class=ArgumentDefaultsAndRawDescriptionHelpFormatter)


    # 2nd level subparsers
    ndx_subparsers = ndx_parser.add_subparsers(
        help='sub-command', dest='subcommand', title='ndx subcommands')
    mdp_subparsers = mdp_parser.add_subparsers(
        help='sub-command', dest='subcommand', title='mdp subcommands')
    make_subparsers = make_parser.add_subparsers(
        help='sub-command', dest='subcommand', title='make subcommands')

    subparser_extend_ndx_by_inversion_group(ndx_subparsers)
    subparser_extend_ndx_by_per_atom_groups(ndx_subparsers)

    subparser_add_pull_groups_to_mdp(mdp_subparsers)
    subparser_make_per_atom_pull_groups(make_subparsers)

    try:
        import argcomplete
        argcomplete.autocomplete(parser)
        # This supports bash autocompletion. To enable this, pip install
        # argcomplete, activate global completion
    except ImportError:
        pass

    args = parser.parse_args()

    # simple log format for non-debug output
    logformat = "%(levelname)s: %(message)s"
    if args.debug:
        # detailed log format for debug output
        logformat = (
            "[%(asctime)s-%(funcName)s()-%(filename)s:%(lineno)s]"
            " %(levelname)s: %(message)s"
        )
        loglevel = logging.DEBUG
    elif args.verbose:
        loglevel = logging.INFO
    else:
        loglevel = logging.WARNING

    logging.basicConfig(level=loglevel, format=logformat)

    # explicitly modify the root logger (necessary?)
    logger = logging.getLogger()
    logger.setLevel(loglevel)

    # remove all handlers
    for h in logger.handlers:
        logger.removeHandler(h)

    # create and append custom handles

    # only info & debug to stdout
    def stdout_filter(record):
        return record.levelno <= logging.INFO

    stdouth = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter(logformat)
    stdouth.setFormatter(formatter)
    stdouth.setLevel(loglevel)
    stdouth.addFilter(stdout_filter)

    # always log warnings and errors to stderr
    stderrh = logging.StreamHandler(sys.stderr)
    stderrh.setFormatter(formatter)
    stderrh.setLevel(logging.WARNING)

    logger.addHandler(stdouth)
    logger.addHandler(stderrh)

    # log everything to file if desired
    if args.log:
        fh = logging.FileHandler(args.log)
        fh.setFormatter(formatter)
        fh.setLevel(loglevel)
        logger.addHandler(fh)

    logger.info("This is '{}' : '{}'.".format(__file__,__name__))
    for arg, val in sorted(vars(args).items()):
        logging.info("Argument {arg:}: '{val:}' ({type:})".format(arg=arg,val=val,type=type(val)))

    if args.command is None:
        # if no command supplied, print help
        parser.print_help()
    elif 'func' not in args:
        parser.print_help()
    else:
        args.func(args)


if __name__ == '__main__':
    main()
