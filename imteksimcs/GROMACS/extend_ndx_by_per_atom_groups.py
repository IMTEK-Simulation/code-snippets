#!/usr/bin/env python3
"""Generate per-atom index groups (i.e. for pulling).

Example:
    Generate index groups for all C12 atoms in all SDS residues within system
    described in `default.top` and `default.gro`, read `in.ndx` index file and
    write extended index groups to `out.ndx`:

        $ extend_ndx_by_per_atom_groups.py --topology-file default.top --coordinates-file default.gro --residue-name SDS --atom-name C12 in.ndx out.ndx

    Same, but read from and write to stdin and stdout instead:

        $ extend_ndx_by_per_atom_groups.py --topology-file default.top --coordinates-file default.gro --residue-name SDS --atom-name C12 < in.ndx > out.ndx

    Print verbose log messages to stdout and debug log messages to log file `out.log`:

        $ extend_ndx_by_per_atom_groups.py --verbose --log out.log --topology-file default.top --coordinates-file default.gro --residue-name SDS --atom-name C12 in.ndx out.ndx
"""
import os
import logging
import numpy as np   # tested with numpy 1.15.2
import gromacs       # tested with GromacsWrapper 0.8.0
import parmed as pmd # tested with ParmEd 3.1.0

logger = logging.getLogger(__name__)

# Environment variable GMXLIB is exprected to point to GROMACS force field
# library, i.e. `/usr/share/gromacs/top`. This folder contains all 'XYZ.ff'
# force field sub-folders and the 'residuetypes.dat' residue type data base.
pmd.gromacs.GROMACS_TOPDIR = os.environ['GMXLIB']

defaults = {
    'gro' :          'default.gro',
    'top' :          'default.top',
    'ndxin' :        'in.ndx',
    'ndxout' :       'out.ndx',
    'atom_name' :    'C12',
    'residue_name' : 'SDS',
    'group_prefix' : 'pull_group_',
}

def docstring_parameter(*args, **kwargs):
    """Allows for dynamic fields in docstrings"""
    def dec(obj):
        obj.__doc__ = obj.__doc__.format(*args, **kwargs)
        return obj
    return dec

@docstring_parameter(**defaults)
def get_atom_indices(
    gro=defaults['gro'], top=defaults['top'],
    atom_name=defaults['atom_name'], residue_name=defaults['residue_name'] ):
    """
    Returns indices of all atoms `atom_name` within all residues `residue_name`.

    Parameters
    ----------
    gro : str, optional (default: '{gro:}')
        GROMACS .gro coordinates file.
    top : str, optional (default: '{top:}')
        GROMACS .top topology file.
    atom_name : str, optional (default: '{atom_name:}')
    residue_name : str, optional (default: '{residue_name:}')

    Returns
    -------
    indices: (M,) np.ndarray of int
        One-dimensional numpy array of all M identified indices.
    """
    global logger

    logger.debug("Reading topology file '{top:}'...".format(top=top))
    pmd_top_gro = pmd.gromacs.GromacsTopologyFile(top)
    logger.debug("Reading coordinates file '{gro:}'...".format(gro=gro))
    pmd_gro = pmd.gromacs.GromacsGroFile.parse(gro)
    pmd_top_gro.box = pmd_gro.box
    pmd_top_gro.positions = pmd_gro.positions

    # gromacs ndx starts at 1:
    atom_ndx = np.array([
        i+1 for i,a in enumerate(pmd_top_gro.atoms) if a.name == atom_name and a.residue.name == residue_name ], dtype=int)
    logger.debug("Identifiefd {n:} atoms named '{atom_name:}' within residue named '{residue_name:}'.".format(
        n=len(atom_ndx), atom_name=atom_name, residue_name=residue_name))
    return atom_ndx


@docstring_parameter(**defaults)
def extend_ndx_by_per_atom_groups(
    indices, ndxin=defaults['ndxin'], ndxout=defaults['ndxout'],
    group_prefix=defaults['group_prefix']):
    """
    Extends GROMACS .ndx index file by per-atom index groups.

    Parameters
    ----------
    indices: (M,) np.ndarray of int
        One-dimensional numpy array of M atom indices
        One index group is generated for each index.
    ndxin : str, optional (default: '{ndxin:}')
        GROMACS .ndx index input file
    ndxout : str, optional (default: '{ndxout:}')
        GROMACS .ndx index output file
    group_name_prefix: str, optional (default: '{group_prefix:}')
        Generated group names are compositum of `group_name_prefix`
        and zero-padded group index.

    Returns
    -------
    group_names: list of str
        Names of all generated groups,
        compositum of `group_name_prefix` and zero-padded group index.
    """
    global logger

    logger.debug("Reading index file '{ndx:}'...".format(ndx=ndxin))
    groups_ndx_in = gromacs.fileformats.NDX(ndxin)
    groups_ndx_out = gromacs.fileformats.NDX()

    group_names = []
    suffix_width = int(np.ceil(np.log10(len(indices))))
    logger.debug("{n:} groups. Groupd ID zero-padded to field width {width:}...".format(
        n=len(indices),width=suffix_width))

    for i, a in enumerate(indices):
        group_name = '{prefix:s}{suffix:0{suffix_width:d}d}'.format(
            prefix=group_prefix,suffix=i,suffix_width=suffix_width)
        groups_ndx_out[group_name] = a
        logger.debug("Creating group '{group:}' for atom index {atom:}...".format(
            group=group_name,atom=a))
        group_names.append(group_name)

    groups_ndx_in.update(groups_ndx_out)
    logger.debug("Writing index file '{ndx:}'...".format(ndx=ndxout))
    groups_ndx_in.write(ndxout)
    return group_names

def main():
    """Generate per-atom index groups (i.e. for pulling)."""
    import sys
    import argparse
    import tempfile

    # in order to have both:
    # * preformatted help text and ...
    # * automatic display of defaults
    class ArgumentDefaultsAndRawDescriptionHelpFormatter(
        argparse.ArgumentDefaultsHelpFormatter,
        argparse.RawDescriptionHelpFormatter):
        pass

    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class = ArgumentDefaultsAndRawDescriptionHelpFormatter)

    ### positional arguments
    parser.add_argument('infile', metavar='IN', default=None, nargs='?',
                        help='GROMACS .ndx index file. If not specified, read from stdin.')
    parser.add_argument('outfile', metavar='OUT', default=None, nargs='?',
                        help='GROMACS .ndx index file. If not specified, read from stdout.')

    ### options
    parser.add_argument('--topology-file','-t',
                        default=defaults['top'], type=str,
                        metavar='TOP', required=False, dest="topology_file",
                        help='GROMACS .top topology file')

    parser.add_argument('--coordinates-file','-c',
                        default=defaults['gro'], type=str,
                        metavar='GRO', required=False, dest="coordinates_file",
                        help='GROMACS .gro coordinates file')

    parser.add_argument('--atom-name','-a',
                        default=defaults['atom_name'], type=str,
                        metavar='ATM', required=False, dest="atom_name",
                        help='name of atom within residue')

    parser.add_argument('--residue-name','-r',
                        default=defaults['residue_name'], type=str,
                        metavar='RES', required=False, dest="residue_name",
                        help='name of residue containing atom of interest')

    parser.add_argument('--group-prefix','-g',
                        default=defaults['group_prefix'], type=str,
                        metavar='GROUP_PREFIX', required=False, dest="group_prefix",
                        help='name of atom within residue')

    ### logging
    parser.add_argument('--debug', default=False, required=False,
                        action='store_true', dest="debug", help='debug flag')
    parser.add_argument('--verbose', default=False, required=False,
                        action='store_true', dest="verbose", help='verbose flag')
    parser.add_argument('--log', required=False, nargs='?', dest="log",
                        default=None, const='std.log', metavar='LOG',
                        help='Write log file std.log, optionally specify log file name')

    try:
        import argcomplete
        argcomplete.autocomplete(parser)
        # This supports bash autocompletion. To enable this, pip install
        # argcomplete, activate global completion
    except ImportError:
        pass

    args = parser.parse_args()

    if args.debug:
        loglevel = logging.DEBUG
    elif args.verbose:
        loglevel = logging.INFO
    else:
        loglevel = logging.WARNING

    logformat  = "[ %(asctime)s - %(funcName)s() - %(filename)s:%(lineno)s ]: %(levelname)s - %(message)s"

    logging.basicConfig(level=loglevel,
                        format=logformat)

    # explicitly modify the root logger (necessary?)
    logger = logging.getLogger()
    logger.setLevel(loglevel)

    # remove all handlers
    for h in logger.handlers: logger.removeHandler(h)

    # create and append custom handles
    ch = logging.StreamHandler()
    formatter = logging.Formatter(logformat)
    ch.setFormatter(formatter)
    ch.setLevel(loglevel)
    logger.addHandler(ch)

    if args.log:
        fh = logging.FileHandler(args.log)
        fh.setFormatter(formatter)
        fh.setLevel(loglevel)
        logger.addHandler(fh)

    logger.info("This is '{}' : '{}'.".format(__file__,__name__))
    for arg, val in sorted(vars(args).items()):
        logging.info("Argument {arg:}: '{val:}' ({type:})".format(arg=arg,val=val,type=type(val)))

    indices = get_atom_indices(
        gro=args.coordinates_file, top=args.topology_file,
        atom_name=args.atom_name, residue_name=args.residue_name )
    logger.debug("Indices: [{:s}]".format(','.join([ str(i) for i in indices])))

    # Handing file objects instead of file names to ParmEd does not work well.
    # Thus, we allow for reading from and writing to stdin and stdout by
    # rerouting through named temporary files:
    if not args.infile:
        logger.info('Reading from stdin.')
        tmpin = tempfile.NamedTemporaryFile(mode='w+',suffix='.ndx')
        tmpin.write(sys.stdin.read())
        infile = tmpin.name
    else:
        infile = args.infile
        logger.info("Reading from '{:s}'.".format(infile))

    if not args.outfile:
        tmpout = tempfile.NamedTemporaryFile(mode='w+',suffix='.ndx')
        outfile = tmpout.name
        logger.info('Writing to stdout.')
    else:
        outfile = args.outfile
        logger.info("Writing to '{:s}'.".format(outfile))

    group_names = extend_ndx_by_per_atom_groups(
        indices, ndxin=infile, ndxout=outfile, group_prefix=args.group_prefix )

    if not args.outfile: # no explicit output file, write results to stdout
        sys.stdout.write(tmpout.read())

    logger.debug("Index groups: [{:s}]".format(','.join(group_names)))

    logger.info("Done.")
    # Temporary files are released automatically


if __name__ == '__main__':
    main()
