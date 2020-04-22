#!/usr/bin/env python3
"""Generate an index group excluding a specific residue."""
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
    'gro':          'default.gro',
    'top':          'default.top',
    'ndxin':        'in.ndx',
    'ndxout':       'out.ndx',
    'residue_name': 'AUM',
    'group_name':   'non-Substrate',
}

def docstring_parameter(*args, **kwargs):
    """Allows for dynamic fields in docstrings"""
    def dec(obj):
        obj.__doc__ = obj.__doc__.format(*args, **kwargs)
        return obj
    return dec

@docstring_parameter(**defaults)
def get_non_residue_indices(
    gro=defaults['gro'], top=defaults['top'],
    residue_name=defaults['residue_name']):
    """
    Returns indices of all atoms not within residue `residue_name`.

    Parameters
    ----------
    gro : str, optional (default: '{gro:}')
        GROMACS .gro coordinates file.
    top : str, optional (default: '{top:}')
        GROMACS .top topology file.
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

    atom_ndx = np.array([
        i+1 for i,a in enumerate(pmd_top_gro.atoms) if a.residue.name != residue_name])
    # gromacs ndx starts at 1:
    logger.debug("Identifiefd {n:} atoms not in residue named '{residue_name:}'.".format(
        n=len(atom_ndx), residue_name=residue_name))
    return atom_ndx


@docstring_parameter(**defaults)
def extend_ndx_by_group(
    indices, ndxin=defaults['ndxin'], ndxout=defaults['ndxout'],
    group_name=defaults['group_name']):
    """
    Extends GROMACS .ndx index file by some index group.

    Parameters
    ----------
    indices: (M,) np.ndarray of int
        One-dimensional numpy array of M atom indices
        One index group is generated for each index.
    ndxin : str, optional (default: '{ndxin:}')
        GROMACS .ndx index input file
    ndxout : str, optional (default: '{ndxout:}')
        GROMACS .ndx index output file
    group_name: str, optional (default: '{group_name:}')
        Name for new group.

    Returns
    -------
    group_name: str
        Name of new group.
    """
    global logger

    logger.debug("Reading index file '{ndx:}'...".format(ndx=ndxin))
    groups_ndx_in = gromacs.fileformats.NDX(ndxin)
    groups_ndx_out = gromacs.fileformats.NDX()

    groups_ndx_out[group_name] = indices
    logger.debug("Creating group '{group:}' with {n:d} indices...".format(
            group=group_name, n=len(indices)))

    groups_ndx_in.update(groups_ndx_out)
    logger.debug("Writing index file '{ndx:}'...".format(ndx=ndxout))
    groups_ndx_in.write(ndxout)
    return group_name
