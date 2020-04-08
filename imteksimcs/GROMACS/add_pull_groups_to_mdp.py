#!/usr/bin/env python3
"""Add pull groups to mdp."""
import os
import logging
import gromacs       # tested with GromacsWrapper 0.8.0

defaults = {
    'reference_group_name': 'Surfactant',
    'nsteps' : None,
    'mdpin' : 'in.mdp',
    'mdpout' : 'out.mdp',
    'k' :      1000,
    'rate':    -0.1,
    '
}

def docstring_parameter(*args, **kwargs):
    """Allows for dynamic fields in docstrings"""
    def dec(obj):
        obj.__doc__ = obj.__doc__.format(*args, **kwargs)
        return obj
    return dec


@docstring_parameter(**defaults)
def add_pull_groups_to_mdp(
    pull_group_names,
    reference_group_name=defaults['reference_group_name'],
    nsteps=defaults['nsteps'],
    k=defaults['k'],
    rate=defaults['rate'],
    mdpin=defaults['mdpin'],
    mdpout=defaults['mdpout']):
    """
    Extends GROMACS .mdp parameter file by pull groups.

    Parameters
    ----------
    pull_group_names: (M,) list of str
        One-dimensional list of M group names.
    reference_group_name: str, optional (default: {reference_group_name:})
        Reference group name.
    nsteps: int, optional (default: {nsteps:})
        Overrides number of steps in .mdp file.
    k: float, optional (default: {k:})
        Harmonic pulling spring constant [kJ mol^-1 nm^-2].
    rate: float, optional (default: {rate:})
        Harmonic pulling rate [nm ps^-1].
    mdpin: str, optional (default: {mdpin:})
        GROMACS .mdp parameter input file.
    mdpout: str, optional (default: {mdpout:})
        GROMACS .mdp parameter output file.
    """
    logger = logging.getLogger(__name__)

    assert isinstance(pull_group_names, list), "'pull_group_names' must be list of str."

    logger.debug("Reading mdp file '{mdp:}'...".format(mdp=mdpin))

    # gromacs wrapper parses mdp files
    pull_mdp = gromacs.fileformats.MDP(mdpin)

    if nsteps:
        pull_mdp['nsteps']  = nsteps

    N_pull_coords = len(pull_group_names)

    pull_mdp['pull_ncoords']  = N_pull_coords
    pull_mdp['pull_ngroups']  = N_pull_coords + 1
    pull_mdp['pull_group1_name'] = reference_group_name # the reference group

    for i, n in enumerate(pull_group_names):
        pull_mdp["pull_group{:d}_name".format(i+2)]     = n
        pull_mdp["pull_coord{:d}_type".format(i+1)]     = 'umbrella'  # harmonic potential
        pull_mdp["pull_coord{:d}_geometry".format(i+1)] = 'distance'  # simple distance increase
        pull_mdp["pull_coord{:d}_dim".format(i+1)]      = 'Y Y Y'     # pull in all directions
        pull_mdp["pull_coord{:d}_groups".format(i+1)]   = "1 {:d}".format(i+2) # groups 1 (Chain A) and 2 (Chain B) define the reaction coordinate
        pull_mdp["pull_coord{:d}_start".format(i+1)]    = 'yes'       # define initial COM distance > 0
        pull_mdp["pull_coord{:d}_rate".format(i+1)]     = rate        # 0.1 nm per ps = 10 nm per ns
        pull_mdp["pull_coord{:d}_k".format(i+1)]        = k           # kJ mol^-1 nm^-2

    pull_mdp.write(mdpout)
