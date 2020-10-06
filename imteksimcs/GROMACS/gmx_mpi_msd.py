# -*- coding: utf-8 -*-
#
# gmx_mpi_msd.py
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
"""mpi4py-based embarassingly parallel computation of RMSDs from GROMACS trajectories."""

import MDAnalysis as mda # here used for reading and analyzing gromacs trajectories
import MDAnalysis.analysis.rms as mda_rms
import numpy as np

import datetime
import getpass
import socket

from mpi4py import MPI

import logging

def atom_rmsd(gro, trr, out, atom_name='AU', **kwargs):
    """Computes time resolved rmsd of atom group identified by atom name.

    rmsd(t) = sqrt(1/N*sum_i=1^N w_i*(x_i(t)-x_i^rev)^2)

    see

    https://www.mdanalysis.org/mdanalysis/documentation_pages/analysis/rms.html

    Units in output textfile are default MDAnalysis units.

    https://www.mdanalysis.org/mdanalysis/documentation_pages/units.html

    Parameters
    ----------
        gro: str
            GROMACS gro coordinates file
        trr: str
            GROMACS trr trajectory file with N frames
        out: str
            output text file
        atom_name: str, optional
            defaults: 'AU'
        **kwargs:
            keyword arguments forwarded to  MDAnalysis.analysis.rms.RMSD

    Output
    ------
        out text file contains time [ps] and rmsd [Ang] in column vectors.
    """
    comm = MPI.COMM_WORLD

    size = comm.Get_size()
    rank = comm.Get_rank()

    logger = logging.getLogger("%s:rank[%i/%i]" % (__name__, rank, size))

    mda_trr = mda.Universe(gro, trr)

    atom_group = mda_trr.atoms[mda_trr.atoms.names == atom_name]

    rmsd_atom_group = mda_rms.RMSD(atom_group, ref_frame=0, **kwargs)

    # in the standard case #ranks > #frames,
    N = len(mda_trr.trajectory)
    span = N//size

    # special treatment for more ranks than frames
    if span < 1:  # less frames than ranks
        span = 1

    n1 = rank*span
    n2 = (rank+1)*span

    if rank >= size:  # treatment for rank > size
        n1 = 0
        n2 = 0
        # in this case, just return empty time_resolved_rdf
    elif rank == size-1:  # treatment for last rank if N >= size
        n2 = N

    if n2 > n1:
        logger.info("RMSD for frame %i to %i." % (n1, n2))
        rmsd_atom_group.run(start=n1, stop=n2)
        data = rmsd_atom_group.rmsd[:, 1:]  # time and rmsd in column vectors
    else:
        data = np.array()
    # format of rmsd:
    # rmsdT = rmsd_atom_group.rmsd.T
    # frame = rmsdT[0]
    # time = rmsdT[1]
    # rmsd = rmsdT[2]


    data_list = comm.gather(data, root=0)
    if rank == 0:
        filtered_data_list = [l for l in data_list if len(l) > 0]

        data = np.vstack(filtered_data_list)
        np.savetxt(out, data, fmt='%.8e',
            header='\n'.join((
                '{modulename:s}, {username:s}@{hostname:s}, {timestamp:s}'.format(
                    modulename=__name__,
                    username=getpass.getuser(),
                    hostname=socket.gethostname(),
                    timestamp=str(datetime.datetime.now()),
                ),
                'https://www.mdanalysis.org/mdanalysis/documentation_pages/analysis/rms.html',
                'rmsd(t) = sqrt(1/N*sum_i=1^N w_i*(x_i(t)-x_i^rev)^2)',
                'time [ps], rmsd [Ang]')))
