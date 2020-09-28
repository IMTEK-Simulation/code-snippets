# -*- coding: utf-8 -*-
#
# gmx_mpi_rdf.py
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
"""mpi4py-based embarassingly parallel computation of RDFs from GROMACS trajectories."""

import MDAnalysis as mda # here used for reading and analyzing gromacs trajectories
import MDAnalysis.analysis.rdf as mda_rdf
import numpy as np

import datetime
import getpass
import socket

import logging

from mpi4py import MPI


def atom_atom_rdf(
        gro, trr, out, atom_name_a='AU', atom_name_b='S',
        interval=(0.0, 50.0), **kwargs):
    """Computes time resolved rdf between atom groups identified by atom name.

    https://www.mdanalysis.org/docs/documentation_pages/analysis/rdf.html

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
        atom_name_a, atom_name_b: str, optional
            defaults: 'AU' and 'S'
        interval: tuple or list, optional
            inner and outer cutoff for rdf. default (0.0,80.0)
        **kwargs:
            keyword arguments forwarded to  MDAnalysis.analysis.rdf.InterRDF

    Output
    ------
        out text file contains bins (1st data line), rdf (following data lines)

        bins: (M,) np.ndarray, centers of M bins
        rdf: (M,N) np.ndarray, rdf on M bins for N frames
    """
    comm = MPI.COMM_WORLD

    size = comm.Get_size()
    rank = comm.Get_rank()

    logger = logging.getLogger("%s:rank[%i/%i]" % (__name__, rank, size))

    mda_trr = mda.Universe(gro, trr)

    atom_group_a = mda_trr.atoms[mda_trr.atoms.names == atom_name_a]
    atom_group_b = mda_trr.atoms[mda_trr.atoms.names == atom_name_b]

    time_resolved_rdf = []

    # this assignment is problematic if there are more ranks than data
    N = len(mda_trr.trajectory)
    n1 = rank*(N//size)
    n2 = (rank+1)*(N//size)
    if rank == size-1:  # treatment for last rank if N >= size
        n2 = N

    logger.info("RDF for frame %i to %i." % (n1, n2))

    for i in range(n1, n2):
        rdf = mda_rdf.InterRDF(
            atom_group_a, atom_group_b, range=interval, **kwargs)
        rdf.run(start=i, stop=i+1)
        time_resolved_rdf.append(rdf.rdf.copy())

    # bins is the center of a bin, see
    # https://www.mdanalysis.org/docs/_modules/MDAnalysis/analysis/rdf.html
    # self.bins = 0.5 * (edges[:-1] + edges[1:])
    time_resolved_rdf = np.array(time_resolved_rdf)

    # gathers list of arrays at rank 0
    time_resolved_rdf_list = comm.gather(time_resolved_rdf, root=0)
    if rank == 0:
        time_resolved_rdf = np.vstack(time_resolved_rdf_list)
        # write file
        # 1st dim is time (frame), 2nd dim is bin
        bins = rdf.bins.copy()
        data = np.vstack((bins, time_resolved_rdf))
        np.savetxt(out, data, fmt='%.8e',
            header='\n'.join((
                '{modulename:s}, {username:s}@{hostname:s}, {timestamp:s}'.format(
                    modulename=__name__,
                    username=getpass.getuser(),
                    hostname=socket.gethostname(),
                    timestamp=str(datetime.datetime.now()),
                ),
                'https://www.mdanalysis.org/docs/documentation_pages/analysis/rdf.html',
                'g_ab(r)=(N_a N_b)^-1 sum_i=1^N_a sum_j=1^N_b <delta(|r_i-r_j|-r)>',
                'normalized to g_ab(r) -> 1 for r -> infty',
                'first line: bin centers [Ang], following lines: per-frame rdf')))
