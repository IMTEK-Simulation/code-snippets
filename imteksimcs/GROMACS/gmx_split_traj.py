# -*- coding: utf-8 -*-
#
# gmx_split_traj.py
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
"""Split GROMACS trajectories into chunks, i.e. for parallel processing.

Requires mpi4py, numpy, GromacsWrapper, and MDAnalysis.

Tested with

    - Python==3.6.8
    - GromacsWrapper==0.8.0
    - MDAnalysis==0.20.1
    - mpi4py==3.0.1
    - numpy==1.15.2

Note: mpi4py is required only for split_traj_by_mpi_ranks(...) in order to
      infer the number of chunks dynamically from the total number of available
      MOI ranks. Nothing parallel happens here.
"""

import logging
import os


def get_traj_length(traj_file="default.xtc"):
    """Get number of frames in GROMACS trajectory."""
    import MDAnalysis as mda
    system = mda.Universe(traj_file)
    return len(system.trajectory)


def get_chunk_length(N, size):
    span = N//size
    if span < 1:
        span = 1
    elif N % span != 0:
        span += 1
    return span


def split_traj_in_chunks(max_nchunks, struc_file='default.gro',
                         traj_file='default.xtc', run_file='default.tpr'):
    """Split the trajectory into as many chunks as possible, but at most
    max_nchunks. All chunks except the last one are equally sized.

    Parameters
    ----------
        max_nchunks: int
            maximum number of chunks
        struc_file: str
            GROMACS gro coordinates file
        traj_file: str
            GROMACS trr or xtc trajectory file with N frames
        run_file: str
            GROMACS tpr compiled run file, needed for topology info.
    Returns
    -------
        int: length of trajectory (number of frames)

    Output
    ------
        trajectory files for each rank. For a trajectory 'traj.xtc' on 96 ranks,
        the output files will follow the naming scheme 'traj_[0-9][0-9].xtc',
        i.e. are numbered from 'traj_00.xtc' to 'traj_95.xtc'
    """
    import numpy as np
    import gromacs

    logger = logging.getLogger(__name__)
    traj_file_basename = os.path.basename(traj_file)
    prefix, ext = os.path.splitext(traj_file_basename)

    out_traj_file_gmx_arg = prefix + '_' + ext

    N = get_traj_length(traj_file)
    span = get_chunk_length(N, max_nchunks)
    nchunks = np.ceil(N/span).astype(int)
    width = np.ceil(np.log10(nchunks)).astype(int)
    logger.info("Split {traj_file:} with a total of {N:} frames into {n:} "
                "chunks of length {span:} each (last chunk exempt).".format(
                    traj_file=traj_file, N=N, n=nchunks, span=span))
    # TODO: perform in temp dir
    gromacs.trjconv(f=traj_file, s=run_file, o=out_traj_file_gmx_arg,
                    pbc='res', ur='compact', nzero=width, split=span, input='0')  # 0 selects whole system

    return N


def split_traj_by_mpi_ranks(serial_dim=4, parallel_dim=None, **kwargs):
    """Split the trajectory for distriubiton on a 2d serial-parallel grid.

    Parameters
    ----------
        serial_dim: int
            desired number of chunks along serial dimension.
        parallel_dim: int, default: MPI.COMM_WORLD.Get_size()
            desired number of chunks along parallel dimension.
        struc_file: str
            GROMACS gro coordinates file
        traj_file: str
            GROMACS trr or xtc trajectory file with N frames
        run_file: str
            GROMACS tpr compiled run file, needed for topology info.

    Returns
    -------
        int: length of trajectory (number of frames)

    Output
    ------
        trajectory files for each rank. For a trajectory 'traj.xtc' on 96 ranks,
        the output files will follow the naming scheme 'traj_[0-9][0-9].xtc',
        i.e. are numbered from 'traj_00.xtc' to 'traj_95.xtc'

    Ideally, the total number of chunks serial_dim*parallel_dim should be a
    divisor of the trajectory's length. Otherwise, there may be less chunks
    than desired."""

    from mpi4py import MPI
    comm = MPI.COMM_WORLD

    size = comm.Get_size()
    rank = comm.Get_rank()

    logger = logging.getLogger("%s:rank[%i/%i]" % (__name__, rank, size))

    if rank == 0:
        if parallel_dim is None:
            parallel_dim = size

        max_nchunks = parallel_dim*serial_dim
        ret = split_traj_in_chunks(max_nchunks, **kwargs)
        logger.info("All done.")
        return ret
