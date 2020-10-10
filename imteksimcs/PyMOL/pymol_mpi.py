# -*- coding: utf-8 -*-
#
# pymol_mpi.py
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
"""mpi4py-based embarassingly parallel pymol rendering of GROMACS trajectories."""

import glob
import logging
import os
import re

from mpi4py import MPI
import numpy as np
import pymol2


def sample_pymol_run(pml, struc_file="default.gro", traj_file="default.xtc", out_prefix="default_"):
    pml.cmd.load(struc_file)

    pml.cmd.select('solvent', 'resn SOL')
    pml.cmd.select('counterion', 'resn BR')
    pml.cmd.select('substrate', 'resn AUM')
    # select surfactant, not (solvent or counterion or substrate)
    pml.cmd.select('surfactant', 'resn CTAB')
    pml.cmd.select('hydrogen', 'surfactant and elem H')
    pml.cmd.select('shell', 'solvent within 3.5 of surfactant')

    pml.cmd.bg_color('white')
    pml.cmd.hide('everything', 'solvent')
    pml.cmd.hide('dots', 'counterion')
    pml.cmd.hide('dots', 'substrate')
    pml.cmd.hide('lines', 'surfactant')
    pml.cmd.hide('everything', 'solvent')

    pml.cmd.set_bond('stick_radius', 0.1, 'surfactant')
    pml.cmd.alter('counterion', 'vdw=0.5')
    pml.cmd.alter('surfactant', 'vdw=0.4')
    pml.cmd.alter('hydrogen', 'vdw=0.2')
    pml.cmd.show('spheres', 'surfactant')
    pml.cmd.show('spheres', 'surfactant')
    pml.cmd.show('sticks', 'surfactant')

    pml.cmd.center('surfactant')
    pml.cmd.orient()

    # https://pymolwiki.org/index.php/Defer_builds_mode
    pml.cmd.set('defer_builds_mode', 3)
    pml.cmd.set('async_builds', 1)
    pml.cmd.set('cache_frames', 0)
    pml.cmd.set('ray_trace_frames', 0)

    pml.cmd.load_traj(traj_file, state=1)
    pml.cmd.mpng(out_prefix)


def run_pymol(pml, struc_file="default.gro", traj_file="default.xtc",
              out_prefix="default_", pml_script="default.pml"):
    """Render per-frame png files after custom pml visualization setup script."""
    pml.cmd.load(struc_file)

    # run any external pml script
    pml.cmd.do('@{}'.format(pml_script))

    # https://pymolwiki.org/index.php/Defer_builds_mode
    pml.cmd.set('defer_builds_mode', 3)
    pml.cmd.set('async_builds', 1)
    pml.cmd.set('cache_frames', 0)
    pml.cmd.set('ray_trace_frames', 0)

    pml.cmd.load_traj(traj_file, state=1)
    pml.cmd.mpng(out_prefix)


# PyMOL does not necessarily write out alphanumerically sortable files, i.e.
# "file_1000.png" and "file_20000.png", fix this:
def batch_rename(chunk_files, out_prefix):
    """Rename intermediate output files for alphanumeric sorting."""
    logger = logging.getLogger(__name__)
    nframes = 0
    nframes_per_step = []
    all_frame_files = []
    for chunk_file in chunk_files:
        chunk_file_basename = os.path.basename(chunk_file)
        frame_file_prefix = os.path.splitext(chunk_file_basename)[0] + '_'
        frame_file_glob_pattern = frame_file_prefix + '*.png'
        frame_files = list(glob.glob(frame_file_glob_pattern))
        all_frame_files.append(frame_files)
        nframes_per_step.append(len(frame_files))
        nframes += len(frame_files)

    offsets = [sum(nframes_per_step[:i]) for i in range(len(nframes_per_step))]
    # offsets = [0, *nframes_per_step[:-1]]
    logger.info("Rendered {nframes:}.".format(nframes=nframes))
    logger.debug("Per-step distribution of frames: {}.".format(nframes_per_step))
    logger.debug("Per-step index offsets : {}.".format(offsets))

    frame_file_regex_pattern = r'.*_([0-9]+)[.]png$'
    frame_file_regex = re.compile(frame_file_regex_pattern)
    # field width when paading with zeroes for alphabetical sorting
    width = np.ceil(np.log10(nframes)).astype(int)
    for offset, frame_files in zip(offsets, all_frame_files):
        for f in frame_files:
            logger.debug("Match '{}' against regex '{}'.".format(f, frame_file_regex_pattern))
            res = frame_file_regex.match(f)
            i = int(res.group(1))
            numstr = "{val:0{width:d}d}".format(val=offset+i, width=width)
            g = "{}{}.png".format(out_prefix, numstr)
            logger.debug("Rename {} to {}.".format(f, g))
            os.rename(f, g)


def render_chunks(struc_file='default.gro', chunk_file_glob_pattern='default_*.xtc',
                  out_prefix="default_", pymol_func=run_pymol, **kwargs):
    """Render a GROMACS trajectory efficiently with PyMOL by distributing
       similarly sized chunks across all ranks.

    Parameters
    ----------
        struc_file: str
            GROMACS gro coordinates file
        chunk_file_glob_pattern: str
            glob pattern to match chunked GROMACS trr or xtc trajectory files
            They will be sorted alphabetically and assigned to ranks.
        out_prefix: str
            prefix for output png files
        pymol_func: callable(pml : pymol2.PyMOL, struc_file : str, traj_file : str, out_prefix : str, **kwargs)
            function to render a trajectory (chunk). See run_pymol(...).
        **kwargs: forwarded to pymol_func

    Output
    ------
        png files for each frame. For a trajectory 'traj.xtc' with 3000 frames,
        the output files will follow the naming scheme 'traj_[0-9][0-9][0-9][0-9].png',
        i.e. are numbered from 'traj_0000.png' to 'traj_2999.png'
    """
    comm = MPI.COMM_WORLD

    size = comm.Get_size()
    rank = comm.Get_rank()

    logger = logging.getLogger("%s:rank[%i/%i]" % (__name__, rank, size))

    # sort alphabetically
    chunk_files = list(sorted(glob.glob(chunk_file_glob_pattern)))
    N = len(chunk_files)
    N_serial = np.ceil(N/size).astype(int)

    if rank == 0:

        logger.info("Render {N:} chunks on {r:} ranks.".format(N=N, r=size))
        if N > size:
            # raise ValueError("There are more chunks than ranks!")
            logger.info("Render at most {n:} chunks per rank.".format(n=N_serial))

    for n in range(N_serial):
        cur_chunk = rank*N_serial + n
        if cur_chunk < N:
            chunk_file = chunk_files[cur_chunk]
            print("Render {chunk_file:} on rank {r:}.".format(
                  chunk_file=chunk_file, r=rank))

            chunk_file_basename = os.path.basename(chunk_file)
            frame_file_prefix = os.path.splitext(chunk_file_basename)[0] + '_'

            # TODO: perform in temp dir
            with pymol2.PyMOL() as pml:
                pymol_func(pml,
                           struc_file=struc_file, traj_file=chunk_file,
                           out_prefix=frame_file_prefix, **kwargs)
        else:
            print("Rank {r:} idle.".format(r=rank))

    comm.Barrier()

    if rank == 0:
        batch_rename(chunk_files, out_prefix)
        logger.info("All done.")
