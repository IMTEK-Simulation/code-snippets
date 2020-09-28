# -*- coding: utf-8 -*-
#
# mpi_pool_executor.py
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
"""Run an MPI-parallel function via mpi4py's distribution mechanism."""
from mpi4py import MPI
from mpi4py.futures import MPIPoolExecutor
import numpy as np
import logging


def sample_func():
    """Sample function for MPI-parallel processing."""
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    print("Hello from rank {}/{}.".format(rank, size))
    data = (rank+1)**2
    data = comm.gather(data, root=0)
    if rank == 0:
         print("Gathered {}".format(data))
         N = np.sum(data)
         print(f'We got the magic number {N}')
         return N


def call(func, *args, **kwargs):
    """Distribute function call to MPI processes via mpi4py.futures serialization/deserialization mechanisum."""
    logger = logging.getLogger(__name__)
    comm = MPI.COMM_WORLD

    with MPIPoolExecutor() as executor:
        rank = comm.Get_rank()
        size = comm.Get_size()
        version = MPI.Get_version()
        logger.debug("MPI version: %s" % str(version))
        logger.debug("Current MPI size is: %s" % size)
        universe_size = comm.Get_attr(MPI.UNIVERSE_SIZE)
        logger.debug("MPI universe size is: %s" % universe_size)

        # distribute payload to rank 1 to MPI_UNIVERSE_SIZE
        jobs = [executor.submit(work, *args, **kwargs) for _ in range(universe_size-1)]
        # also run payload on rank 0
        data = work(*args, **kwargs)

    # function may have return value on rank 0
    if data is not None:
        return data

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    call(sample_func)