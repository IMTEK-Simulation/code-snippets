#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2020
# Author: Wolfram Georg Nöhring, wolfram.noehring@imtek.uni-freiburg.de
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
#
"""Snippet: mpi4py worker pool with workers running Lammps

This script demonstrates how to setup a mpi4py worker pool with workers
running Lammps instances that are persistent across invokations of the
work functions. Workers need then read/create the configuration only once
upon startup and can work on the configuration in memory during work.

The MPI.COMM_WORLD communicator needs to be split. Rank 0 of MPI.COMM_WORLD
manages the job queue. The other ranks are group communicators. Rank 0 in
each group accepts work from rank 0 of MPI.COMM_WORLD and broadcasts it to
the workers. 

See the discussion here for how to setup a
pool with multiple ranks per worker thread:
https://bitbucket.org/mpi4py/mpi4py/issues/108/mpi4pyfutures-and-jobs-that-uses-more-than

Rank 0 in each group does not participate in the Lammps calculation and
will be idle. Therefore, this work sharing construct would only make sense
if both the number of groups and the number of ranks per group is large.

The work is done by the Calculator class, which also holds the Lammps 
instance. Adapt the class methods according to your needs. Currently,
it implements methods for calculating efffective stiffness constants
or compliances for a particular atom if another atom is displaced.

Currently, the following work is hardcoded:

`work_list = [(True, (1, i)) for i in range(2, 10)]`

The workers displace atoms `i` and measure the force on 
atom 1 to calculate stiffness constants. 

Example:

.. code:: bash

    mpirun -np 6 --oversubscribe python ./lammps_pool.py 2

"""
import argparse
import numpy as np
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor
from lammps import lammps
from textwrap import dedent


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "num_work_groups",
        type=int,
        help="Number of work groups. The total number of ranks must be greater or equal to 2*num_work_groups+1, so that there are at least 2 ranks per group, and one queue manager",
    )
    args = parser.parse_args()

    # Commands for initializing the Lammps instance
    # Could also read a restart file
    #  TODO: remove hardcoding
    initialization = dedent(
        """\
    atom_style atomic
    units metal
    boundary p p p
    lattice fcc 3.615
    region box_region block 0 3 0 3 0 3
    create_box 1 box_region
    create_atoms 1 box
    pair_style eam/alloy
    pair_coeff * * Cu_mishin1.eam.alloy Cu
    """
    )
    # Make a list of work for the worker groups Each entry of
    # the list is a tuple, and the first element must be 'True'.
    # The manager of each group can later tell the workers to
    # stop waiting for work by broadcasting (False, None).
    #  TODO: remove hardcoding
    work_list = [(True, (1, i)) for i in range(2, 10)]

    if MPI.COMM_WORLD.size < args.num_work_groups + 1:
        message = (
            f"number of ranks in MPI.COMM_WORLD is {MPI.COMM_WORLD.size}. Need at least {args.num_work_groups+1} "
            f"ranks for {args.num_work_groups} worker groups and one manager"
        )
        raise ValueError(message)
    group = (MPI.COMM_WORLD.rank - 1) % args.num_work_groups
    rank_in_group = (MPI.COMM_WORLD.rank - 1) // args.num_work_groups

    # Setup communicators for the workers
    if (
        0
        < MPI.COMM_WORLD.rank
        < (MPI.COMM_WORLD.size - ((MPI.COMM_WORLD.size - 1) % args.num_work_groups))
        and rank_in_group != 0
    ):
        MPI.COMM_WGROUP = MPI.COMM_WORLD.Split(group, rank_in_group)
        MPI.COMM_WGROUP.name = "comm_group_{}".format(group)
    else:
        MPI.COMM_WGROUP = MPI.COMM_WORLD.Split(MPI.UNDEFINED, MPI.COMM_WORLD.rank)

    if (
        0
        < MPI.COMM_WORLD.rank
        < (MPI.COMM_WORLD.size - ((MPI.COMM_WORLD.size - 1) % args.num_work_groups))
    ):
        MPI.COMM_GROUP = MPI.COMM_WORLD.Split(group, rank_in_group)
        MPI.COMM_GROUP.name = "comm_group_{}".format(group)
    else:
        MPI.COMM_GROUP = MPI.COMM_WORLD.Split(MPI.UNDEFINED, MPI.COMM_WORLD.rank)

    # Setup the job communicator, which will be used by MPICommExecutor to assign work
    # Only one rank per worker group participates
    if MPI.COMM_WORLD.rank == 0 or (MPI.COMM_GROUP and MPI.COMM_GROUP.rank == 0):
        MPI.COMM_JOB = MPI.COMM_WORLD.Split(0, MPI.COMM_WORLD.rank)
    else:
        MPI.COMM_JOB = MPI.COMM_WORLD.Split(MPI.UNDEFINED, MPI.COMM_WORLD.rank)

    # Create workers' Lammps instances
    if MPI.COMM_WORLD.rank != 0 and MPI.COMM_WGROUP:
        calculator = Calculator(comm=MPI.COMM_WGROUP, identifier=group)
        for line in initialization.splitlines():
            calculator.lammps.command(line)

    MPI.COMM_WORLD.Barrier()
    if MPI.COMM_JOB != MPI.COMM_NULL:
        with MPICommExecutor(
            MPI.COMM_JOB, root=0, max_workers=args.num_work_groups
        ) as executor:
            _ = executor.map(_work, work_list)
        if MPI.COMM_WORLD.rank == 0:
            pass
        else:
            MPI.COMM_GROUP.bcast((False, None), root=0)

    elif MPI.COMM_WGROUP != MPI.COMM_NULL:
        still_running = True
        while still_running:
            # synchronise the workers with the "main worker"
            still_running, work = MPI.COMM_GROUP.bcast(None, root=0)
            if work is not None:
                m, n = work
                matrix = calculator.calculate_stiffness_matrix(m, n)
                print(matrix)


def _work(work):
    MPI.COMM_GROUP.bcast(work, root=0)  # on entry, root will broadcast work


class Calculator(object):
    """Lammps calculator
    Holds a Lammps instance and provides methods that should be used by workers.
    """

    def __init__(self, comm, identifier):
        self.lammps = lammps(
            comm=comm,
            cmdargs=[
                "-echo",
                "none",
                "-screen",
                "none",
                "-log",
                f"log.lammps.{identifier}",
            ],
        )
        self._etol = 0
        self._ftol = 1e-6
        self._maxiter = 20000
        self._maxeval = 20000
        self._min_style = "fire"
        self._minimization_command = (
            f"minimize {self._etol} {self._ftol} {self._maxiter} {self._maxeval}"
        )

    def calculate_stiffness_matrix(self, m, n):
        matrix = np.zeros((3, 3), float)
        for j in range(3):
            matrix[:, j] = self._calculate_stiffness_column(m, n, j)
        return matrix

    def _calculate_stiffness_column(self, m, n, j, displacement=1e-5):
        """Calculate the stiffness constants of atom m for a displacement of atom n in j-direction."""
        self.lammps.command(f"group tagged id {m}")
        self.lammps.command(f"compute force tagged reduce sum fx fy fz")
        self.lammps.command(f"thermo 1")
        self.lammps.command(
            f"thermo_style custom etotal c_force[1] c_force[2] c_force[3]"
        )
        vector = np.zeros(3, float)
        vector[j] = displacement
        self.lammps.command(f"group displaced id {n}")
        self.lammps.command(
            f"displace_atoms displaced move {vector[0]} {vector[1]} {vector[2]} units box"
        )
        self.lammps.command(f"run 0")
        f = self.lammps.extract_compute("force", 0, 1)
        force_1 = np.array((f[0], f[1], f[2]))
        vector[j] = -2.0 * displacement
        self.lammps.command(
            f"displace_atoms displaced move {vector[0]} {vector[1]} {vector[2]} units box"
        )
        self.lammps.command(f"run 0")
        f = self.lammps.extract_compute("force", 0, 1)
        force_2 = np.array((f[0], f[1], f[2]))
        tangent = (force_1 - force_2) / 2.0 / displacement
        vector[j] = 1.0 * displacement
        # Cleanup
        self.lammps.command(
            f"displace_atoms displaced move {vector[0]} {vector[1]} {vector[2]} units box"
        )
        self.lammps.command(f"uncompute force")
        self.lammps.command(f"thermo_style custom step etotal")
        self.lammps.command(f"group displaced delete")
        self.lammps.command(f"group tagged delete")
        return tangent

    def calculate_compliance_matrix(self, m, n):
        matrix = np.zeros((3, 3), float)
        for j in range(3):
            matrix[:, j] = self._calculate_compliance_column(m, n, j)
        return matrix

    def _calculate_compliance_column(self, m, n, j, force=1e-3):
        """Calculate compliance constants of atom m 

        Apply a force on atom n in positive and negative
        j-direction and measure the displacement of mafter
        energy minimization. Calculate the tangen (compliance)

        Parameters
        ----------



        Returns 
        -------
        tangent : array-like
            Compliance in the three Cartesian directions.
        """
        self.lammps.command(f"min_style {self._min_style}")
        if self._min_style == "fire":
            self.lammps.command(f"min_modify integrator verlet")
        self.lammps.command(f"group tagged id {m}")
        self.lammps.command(f"group forced id {n}")
        displacements = np.zeros((2, 3), float)
        image_flags = np.asarray(self.lammps.gather_atoms("image", 0, 1))
        positions = np.asarray(self.lammps.gather_atoms("x", 1, 3))
        for i in range(2):
            sign = (-1.0) ** i
            vector = np.zeros(3, float)
            vector[j] = sign * force
            self.lammps.command(
                f"fix addforce forced addforce {vector[0]} {vector[1]} {vector[2]}"
            )
            self.lammps.command(f"compute displacement_atom tagged displace/atom")
            self.lammps.command(
                f"compute displacement tagged reduce sum c_displacement_atom[1] c_displacement_atom[2] c_displacement_atom[3]"
            )
            self.lammps.command(f"thermo 1")
            self.lammps.command(
                f"thermo_style custom etotal fnorm c_displacement[1] c_displacement[2] c_displacement[3]"
            )
            self.lammps.command(self._minimization_command)
            xnew = np.asarray(self.lammps.gather_atoms("x", 1, 3))
            mm = 3 * (m - 1)
            xdiff = xnew[mm + 0] - positions[mm + 0]
            ydiff = xnew[mm + 1] - positions[mm + 1]
            zdiff = xnew[mm + 2] - positions[mm + 2]
            if n == 3:
                print(xdiff, ydiff, zdiff, m, n, j)
            u = self.lammps.extract_compute("displacement", 0, 1)
            # for c in range(3):
            #    displacements[i, c] = u[c]
            displacements[i, 0] = xdiff
            displacements[i, 1] = ydiff
            displacements[i, 2] = zdiff
            self.lammps.command("unfix addforce")
            self.lammps.command("uncompute displacement")
            self.lammps.command("uncompute displacement_atom")
            # Reset positions
            self.lammps.scatter_atoms("x", 1, 3, np.ctypeslib.as_ctypes(positions))
            self.lammps.scatter_atoms(
                "image", 0, 1, np.ctypeslib.as_ctypes(image_flags)
            )
        if n == 3:
            print(displacements, "displacements")
        tangent = (displacements[0, :] - displacements[1, :]) / 2.0 / force
        # Cleanup
        self.lammps.command(f"thermo_style custom step etotal")
        self.lammps.command(f"group forced delete")
        self.lammps.command(f"group tagged delete")
        return tangent

    def __del__(self):
        self.lammps.close()


if __name__ == "__main__":
    main()
