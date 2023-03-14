# modifier/StructureFactor2D.ovito-3.8.0.py
#
# Copyright (C) 2023 IMTEK Simulation
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
"""
Computes structure factor Sjk in qj, qk plane from x, y coordinates of all
selected atoms and stores them in tables

   CenterOfMass.[X|Y|Z]

and

   CenterOfMassVelocity.[X|Y|Z]

Values previously stored in these attributes are moved to

   OriginalCenterOfMass.[X|Y|Z]

and

   OriginalCenterOfMassVelocity.[X|Y|Z]
"""

from ovito.data import *
from ovito.vis import *
import numpy as np


def compute_structure_factor_2d(rxy, Lx, Ly, nj=100, nk=100):
    """Compute 2d structure factor from coordinates in rxy.

    Parameters
    ----------
    rxy: N x 2 np.ndarray
        coordinates  of N paricles in the xy-plane
    Lx: float
    Ly: float
        xy plane dimension
    """
    N = rxy.shape[0]

    print(f"Lx, Ly: {Lx, Ly}")

    x = np.linspace(0, Lx / 2, nj, endpoint=True)
    y = np.linspace(0, Ly / 2, nk, endpoint=True)

    qjk = np.array([2. * np.pi * x / Lx, 2. * np.pi * y / Ly]).T

    Qj, Qk = np.meshgrid(qjk[:, 0], qjk[:, 1])
    print(f"qjk shape: {qjk.shape}")
    print(f"Qj, Qk shape: {Qj.shape}, {Qk.shape}")

    Qjk = np.array([Qj, Qk])
    print(f"Qjk shape: {Qjk.shape}")
    Ql = Qjk.reshape(2, nj * nk)
    print(f"Ql shape: {Ql.shape}")

    Ql_dot_rxy = np.dot(rxy, Ql)
    print(f"Ql.rxy shape: {Ql_dot_rxy.shape}")

    Sl = np.square(np.abs(np.sum(np.exp(1.j * Ql_dot_rxy), axis=0))) / N
    print(f"Sl shape: {Sl.shape}")
    Sjk = Sl.reshape(nj, nk)
    print(f"Sjk shape: {Sjk.shape}")

    table = np.array([Qj.reshape(nj * nk), Qk.reshape(nj * nk), Sjk.reshape(nj * nk)]).T
    print(f"qj, qk, Sjk table shape: {table.shape}")

    return table


def modify(frame, data):
    """
    This modifier computes the 2d structure factor from positions of all selected particles in the xy plane.

    Parameters
    ----------
    frame : int
        current animation frame number at which the pipeline is being evaluated.
    data: DataCollection
        passed in from the pipeline system.
    """

    Lx = data.cell[0, 0]
    Ly = data.cell[1, 1]

    if "Selection" in data.particles:
        selection = data.particles["Selection"]
    else:
        selection = np.ones(data.particles.count)

    print(f"Selection shape: {selection.shape}")

    index_selection = np.nonzero(selection)[0]
    print(f"Index selection shape: {index_selection.shape}")
    rxy = data.particles['Position'][index_selection, :2]
    print(f"rxy shape: {rxy.shape}")

    nj = 100
    nk = 100

    qj, qk, sjk = np.split(compute_structure_factor_2d(rxy, Lx, Ly, nj=nj, nk=nk), 3, 1)

    Sjk = sjk.reshape(nj, nk)
    # table = DataTable(identifier='structure-factord-2d', title='Structure factor Sjk', plot_mode=DataTable.PlotMode.NoPlot)
    # table.x = table.create_property('wave vector qj', data=Qj.flatten())
    # table.y = table.create_property('wave vector qk', data=Qk.flatten())
    # table.z = table.create_property('structure factor Sjk', data=Sjk.flatten())

    # Starting with an empty DataCollection:
    dummy_data = DataCollection()

    # Create a new SimulationCell object defining the outer spatial dimensions
    # of the grid and the boundary conditions, and add it to the DataCollection:
    dummy_data.cell = SimulationCell(pbc=(True, True, False),
                                     vis=SimulationCellVis(line_width=0.03),
                                     is2D=True)
    dummy_data.cell_[:, :3] = [[np.pi, 0, 0], [0, np.pi, 0], [0, 0, 1]]

    # Create the VoxelGrid object and give it a unique identifier by which it can be referred to later on.
    # Link the voxel grid to the SimulationCell object created above, which defines its spatial extensions.
    # Specify the shape of the grid, i.e. the number of cells in each spatial direction.
    # Finally, assign a VoxelGridVis visual element to the data object to make the grid visible in the scene.
    grid = VoxelGrid(
        identifier='structure-factor-2d',
        domain=dummy_data.cell,
        shape=(*Sjk.shape,1),
        grid_type=VoxelGrid.GridType.CellData,
        vis=VoxelGridVis(enabled=True, transparency=0.6))

    # Associate a new property with the voxel grid cells and initialize it with the data from the Numpy array.
    # Note that the data must be provided as linear (1-dim.) array with the following type of memory layout:
    # The first grid dimension (x) is the fasted changing index while the third grid dimension (z) is the
    # slowest varying index. In this example, this corresponds to the "Fortran" memory layout of Numpy.
    grid.create_property('Structure Factor 2D', data=Sjk.flatten(order='F'))

    # Insert the VoxelGrid object into the DataCollection.
    data.objects.append(grid)


