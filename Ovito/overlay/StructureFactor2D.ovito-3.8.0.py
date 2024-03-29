# overlay/StructureFactor2D.ovito-3.8.0.py
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
Render overlay from Structure Factor 2D voxel grid created by Structure Factor modifier.

Adapted from https://www.ovito.org/docs/3.8.0/python/introduction/examples/overlays/data_plot.html
"""
import matplotlib

matplotlib.use('Agg')  # Optional: Activate 'Agg' backend for off-screen plotting.
import matplotlib.pyplot as plt
import numpy as np
from ovito.qt_compat import QtGui


def render(args):
    """Render overlay from Structure Factor 2D voxel grid created by Structure Factor modifier."""

    data = args.scene.selected_pipeline.compute(args.frame)

    if 'structure-factor-2d' not in data.grids:
        raise RuntimeError('No 2D structure factor data found')

    structure_factor_2d_grid = data.grids["structure-factor-2d"]

    lj = structure_factor_2d_grid.domain[0, 0]
    lk = structure_factor_2d_grid.domain[1, 1]
    nk, nj, _ = structure_factor_2d_grid.shape

    qj = np.linspace(0, lj, nj)
    qk = np.linspace(0, lk, nk)

    print(f"qj, qk shape: {qj.shape}, {qk.shape}")
    Qj, Qk = np.meshgrid(qj, qk)
    print(f"Qj, Qk shape: {Qj.shape}, {Qk.shape}")

    # Reshape the flat stucture factor into mesh
    Sjk = structure_factor_2d_grid['Structure Factor 2D'].array.reshape((nk, nj), order='F')
    print(f"Sjk shape: {Sjk.shape}")

    # Compute plot size in inches (DPI determines label size)
    dpi = 300
    plot_width = 0.6 * args.size[0] / dpi
    plot_height = 0.6 * args.size[1] / dpi

    # Create matplotlib figure:
    fig, ax = plt.subplots(figsize=(plot_width, plot_height), dpi=dpi)
    fig.patch.set_alpha(1.0)
    plt.title('Structure Factor 2D')

    shadeopts = {'cmap': 'viridis', 'shading': 'gouraud'}
    lnrwidth = 0.5

    norm = matplotlib.colors.SymLogNorm(linthresh=lnrwidth, linscale=1,
                                        vmin=Sjk.min(), vmax=Sjk.max(), base=10)

    im = ax.pcolormesh(Qj, Qk, Sjk, norm=norm, **shadeopts)
    ax.set_xlabel(r"qj (1/Ang)")
    ax.set_ylabel("qk (1/Ang)")
    plt.tight_layout()
    fig.colorbar(im, ax=ax)

    # Render figure to an in-memory buffer.
    buf = fig.canvas.print_to_buffer()
    plt.close(fig)

    # Create a QImage from the memory buffer
    res_x, res_y = buf[1]
    img = QtGui.QImage(buf[0], res_x, res_y, QtGui.QImage.Format_RGBA8888)

    # Paint QImage onto viewport canvas
    args.painter.drawImage(0, 0, img)
