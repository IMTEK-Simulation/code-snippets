# MeanOrientation.ovito-3.0.0-dev234.py
#
# Copyright (C) 2019 IMTEK Simulation
# Author: Johannes Hoermann, johannes.hoermann@imtek.uni-freiburg.de
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
Computes mean orientation of (a selection of) atoms and writes it to global
attributes 'Mean Quaternion Orientation' and 'Mean Rodrigues Orientation'.

Aplly after per-atom property 'Orientation' has been determined with the
'Polyhedredral Template Matching' modifier.

Largely a copy of the Ovito scripting example
"Example M3: Color mapping to visualize local lattice orientation" at
https://ovito.org/manual_testing/python/introduction/examples.html, thus GPLv3.
"""
from ovito.data import *
import numpy as np

#source: https://www.ovito.org/manual/python/introduction/examples.html#example-visualize-local-lattice-orientation
def modify(frame, input, output):
    orientation = input.particles["Orientation"]
    selection = input.particles["Selection"]

    qs = orientation[np.nonzero(selection)]

    if len(qs.shape) != 2:
        raise RuntimeError("qs must be a 2-dimensional array")

    if qs.shape[1] != 4:
        raise RuntimeError("qs must be a n x 4 dimensional array")

    # Project quaternions into Rodrigues space:
    # rs = (qs.X/qs.W, qs.Y/qs.W, qs.Z/qs.W)
    # Note that the qs.W may be zero for particles for which no lattice
    # orientation could be computed by the PTM modifier.
    rs = np.zeros_like(qs[:,:3])
    np.divide(qs[:,0], qs[:,3], out=rs[:,0], where=qs[:,3] != 0)
    np.divide(qs[:,1], qs[:,3], out=rs[:,1], where=qs[:,3] != 0)
    np.divide(qs[:,2], qs[:,3], out=rs[:,2], where=qs[:,3] != 0)

    # Compute vector lengths rr = norm(rs)
    rr = np.linalg.norm(rs, axis=1)
    rr = np.maximum(rr, 1e-9) # hack

    # Normalize Rodrigues vectors.
    rs[:,0] /= rr
    rs[:,1] /= rr
    rs[:,2] /= rr

    output.attributes["Mean Quaternion Orientation"] = np.mean( qs, axis=0 )
    output.attributes["Mean Rodrigues Orientation"] = np.mean( rs, axis=0 )
    print("Mean Orientation (X,Y,Z,W): {}".format(
        output.attributes["Mean Quaternion Orientation"] ) )
    print("Mean Orientation (x,y,z): {}".format(
        output.attributes["Mean Rodrigues Orientation"] ) )
