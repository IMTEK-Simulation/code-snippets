# H2OMolecularOrientation.ovito-3.0.0-dev794.py
#
# Copyright (C) 2020 IMTEK Simulation
# Authors: Christian Seidl, Johannes Hoermann, johannes.hoermann@imtek.uni-freiburg.de
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
Compute molecular orientations of selected water molecules as unit vectors
pointing along the HOH angular bisector from O towards HH. Stored as per-atom
property 'Molecular Orientation', same for all three atoms within an H2O
molecule. Set 'Molecular Orientation' to [0,0,0] for all other unaffected atoms.

Set oxygen and hydrogen atom types in modifier header.
"""
from ovito.data import *
import numpy as np

o_type = 7
h_type = 8

def modify(frame, data):

    if "Selection" in data.particles:
        selection = data.particles["Selection"] != 0
    else:
        selection = np.ones(data.particles.count, dtype=np.bool)

    molecules_in_selection = data.particles['Molecule Identifier'][selection]
    unique_molecule_ids = np.unique(molecules_in_selection)

    count = len(unique_molecule_ids)
    print("There are %i molecules within the selection." % count)

    molecular_orientation = np.zeros((data.particles.count, 3))

    o_sel = (data.particles['Particle Type'] == o_type) & selection
    o_sort = np.argsort(data.particles['Molecule Identifier'][o_sel])
    h_sel = (data.particles['Particle Type'] == h_type) & selection
    h_sort = np.argsort(data.particles['Molecule Identifier'][h_sel])
    # print(data.particles['Molecule Identifier'][h_sel][h_sort])
    # print(data.particles['Molecule Identifier'][h_sel][h_sort].reshape(count,2))

    o = data.particles['Position'][o_sel][o_sort]
    assert len(o) == count, "Each molecule must have exactly one oxygen atom."

    h = data.particles['Position'][h_sel][h_sort]
    assert len(h) == 2*count, "Each molecule must have exactly two hydrogen atoms."
    o = np.repeat(o, 2, axis=0)
    flat_vecs = -o + h
    vecs = flat_vecs.reshape((count,2,3))
    vec = vecs[:,0] + vecs[:,1]
    vec_norm = np.linalg.norm(vec, axis=1)
    unit_vec = vec.T / vec_norm

    rev_o_sort = np.argsort(o_sort)
    rev_h_sort = np.argsort(h_sort)
    molecular_orientation[o_sel] = unit_vec.T[rev_o_sort]
    molecular_orientation[h_sel] = np.repeat(unit_vec.T, 2, axis=0)[rev_h_sort]

    data.particles_.create_property('Molecular Orientation', data=molecular_orientation)
