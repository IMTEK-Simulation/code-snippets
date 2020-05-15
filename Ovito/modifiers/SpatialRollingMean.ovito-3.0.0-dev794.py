# SpatialRollingMean.ovito-3.0.0-dev794.py
#
# Copyright (C) 2020 IMTEK Simulation
# Authors: Johannes Hoermann, johannes.hoermann@imtek.uni-freiburg.de
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
Sweep slice over selection and compute an elemnt-wise spatial rolling mean...

...on a specific per-atom vector property. Only for visualization purposes.
The total volume covere is determined from minimum/maximum coordinates within
selection +/- tolerance. Expects 'Atomic Volume' property, i.e. from Voronoi
analysis.

For i = 0 .. N particles in total, j = 0 .. M_i slabs that contain particle i,
and k = 0 .. N_j particles within each slice, V_j = sum_k=0 V_k

    <p_i> = (sum_j=0..M sum_k=0..N_j V_k * p_k) / (sum_j=0..M V_j)
"""
from ovito.data import *
import numpy as np

# TODO: make parameters tunable elsewhere

dim = 3  # system dimensionality
slice_thickness = 5.  # thickness of volume slab in sweep direction
positional_step = 1.  # shift of volume slab per evaluation
slice_normal = np.array([0., 0., 1.])  # spatial direction of sweep)
input_property_name = 'Molecular Orientation'
output_property_name = 'Molecular Orientation Rolling Mean'
tolerance = 2.0


def modify(frame, data):
    position = data.particles["Position"]
    atomic_volumes = data.particles["Atomic Volume"]

    # process selection only, otherwise whole system
    if "Selection" in data.particles:
        global_selection = data.particles["Selection"] != 0
    else:
        global_selection = np.ones(data.particles.count, dtype=np.bool)

    prop = data.particles[input_property_name]
    global_prop = prop[global_selection]

    global_natoms = np.count_nonzero(global_selection)

    global_max_pos = np.max(position[global_selection], axis=0) + tolerance
    global_min_pos = np.min(position[global_selection], axis=0) - tolerance

    global_measure = global_max_pos - global_min_pos

    # vector spanning the slice surface:
    slice_surface_diagonal = global_measure * (np.ones(dim) - slice_normal)

    # at given slab thickness and step, that many slices fit into selected vol:
    slice_count = int(np.floor(
        np.dot(global_measure + slice_normal * slice_thickness, slice_normal)
        / positional_step)) + 1

    # per-atom properties for intermmeditae results
    slice_volume_sum = data.particles_.create_property(
        'Slice Volume Sum', dtype=float, components=1)
    local_cumulative_prop = data.particles_.create_property(
        'Local Cumulative Property',
        dtype=float, components=global_prop.shape[1])

    with slice_volume_sum:
        slice_volume_sum[global_selection] = np.zeros(global_natoms)

    msg = """Sweeping selection of {} particles with
        extreme coordinates
            {}
        and
            {}
        by {} slices of {} [length units] thickness
        at steps of {} [length units] in direction ({})""".format(
            global_natoms,
            global_min_pos,
            global_max_pos,
            slice_count,
            slice_thickness,
            positional_step,
            slice_normal)

    print(msg)
    yield msg

    start_pos = global_min_pos - slice_thickness * slice_normal

    step_vec = positional_step * slice_normal

    # sweep "representative volume" and an according selection across system
    for i in range(slice_count):
        yield (i / slice_count)  # ovito progress bar

        cur_min_pos = start_pos + i * step_vec
        cur_max_pos = start_pos + i * step_vec \
            + slice_thickness * slice_normal + slice_surface_diagonal

        print("------------------------------------------")
        print("""slice #{} of {}, spanned between corners
                {}
            and
                {}""".format(i+1, slice_count, cur_min_pos, cur_max_pos))

        selection = (
              np.greater_equal(position, cur_min_pos).all(axis=1)
            & np.less(position, cur_max_pos).all(axis=1)
            & global_selection
        )

        natoms = np.count_nonzero(selection)

        if natoms < 1:
            continue

        # <p_i> = (sum_j=0..M sum_k=0..N_j V_k * p_k) / (sum_j=0..M V_j)
        V_j = np.sum(atomic_volumes[selection])
        p_k_V_k = np.sum(
            np.atleast_2d(atomic_volumes[selection]).T*prop[selection], axis=0)



        print("  #selected particles in slab: {}".format(natoms))
        print("  cumulative property p_k*V_k: {}".format(p_k_V_k))
        # print("  mean property in slab      : {}".format(cumulative_prop))
        print("  slab volume V_j            : {}".format(V_j))
        print("  mean sum_k(p_k*V_k)/V_j    : {}".format(p_k_V_k/V_j))

        # sum up slice stresses:
        mean_prop_outer = np.outer(p_k_V_k, np.ones(natoms))
        with local_cumulative_prop:
            local_cumulative_prop[selection] += mean_prop_outer.T
        with slice_volume_sum:
            slice_volume_sum[selection] += V_j

    local_mean_prop = data.particles_.create_property(
        output_property_name, dtype=float,
        components=global_prop.shape[1])

    # weighted mean for any atom part of one or several overlaps
    # <p_i> = sum(p_{i,j}, j) / sum(V_j, j)
    with local_mean_prop:
        local_mean_prop[global_selection] = (
            local_cumulative_prop[global_selection]
            / np.atleast_2d(slice_volume_sum[global_selection]).T
        )
