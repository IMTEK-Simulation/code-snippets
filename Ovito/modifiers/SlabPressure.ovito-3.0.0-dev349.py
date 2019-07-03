# SlabPressure.ovito-3.0.0-dev349.py
#
# Copyright (C) 2019 IMTEK Simulation
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
Sweeps representative volume element along cartesian direction, estimates the
true volume V of a material slab contained within each volume based on the
extreme coordinates of all encompassed N atoms, computes "anisotropic pressures"
P_ii and hydrostatic pressure P with respect to each volume element, and assigns
these to the encompassed N atoms. If an atom n is contained in multiple
representative volume elements, then the effective pressure values are assigned
as a mean with each contribution weighted by its respective volume.
To avoid overlap, set positional_step >= slice_thickness.

The "anisotropic pressures" are

    P_ii = - 1/V*sum(S_{n,ii},n={0,...,N})

with i=x,y,z, V the representative volume element, N encompassed atoms, and
S_{n,ij} the per-atom stress tensor expected as written by LAMMPS in array
"c_peratom_stress" with six independent spatial components per frame per atom
(refer to https://lammps.sandia.gov/doc/compute_stress_atom.html) and its
diagonal entries in the first three fields (c_peratom_stress[0:3]). Its unit is
[S] = [P]*[V]. The hydrostatic pressure is

    P = sum(P_ii, i=x,y,z) / dim

with dim = 3 the system's dimensionality.
"""
from ovito.data import *
import numpy as np

def modify(frame, data):
    # TODO: make parameters tunable elsewhere
    dim = 3 # system dimensionality
    slice_thickness = 50. # thickness of volume slab in sweep direction
    positional_step = 50. # shift of volume slab per evaluation
    slice_normal = np.array([0.,0.,1.]) # spatial direction of sweep

    # ignore off-diagonal stresses for now
    peratom_stress = data.particles["c_peratom_stress"][:,0:3]
    position = data.particles["Position"]

    # process selection only, otherwise whole system
    if "Selection" in data.particles:
        global_selection = data.particles["Selection"]
    else:
        global_selection = np.ones(data.particles.count)

    global_peratom_stress = peratom_stress[ np.nonzero(global_selection) ]
    global_position = position[ np.nonzero(global_selection) ]

    global_natoms = np.count_nonzero( global_selection )

    global_max_pos = np.max( position[ np.nonzero(global_selection) ], axis=0 )
    global_min_pos = np.min( position[ np.nonzero(global_selection) ], axis=0 )

    global_measure = global_max_pos - global_min_pos

    # vector spanning the slice surface:
    slice_surface_diagonal = global_measure * ( np.ones(dim) - slice_normal)

    # at given slab thcikness and step, that many slices fit into selected vol:
    slice_count = int ( np.floor(
        np.dot(global_measure + slice_normal * slice_thickness, slice_normal) \
            / positional_step ) ) + 1

    # per-atom properties for intermmeditae results
    slice_overlap_count = data.particles_.create_property('Slice Overlap Count',
        dtype=int, components=1)
    slice_volume_sum    = data.particles_.create_property('Slice Volume Sum',
        dtype=float, components=1)
    local_cumulative_stress_tensor_diagonal = data.particles_.create_property(
        'Local Cumulative Stress', dtype=float, components=3)

    with slice_overlap_count:
        slice_overlap_count[ np.nonzero(global_selection) ] = np.zeros(
            global_natoms )
    with slice_volume_sum:
        slice_volume_sum[ np.nonzero(global_selection) ] = np.zeros(
            global_natoms )

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

    step_vec  = positional_step * slice_normal

    # sweep "representative volume" and an according selection across system
    for i in range(slice_count):
        yield (i / slice_count) # ovito progress bar

        cur_min_pos = start_pos + i * step_vec
        cur_max_pos = start_pos + i * step_vec \
            + slice_thickness * slice_normal + slice_surface_diagonal

        print("------------------------------------------")
        print("""slice #{} of {}, spanned between corners
                {}
            and
                {}""".format(i+1, slice_count, cur_min_pos, cur_max_pos))

        selection = \
                np.greater_equal(
                    position,
                    cur_min_pos ).all(axis=1) \
            &   np.less(
                    position,
                    cur_max_pos ).all(axis=1) \
            &   global_selection

        natoms  = np.count_nonzero( selection )
        print("  #selected particles        : {}".format( natoms ) )

        if natoms < 1:
            continue

        stress  = np.sum( peratom_stress[ np.nonzero(selection) ] , axis=0 )
        max_pos = np.max( position[ np.nonzero(selection) ], axis=0 )
        min_pos = np.min( position[ np.nonzero(selection) ], axis=0 )
        measure = max_pos - min_pos
        volume  = np.product( measure )

        pressure_tensor_diagonal = - stress / volume
        pressure_tensor_trace   = np.sum( pressure_tensor_diagonal ) / dim

        print("  cumulative stress   (X,Y,Z): {}".format( stress ) )
        print("  maximum coordinates (X,Y,Z): {}".format( max_pos ) )
        print("  minimum coordinates (X,Y,Z): {}".format( min_pos ) )
        print("  slab measures       (X,Y,Z): {}".format( measure ) )
        print("  slab volume                : {}".format( volume ) )
        print("  pressure tensor diagonal   : {}".format(
            pressure_tensor_diagonal ) )
        print("  pressure tensor trace      : {}".format(
            pressure_tensor_trace ) )

        # sum up slice stresses:
        stress_tensor_diagonal_outer = np.outer( stress, np.ones( natoms ) )
        with local_cumulative_stress_tensor_diagonal:
            local_cumulative_stress_tensor_diagonal[ np.nonzero(selection) ] \
                += stress_tensor_diagonal_outer.T
        with slice_overlap_count:
            slice_overlap_count[ np.nonzero(selection) ] += 1
        with slice_volume_sum:
            slice_volume_sum[ np.nonzero(selection) ] += volume

    local_mean_pressure_tensor_diagonal = \
        data.particles_.create_property('Local Mean Pressure Tensor Diagonal',
            dtype=float, components=3)
    local_mean_pressure_tensor_trace = \
        data.particles_.create_property('Local Mean Pressure Tensor Trace',
            dtype=float, components=1)

    # in case of overlap between volume elements:
    # weighted mean for any atom part of one or several overlaps
    # <p_ii> = sum(p_{ii,j} * V_j, j) / sum(3 V_j, j)
    with local_mean_pressure_tensor_diagonal:
        local_mean_pressure_tensor_diagonal[ np.nonzero(global_selection) ] = \
            - local_cumulative_stress_tensor_diagonal[
                    np.nonzero(global_selection) ] / \
                np.atleast_2d( dim * slice_volume_sum[
                    np.nonzero(global_selection) ] ).T
    with local_mean_pressure_tensor_trace:
        local_mean_pressure_tensor_trace[ np.nonzero(global_selection) ] = \
            np.sum( local_mean_pressure_tensor_diagonal[ \
                np.nonzero(global_selection) ], axis=1 )
