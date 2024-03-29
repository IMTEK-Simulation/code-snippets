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
selected atoms and stores them in a voxel grid 'structure-factor-2d'.
"""

from ovito.data import *
from ovito.vis import *
import numpy as np

from ovito.data import DataCollection
from ovito.pipeline import ModifierInterface, FileSource
from ovito.traits import OvitoObjectTrait
from ovito.vis import VectorVis, ParticlesVis

PRINCIPAL_AXIS_LABELS = ["1st Principal Axis", "2nd Principal Axis", "3rd Principal Axis"]
PRINCIPAL_MOMENT_LABELS = ["1st Principal Moment", "2nd Principal Moment", "3rd Principal Moment"]


class GyrationTensor(ModifierInterface):
    """
    Visualize gyration tensor and compute secondary properties.

    Needs to follow a ClusterAnalysis modifier that computes the gyration tensor and center of masses.

    This modifier transforms the center of mass and gyration tensor into pseudo particles and vectors for display.

    Creates new pseudo particles 'COM' with next higher particle type id to maximum id in system.
    """
    particle_vis_com = OvitoObjectTrait(ParticlesVis,
                                        shape=ParticlesVis.Shape.Sphere)

    vector_vis_1st_principal_axis = OvitoObjectTrait(VectorVis,
                                                     alignment=VectorVis.Alignment.Base,
                                                     flat_shading=False,
                                                     color=(1., 0., 0.),
                                                     width=0.2,
                                                     scaling=5.,
                                                     title='1st Principal Axis of Gyration Tensor')

    vector_vis_2nd_principal_axis = OvitoObjectTrait(VectorVis,
                                                     alignment=VectorVis.Alignment.Base,
                                                     flat_shading=False,
                                                     color=(0., 1., 0.),
                                                     width=0.2,
                                                     scaling=5.,
                                                     title='2nd Principal Axis of Gyration Tensor')

    vector_vis_3rd_principal_axis = OvitoObjectTrait(VectorVis,
                                                     alignment=VectorVis.Alignment.Base,
                                                     flat_shading=False,
                                                     color=(0., 0., 1.),
                                                     width=0.2,
                                                     scaling=5.,
                                                     title='3rd Principal Axis of Gyration Tensor')

    def modify(self, data: DataCollection, **kwargs):
        """
        This modifier visualizes COM and gyration tensors of clusters.


        Parameters
        ----------
        frame : int
            current animation frame number at which the pipeline is being evaluated.
        data: DataCollection
            passed in from the pipeline system.
        """

        vector_vis_principal_axes = [self.vector_vis_1st_principal_axis,
                                     self.vector_vis_2nd_principal_axis,
                                     self.vector_vis_3rd_principal_axis]

        if "Selection" in data.particles:
            selection = data.particles["Selection"]
        else:
            selection = np.ones(data.particles.count)

        print(f"Selection shape: {selection.shape}")

        index_selection = np.nonzero(selection)[0]
        print(f"Index selection shape: {index_selection.shape}")

        # Extract data from previously run CluserAnalysis modifiers
        cluster_table = data.tables["clusters"]

        unique_cluster_ids = cluster_table['Cluster Identifier']
        print(f"unique_cluster_ids shape: {unique_cluster_ids.shape}")

        coms = cluster_table['Center of Mass']
        gyr_tensors = cluster_table['Gyration Tensor']

        com_particle_ids = np.array([data.particles_.add_particle(com) for com in coms], dtype=int)

        data.particles_['Cluster'][com_particle_ids] = unique_cluster_ids

        print(f"gyr_tensors shape: {gyr_tensors.shape}")

        n_particles = data.particles['Position'].shape[0]
        print(f"n_particles: {n_particles}")
        principal_axes = np.zeros((n_particles, 3, 3))
        principal_momenta = np.zeros((n_particles, 3))
        tilt = np.zeros((n_particles, 3))

        # tensor: vector of 6 components [XX, YY, ZZ, XY, XZ, YZ].
        for com_particle_id, gyr_tensor_vector in zip(com_particle_ids, gyr_tensors):
            xx, yy, zz, xy, xz, yz = gyr_tensor_vector
            gyr_tensor_matrix = np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])

            eigenvalues, eigenvectors = np.linalg.eigh(gyr_tensor_matrix)

            # store eigenvector*eigenvalue, not unit vector
            for i in range(3):
                principal_axes[com_particle_id, :, i] = eigenvectors[:, i]
                principal_momenta[com_particle_id, i] = eigenvalues[i]

                # tilt angle cosne of 3rd principal axis with x,y, or z axis
                tilt[com_particle_id, i] = principal_axes[com_particle_id, -1, i]

        print(f"principal_axes shape: {principal_axes.shape}")
        print(f"principal_momenta shape: {principal_momenta.shape}")

        for i, (principal_axis_label, principal_moment_label) in enumerate(
                zip(PRINCIPAL_AXIS_LABELS, PRINCIPAL_MOMENT_LABELS)):
            data.particles_.create_property(principal_axis_label, data=principal_axes[:, :, i]).vis = \
                vector_vis_principal_axes[i]
            data.particles_.create_property(principal_moment_label, data=principal_momenta[:, i])

        # Create new particle type for pseudo COM particles
        type_prop = data.particles['Particle Type']
        type_name_dict = {pt.id: pt.name for pt in type_prop.types}
        print(f"Current particle types: {type_name_dict}")
        com_particle_type_id = np.max(list(type_name_dict.keys())) + 1
        print(f"COM particle type id: {com_particle_type_id}")
        com_particle_type = ParticleType(
            id=com_particle_type_id, name='COM', color=(0.1, 0.9, 0.0), shape=self.particle_vis_com.shape)

        type_prop.types.append(com_particle_type)
        data.particles_['Particle Type'][com_particle_ids] = com_particle_type_id
        data.particles_['Particle Identifier'][com_particle_ids] = com_particle_ids

        # Display center of mass as ellipsoids
        # https://www.ovito.org/manual/advanced_topics/aspherical_particles.html#ellipsoids

        # derive quaternion from principal axes of gyration tensor
        # https://www.ovito.org/forum/topic/quaternion-representation-of-aspherical-particels
        orientation = np.zeros((n_particles, 4))
        # Attention: reference orientation for aspherical particles is along z axis
        # https://www.ovito.org/forum/topic/quaternion-representation-of-aspherical-particels/
        orientation[com_particle_ids, :] = [
            [*np.cross([0, 0, 1], rot_mat[:, -1]), 1 + np.dot([0, 0, 1], rot_mat[:, -1])] for rot_mat in
            principal_axes[com_particle_ids]]

        data.particles_.create_property("Aspherical Shape")
        data.particles_.create_property("Orientation")
        data.particles_['Aspherical Shape'][com_particle_ids] = np.sqrt(principal_momenta[com_particle_ids, :])
        data.particles_['Orientation'][com_particle_ids] = orientation[com_particle_ids, :]

        # tilt cosine (of angle between 3rd principal axis with x, y, and z unit vectors)
        data.particles_.create_property("Tilt Cosine", data=tilt)