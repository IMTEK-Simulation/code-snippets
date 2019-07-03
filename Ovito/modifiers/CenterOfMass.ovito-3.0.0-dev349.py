# CenterOfMass.ovito-3.0.0-dev349.py
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
Computes center of mass position and velocity of all selected atoms
and stores them in global attributes

   CenterOfMass.[X|Y|Z]

and

   CenterOfMassVelocity.[X|Y|Z]

Values previously stored in these attributes are moved to

   OriginalCenterOfMass.[X|Y|Z]

and

   OriginalCenterOfMassVelocity.[X|Y|Z]
"""

from ovito.data import *
import numpy as np

# assumes all atoms to be of uniform weight
def modify(frame, input, output):
    attribute_names = [ "CenterOfMass.X","CenterOfMass.Y","CenterOfMass.Z",
                        "CenterOfMassVelocity.X","CenterOfMassVelocity.Y","CenterOfMassVelocity.Z"]

    previous_attribute_names = ["OriginalCenterOfMass.X","OriginalCenterOfMass.Y","OriginalCenterOfMass.Z",
                                "OriginalCenterOfMassVelocity.X","OriginalCenterOfMassVelocity.Y","OriginalCenterOfMassVelocity.Z"]

    position = input.particles["Position"]
    velocity = input.particles["Velocity"]
    if "Selection" in input.particles:
        selection = input.particles["Selection"]
    else:
        selection = np.ones(input.particles.count)


    com     = np.mean( position[np.nonzero(selection) ] , axis=0 )
    com_vel = np.mean( velocity[np.nonzero(selection) ] , axis=0 )

    for a,p in zip(attribute_names,previous_attribute_names):
        if a in input.attributes:
            output.attributes[p] = input.attributes[a]

    (   output.attributes["CenterOfMass.X"],
        output.attributes["CenterOfMass.Y"],
        output.attributes["CenterOfMass.Z"] ) = com

    (  output.attributes["CenterOfMassVelocity.X"],
        output.attributes["CenterOfMassVelocity.Y"],
        output.attributes["CenterOfMassVelocity.Z"] ) = com_vel

    print("COM position (X,Y,Z): {}".format( com ) )
    print("COM velocity (X,Y,Z): {}".format( com_vel ) )
