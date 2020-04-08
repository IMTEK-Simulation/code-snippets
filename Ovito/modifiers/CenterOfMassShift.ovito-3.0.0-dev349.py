# CenterOfMassShift.ovito-3.0.0-dev349.py
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
Displays difference between global "vector" attributes, per defaulf

   CenterOfMass.[X|Y|Z] - OriginalCenterOfMass.[X|Y|Z]

and

   CenterOfMassVelocity.[X|Y|Z] -  OriginalCenterOfMassVelocity.[X|Y|Z]

Use after two CenterOfMass modifiers at different positions in pipeline.
"""

from ovito.data import *

def modify(frame, input, output):
    shift_attribute_names = ["comShift.X","comShift.Y","comShift.Z"]
    com_attribute_names = ["CenterOfMass.X","CenterOfMass.Y","CenterOfMass.Z"]
    previous_com_attribute_names = ["OriginalCenterOfMass.X","OriginalCenterOfMass.Y","OriginalCenterOfMass.Z"]

    vel_shift_attribute_names = ["comVelocityShift.X","comVelocityShift.Y","comVelocityShift.Z"]
    com_vel_attribute_names = ["CenterOfMassVelocity.X","CenterOfMassVelocity.Y","CenterOfMassVelocity.Z"]
    previous_com_vel_attribute_names = ["OriginalCenterOfMassVelocity.X","OriginalCenterOfMassVelocity.Y","OriginalCenterOfMassVelocity.Z"]

    com_shift = []
    try:
        for shift,com,previous_com in zip(shift_attribute_names,com_attribute_names,previous_com_attribute_names):
            output.attributes[shift] = output.attributes[com] - output.attributes[previous_com]
            com_shift.append(output.attributes[shift])
    except:
        print("Warning: error evaluating positional shift.")

    com_vel_shift = []
    try:
        for shift,com_vel,previous_com_vel in zip(vel_shift_attribute_names,com_vel_attribute_names,previous_com_vel_attribute_names):
            output.attributes[shift] = output.attributes[com_vel] - output.attributes[previous_com_vel]
            com_vel_shift.append(output.attributes[shift])
    except:
        print("Warning: error evaluating velocity shift.")

    print("COM positional shift - (X,Y,Z): {}".format( com_shift) )
    print("COM velocity shift   - (X,Y,Z): {}".format( com_vel_shift) )
