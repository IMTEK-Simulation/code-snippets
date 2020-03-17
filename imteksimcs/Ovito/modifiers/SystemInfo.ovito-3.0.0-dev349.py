# SystemInfo.ovito-3.0.0-dev349.py
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
Displays system information and changes the working directory to the directory
containing "SourceFile", usually the last datafile or trajectory loaded in the
course of pipeline execution. Useful for opening files lateron in the pipeline.
"""
from ovito.data import *

def modify(frame, data):

    # an additional attribute holding the source directory of the loaded file:
    if "SourceFile" in data.attributes:
        from os.path import dirname
        data.attributes["SourceDir"] = dirname( data.attributes["SourceFile"] )
        from os import chdir
        chdir(data.attributes["SourceDir"])

    print("There are %i atrributes with the following values:" % len(data.attributes))

    for attribute_name, attribute_value in data.attributes.items():
        print("  '{:24s}: {}'".format(attribute_name, attribute_value))

    print("")

    if data.particles != None:
        print("There are %i particles with the following properties:" % data.particles.count)
        for property_name in data.particles.keys():
            print("  '%s'" % property_name)
