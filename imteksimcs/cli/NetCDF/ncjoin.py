#! /usr/bin/env python

# =========================================================================
# Copyright (2011-2015) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your argument) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# =========================================================================

"""
Join a individual NetCDF trajectory files into a single one.
"""

from __future__ import print_function, division

import os
import sys
from argparse import ArgumentParser

import numpy as np

from netCDF4 import Dataset

from imteksimcs.NetCDF.ncjoin import open_trajs, filter_trajs

###

FRAME_DIM = 'frame'
ATOM_DIM = 'atom'

###

def main():

    parser = ArgumentParser(description='Join multiple NetCDF trajectory files '
                                        'into a single one. The code uses some '
                                        'heuristics to determine if files are '
                                        'consecutive.')
    parser.add_argument('filenames', metavar='FN', nargs='+',
                        help='file to concatenate')
    parser.add_argument('-e', '--every', dest='every', type=float,
                        help='copy only frames at times that are multiples of '
                             'EVERY',
                        metavar='EVERY')
    parser.add_argument('-v', '--test_var', dest='test_var', default='coordinates',
                        help="use variable VAR to test if two frames are "
                             "identical (default: 'coordinates')",
                        metavar='VAR')
    parser.add_argument('-t', '--test_tol', dest='test_tol',
                        default='1e-6',
                        help='use tolerance TOL to test if two frames are '
                             'identical',
                        metavar='TOL')
    parser.add_argument('-k', '--netcdf_format', dest='netcdf_format',
                        default='NETCDF3_64BIT',
                        help="use NetCDF format KIND; available formats are "
                             "'NETCDF3_CLASSIC', 'NETCDF3_64BIT' (default), "
                             "'NETCDF4_CLASSIC' and 'NETCDF4'",
                        metavar='KIND')
    parser.add_argument('-x', '--exclude', dest='exclude',
                        nargs='?', const=None, default='id',
                        help="exclude variables EXCLUDE (comman separated list) "
                             "from being written to the output file (default: "
                             "EXCLUDE='id'). Specify --exclude without any "
                             "following argument to exclude nothing.",
                        metavar='EXCLUDE')
    parser.add_argument('-i', '--index', dest='index', default='id',
                        help="variable INDEX contains particle ids (default: "
                             "INDEX='id')",
                        metavar='INDEX')
    parser.add_argument('-o', '--index-offset', dest='index_offset', type=int,
                        default=-1,
                        help='OFFSET will be added to particle id to get zero-based'
                             ' particle index (default: OFFSET=-1)',
                        metavar='OFFSET')
    arguments = parser.parse_args()
    if ',' in arguments.test_tol:
        arguments.test_tol = np.array([float(x) for x in arguments.test_tol.split(',')])
    else:
        arguments.test_tol = float(arguments.test_tol)
    print('every =', arguments.every, ', test_var =', arguments.test_var, \
          ', test_tol =', arguments.test_tol, ', exclude =', arguments.exclude, \
          ', index =', arguments.index, ', index_offset =', arguments.index_offset)


    ### Sanity check

    if os.path.exists('traj.nc'):
        raise RuntimeError('traj.nc exists already.')


    ### Open input files and filter if requested

    print('Opening files and checking file order...')
    idata_f = open_trajs(arguments.filenames, test_var=arguments.test_var,
                         test_tol=arguments.test_tol)
    if arguments.every is not None:
        idata_f = filter_trajs(idata_f, arguments.every)


    # Report which time slots will be used
    print('The final merged trajectory will contains the following time slots...')
    for trajfn, idata, data_slice, time in idata_f:
        print('--- {0}:'.format(trajfn))
        print(time)


    # Create output file
    odata = Dataset('traj.nc', 'w', clobber=False, format=arguments.netcdf_format)


    ### Copy global attributes

    for attr_str in idata_f[0][1].ncattrs():
        print("> creating global attribute '%s'..." % attr_str)
        odata.setncattr(attr_str, idata_f[0][1].getncattr(attr_str))


    ### Prepare exclude list

    exclude_list = set()
    if arguments.exclude is not None:
        exclude_list = exclude_list.union(set(arguments.exclude.split(',')))


    ### Copy everything else

    cursor = 0
    last_data = None
    for trajfn, idata, data_slice, time in idata_f:
        print("Appending '%s' starting at frame %i..." % ( trajfn, cursor ))
        print('File contains %i relevant time slots: ' % len(time), time)

        index = None
        if arguments.index in idata.variables:
            index = idata.variables[arguments.index]
            if len(index.dimensions) != 2 or index.dimensions[0] != FRAME_DIM or \
                index.dimensions[1] != ATOM_DIM:
                raise RuntimeError('*INDEX* variable must have dimensions (frame, '
                                   'atom).')

        for var_str, var in idata.variables.items():
            if var_str in exclude_list:
                continue

            if var_str not in odata.variables:
                if arguments.netcdf_format == 'NETCDF4' and var_str in idata.dimensions:
                    # In NETCDF4 (HDF5) there cannot be dimension and variable of
                    # the same name
                    print("= skipping variable '%s' because there is a dimension " \
                          "of the same name" % var_str)
                else:
                    print("> creating variable '%s'..." % var_str)
                    for dim_str in var.dimensions:
                        if dim_str not in odata.dimensions:
                            print("> creating dimension '%s'......" % dim_str)
                            dim = idata.dimensions[dim_str]
                            if dim.isunlimited():
                                odata.createDimension(dim_str)
                            else:
                                odata.createDimension(dim_str, len(dim))
                    odata.createVariable(var_str, var.dtype, var.dimensions)
                    ovar = odata.variables[var_str]
                    for attr_str in var.ncattrs():
                        print("> creating attribute '%s' of variable '%s'..." % \
                              ( attr_str, var_str ))
                        ovar.setncattr(attr_str, var.getncattr(attr_str))

        for var_str, var in idata.variables.items():
            if var_str in exclude_list:
                continue

            if var_str not in idata.dimensions:
                if var.dimensions[0] == FRAME_DIM:
                    print("Copying variable '%s'..." % var_str)
                    if var_str == 'time':
                        odata.variables[var_str][cursor:] = time
                    else:
                        n = var.shape[0]
                        for oframe, iframe in enumerate(np.arange(n)[data_slice]):
                            sys.stdout.write('=== {0}/{1} -> {2} ===\r' \
                                .format(iframe+1, n, cursor+oframe+1))
                            var_data = np.array(var[iframe])
                            if not np.isfinite(var_data).all():
                                print("Data is nan or inf in variable '{}' at " \
                                      "frame {}.".format(var_str, cursor))
                            # Reorder atoms by index if exists
                            if index is not None and len(var.dimensions) > 1 and \
                                var.dimensions[1] == ATOM_DIM:
                                sorted_var_data = np.zeros_like(var_data)
                                sorted_var_data[index[iframe] + \
                                                arguments.index_offset] = var_data
                                var_data = sorted_var_data
                            odata.variables[var_str][cursor+oframe] = var_data
                else:
                    if not last_data or var_str not in last_data.variables:
                        print("Copying variable '%s'..." % var_str)
                        odata.variables[var_str][:] = var[:]
                    else:
                        print("Checking variable '%s' for consistency across files..."%\
                              var_str)
                        if np.any(last_data.variables[var_str][:] != var[:]):
                            raise RuntimeError("Data for per-file variable '%s' "
                                               "differs in '%s' and '%s'." % 
                                               ( var_str, trajfns[i-1], trajfns[i] ))

        cursor += len(time)
        last_time = time[-1]
        if last_data:
            last_data.close()
        last_data = idata
    last_data.close()
    odata.close()
    print('Successfully wrote {0} frames.'.format(cursor))

if __name__ == '__main__':
  main()

