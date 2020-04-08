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

import numpy as np

from netCDF4 import Dataset


def get_nearest_indices(time, every):
    """
    Given a list of times, return indices that are most closely to evenly
    spaced time intervals of length *every*.
    """

    time = np.asarray(time)
    r = [ ]
    n = int(np.min(time)/every)
    m = int(np.max(time)/every)+1
    last_j = -1
    for i in range(n, m):
        difftime = np.abs(time-i*every)
        j = np.argmin(difftime)
        if difftime[j] < 1.0 and j != last_j:
            r += [ j ]
            last_j = j
    return np.array(r)

###

def strip_fn(fn):
    if len(fn) > 0 and fn[0] == '+':
        return fn[1:]
    return fn

###

def fix_time(time, dtime=None):
    time = time[:].copy()
    if len(time) > 2:
        dtime = time[-1]-time[-2]
    if len(time) > 1 and abs(dtime-(time[1]-time[0])) > 1e-3:
        time = np.array(time)
        time[0] = time[1]-dtime
    return time, dtime


def open_trajs(trajfns, time_var='time', test_var='coordinates', test_tol=1e-6):
    """
    Open multiple NetCDF trajectory files and check that they are in order.
    Remove overlap if files overlap. Returns list of tuples. Tuple contains
        filename, netCDF4 Dataset object, first frame, last frame
    """

    test_index = None
    i = test_var.find('[')
    j = test_var.find(']')
    if i > 0 and j > 0 and i < j:
        test_index = int(test_var[i+1:j])
        test_var = test_var[:i]

    test_tol = np.asarray(test_tol)

    data_f = list(zip(trajfns, map(Dataset, map(strip_fn, trajfns))))
    filtered_data_f = [ ]

    fn2, data2 = data_f[0]
    last_time = None
    dtime = None
    first1 = 0
    for i in range(len(data_f)-1):
        fn1, data1 = data_f[i]
        fn2, data2 = data_f[i+1]

        print('... %s and %s ...' % ( fn1, fn2 ))

        max_maxdiff = np.zeros_like(test_tol)
        min_maxdiff = test_tol+1e12

        if fn2[0] == '+':
            # Skip test and use all frames, including last
            first2 = 0
            last1 = data1.variables[test_var].shape[0]
        else:
            test1 = data1.variables[test_var]
            test2 = data2.variables[test_var]
            if test_var == time_var:
                test1, dummy = fix_time(test1)
                test2, dummy = fix_time(test2)
            if test_index is not None:
                test1 = test1[:, test_index]
                test2 = test2[:, test_index]

            maxdiff = test_tol+1.0
            first2 = -1
            while first2 < min(test2.shape[0]-1, 5) and \
                      np.any(maxdiff > test_tol):
                first2 += 1
                # Last element in previous file
                last1 = test1.shape[0]-1
                # Maximum difference in test variable
                maxdiff = np.abs(test1[last1] - test2[first2])
                if test_tol.shape == ():
                    maxdiff = np.max(maxdiff)
                while last1 >= 0 and np.any(maxdiff > test_tol):
                    #print 'Frame %i of %s does not agree with frame %i of %s ' \
                    #      '(maximum difference %e). Checking frame %i.' % \
                    #      ( last1, fn1, first2, fn2, maxdiff, last1-1 )
                    max_maxdiff = np.maximum(maxdiff, max_maxdiff)
                    min_maxdiff = np.minimum(maxdiff, min_maxdiff)
                    last1 -= 1
                    maxdiff = np.max(np.abs(test1[last1] - test2[first2]))
                    if test_tol.shape == ():
                        maxdiff = np.max(maxdiff)
                max_maxdiff = np.maximum(maxdiff, max_maxdiff)
                min_maxdiff = np.minimum(maxdiff, min_maxdiff)

        # Sanity check. Frame *last1* of previous file should be identical to
        # frame 0 of current file.
        if last1 < 0:
            raise RuntimeError('{} and {} are not consecutive. Minimum '
                               'residual found was {}, maximum residual {}. '
                               'It may help to increase *test_tol*.' \
                               .format(fn1, fn2, min_maxdiff, max_maxdiff))

        # Retrieve time information. If file has no time information number
        # the individual frames consecutively.
        if time_var in data1.variables:
            # Some files have a bug where the first time slot is zero. Fix by
            # assuming constant time offset between frames.
            time1, dtime = fix_time(data1.variables[time_var], dtime=dtime)
        else:
            time1 = np.arange(data1.variables[test_var].shape[0])
        data_slice = slice(first1, last1)
        time = time1[data_slice]
 
        if last_time is not None:
            # These files are consecutive in the test_var, but may not be
            # consecutive in time. Add an offset to the time.
            time += last_time - time[0]
        filtered_data_f += [ ( fn1, data1, data_slice, time ) ]

        # This becomes the last time of the previous file when in the loop
        if last1 >= len(time):
            if len(time) >= 2:
                last_time = time[-1]+(time[-1]-time[-2])*(len(time)-last1+1)
            else:
                if dtime is None:
                    raise RuntimeError('Cannot guess correct time for file'
                                       ' {}.'.format(fn2))
                else:
                    last_time = time[-1]+dtime*(len(time)-last1+1)
        else:
            last_time = time[last1]

        # Start next file where we think it should start
        first1 = first2

    if time_var in data2.variables:
        # Some files have a bug where the first time slot is zero. Fix by
        # assuming constant time offset between frames.
        time, dtime = fix_time(data2.variables[time_var][:], dtime)
    else:
        time = np.arange(data2.variables[test_var].shape[0])

    if last_time is not None:
        # These files are consecutive in the test_var, but may not be
        # consecutive in time. Add an offset to the time.
        time += last_time - time[0]
    filtered_data_f += [ ( fn2, data2, slice(first1, len(time)), time ) ]

    return filtered_data_f

###

def filter_trajs(idata_f, every):
    # Create a list of all frames
    idata_oframe = reduce(lambda a,b: a+b,
                          [ [ ( fn, data, i, time[i] )
                              for i in range(len(time))[data_slice] ]
                            for fn, data, data_slice, time in idata_f ],
                          [])

    # Get indices that corresponds to roughly equally spaced time slots
    i = get_nearest_indices([time for fn, data, i, time in idata_oframe], every)
    if len(i) == 0:
        raise RuntimeError('No frames left after filtering.')
    else:
        idata_oframe = [ idata_oframe[_i] for _i in i ]

    # Consolidate into a list that contains per-file information
    filtered_idata_f = [ ]
    cur_fn, cur_data, cur_slice, cur_time = idata_oframe[0]
    cur_slice = [ ]
    cur_time = [ ]
    for fn, data, i, time in idata_oframe:
        if data is not cur_data:
            filtered_idata_f += [ ( cur_fn, cur_data, np.array(cur_slice),
                                    np.array(cur_time) ) ]
            cur_fn = fn
            cur_data = data
            cur_slice = [ ]
            cur_time = [ ]
        cur_slice += [ i ]
        cur_time += [ time ]
    filtered_idata_f += [ ( cur_fn, cur_data, np.array(cur_slice),
                            np.array(cur_time) ) ]
    return filtered_idata_f
