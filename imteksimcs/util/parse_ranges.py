#!/usr/bin/env ovitos
#
# parse_ranges.py
#
# Copyright (C) 2018, 2019 IMTEK Simulation
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


def parseRanges(s):
    """ Parses ranges specified in a string into lists:

    Examples
    --------

    In [3]: s = '1, 2, 3, 4, 5, 6'

    In [4]: parseRanges(s)
    Out[4]: [1, 2, 3, 4, 5, 6]

    In [5]: parseRanges('3-6')
    Out[5]: [3, 4, 5]

    In [6]: parseRanges('3-16-3')
    Out[6]: [3, 6, 9, 12, 15]

    In [6]: parseRanges('1,20-22,3-16-3')
    Out[6]: [1, 3, 6, 9, 12, 15, 20, 21, 22]

    Adapted from
        https://stackoverflow.com/questions/4726168/parsing-command-line-input-for-numbers
    """
    result = set()
    for part in s.split(','):
        #x = part.split('-')
        if '-' in part:
            i = part.index('-')
            result.update( range(*map(int, part.split('-'))) )
        else:
            result.add( int(part) )
        #result.update(range(int(x[0]), int(x[-1]) + 1))
    return sorted(result)

    # if '-' in s:
    #     i = s.index('-')
    #     return range(*map(int, s.split('-')))
    # return map(int, s.split(','))
