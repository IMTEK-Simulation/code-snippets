# Copyright (c) 2015 Lars Pastewka
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do
# so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
# FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#
# This code uses a Fourier-filtering algorithm to create the self-affine field.
# For the two-dimensional case, the algorithm is briefly described in sections 3 of the following papers:
# T. D. B. Jacobs, T. Junge, and L. Pastewka, Surf. Topogr. Metrol. Prop. 5, 13001 (2017).
# S. B. Ramisetti, C. Campañá, G. Anciaux, J.-F. Molinari, M. H. Müser, and M. O. Robbins, J. Phys. Condens. Matter 23, 215004 (2011).
#

from __future__ import division, print_function

import sys
from argparse import ArgumentParser, ArgumentTypeError

import numpy as np

from imteksimcs.Random_Fields.make_self_affine import Fourier_synthesis


###

def tuplen(s, type=float):
    return tuple(type(x) for x in s.split(','))


###

def main():
    commandline = ' '.join(sys.argv[:])
    parser = ArgumentParser(description='Create a self-affine, randomly rough'
                                        'surface using a Fourier-filtering '
                                        'algorithm.')
    parser.add_argument('filename', metavar='FILENAME',
                        help='name of field output file')
    parser.add_argument('--rms-amplitude', dest='rms_amplitude',
                        type=float,
                        help='create field with root mean square amplitude AMPLITUDE',
                        metavar='AMPLITUDE')
    parser.add_argument('--rms-derivative', dest='rms_derivative',
                        type=float,
                        help='create field with root mean square derivative DERIVATIVE',
                        metavar='DERIVATIVE')
    parser.add_argument('--Hurst', dest='Hurst', type=float, default=0.8,
                        help='self-affine scaling with Hurst exponent HURST',
                        metavar='HURST')
    parser.add_argument('--short-cutoff', dest='short_cutoff', type=float,
                        help='use short wavelength cutoff SHORTCUTOFF',
                        metavar='SHORTCUTOFF')
    parser.add_argument('--long-cutoff', dest='long_cutoff', type=float,
                        help='use long wavelength cutoff LONGCUTOFF',
                        metavar='LONGCUTOFF')
    parser.add_argument('--rolloff', dest='rolloff', type=float, default=1.0,
                        help='rolloff to ROLLOFF * power at LONGCUTOFF',
                        metavar='ROLLOFF')
    parser.add_argument('--nb-grid-pts', dest='nb_grid_pts',
                        type=lambda x: tuplen(x, int), default=(128, 128),
                        help='output field has number of grid points NBGRIDPTS',
                        metavar='NBGRIDPTS')
    parser.add_argument('--size', dest='size',
                        type=tuplen, default=None,
                        help='physical size of output field is SIZE',
                        metavar='SIZE')
    parser.add_argument('--unit', dest='unit',
                        type=str, default='m',
                        help='length unit is UNIT (this information is simply '
                             'dumped to the output file)',
                        metavar='UNIT')
    arguments = parser.parse_args()
    print('filename = {}'.format(arguments.filename))
    print('rms-amplitude = {}'.format(arguments.rms_amplitude))
    print('rms-derivative = {}'.format(arguments.rms_derivative))
    print('Hurst = {}'.format(arguments.Hurst))
    print('short-cutoff = {}'.format(arguments.short_cutoff))
    print('long-cutoff = {}'.format(arguments.long_cutoff))
    print('rolloff = {}'.format(arguments.rolloff))
    print('nb-grid-pts = {}'.format(arguments.nb_grid_pts))
    if arguments.size is None:
        arguments.size = (1.,)*len(arguments.nb_grid_pts)
    print('size = {}'.format(arguments.size))
    print('unit = {}'.format(arguments.unit))

    field = Fourier_synthesis(arguments.nb_grid_pts,
                              arguments.size,
                              arguments.Hurst,
                              rms_amplitude=arguments.rms_amplitude,
                              rms_derivative=arguments.rms_derivative,
                              short_cutoff=arguments.short_cutoff,
                              long_cutoff=arguments.long_cutoff,
                              rolloff=arguments.rolloff)
    field -= np.mean(field)

    print('Field creation done. Measured properties:')
    print('- amplitude:', np.std(field))
    print('- minimum:  ', np.min(field))
    print('- maximum:  ', np.max(field))
    print('- max-min:  ', np.max(field)-np.min(field))

    if len(arguments.nb_grid_pts) == 2:        
        print('Field is two-dimensional. Saving topography to a text file...')
        sx, sy = arguments.size
        np.savetxt(arguments.filename, field, header='{}\nWidth: {} {}\n'
                   'Height: {} {}\nValue units: {}'.format(commandline,
                                                           sx, arguments.unit,
                                                           sy, arguments.unit,
                                                           arguments.unit))
        print('...done')
    else:
        print('Field is {}-dimensional. Saving to a npy file...'.format(len(arguments.nb_grid_pts)))
        np.save(arguments.filename, field)
        print('...done')

if __name__ == '__main__':
    main()
