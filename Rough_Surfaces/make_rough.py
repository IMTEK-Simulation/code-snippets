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
# This code uses a Fourier-filtering algorithm to create the rough surfaces.
# The algorithm is briefly described in sections 3 of the following papers:
# T. D. B. Jacobs, T. Junge, and L. Pastewka, Surf. Topogr. Metrol. Prop. 5, 13001 (2017).
# S. B. Ramisetti, C. Campañá, G. Anciaux, J.-F. Molinari, M. H. Müser, and M. O. Robbins, J. Phys. Condens. Matter 23, 215004 (2011).
#

from __future__ import division, print_function

import sys
from argparse import ArgumentParser, ArgumentTypeError

import numpy as np

###

def progressbar(x, maxx, len=60):
    xx = x/maxx
    return '|'+int(xx*len)*'#'+((len-int(xx*len))*'-')+'| {:>3}% ' \
           '({}/{})'.format(int(xx*100), x, maxx)

###

def irfft2(karr, rarr):
    nrows, ncolumns = karr.shape
    for i in range(ncolumns):
        sys.stdout.write('{}\r'.format(progressbar(i, ncolumns+nrows-1)))
        sys.stdout.flush()
        karr[:, i] = np.fft.ifft(karr[:, i])
    for i in range(nrows):
        sys.stdout.write('{}\r'.format(progressbar(i+ncolumns, ncolumns+nrows-1)))
        sys.stdout.flush()
        rarr[i, :] = np.fft.irfft(karr[i, :])

###

def self_affine_prefactor(nx, ny, sx, sy, Hurst, rms_height=None,
                          rms_slope=None, short_cutoff=None, long_cutoff=None):
    if short_cutoff is not None:
        q_max = 2*np.pi/short_cutoff
    else:
        q_max = np.pi*min(nx/sx, ny/sy)

    if long_cutoff is not None:
        q_min = 2*np.pi/long_cutoff
    else:
        q_min = 2*np.pi*max(1/sx, 1/sy)

    area = sx*sy
    if rms_height is not None:
        fac = 2*rms_height/np.sqrt(q_min**(-2*Hurst)-
                                   q_max**(-2*Hurst))*np.sqrt(Hurst*np.pi)
    elif rms_slope is not None:
        fac = 2*rms_slope/np.sqrt(q_max**(2-2*Hurst)-
                                  q_min**(2-2*Hurst))*np.sqrt((1-Hurst)*np.pi)
    else:
        raise ValueError('Neither rms height nor rms slope is defined!')
    return fac * nx*ny/np.sqrt(area)

def Fourier_synthesis(nx, ny, sx, sy, Hurst, rms_height=None, rms_slope=None,
                      short_cutoff=None, long_cutoff=None, rolloff=1.0,
                      rfn='rarr.numpy', kfn='karr.numpy'):
    if short_cutoff is not None:
        q_max = 2*np.pi/short_cutoff
    else:
        q_max = np.pi*min(nx/sx, ny/sy)

    if long_cutoff is not None:
        q_min = 2*np.pi/long_cutoff
    else:
        q_min = None

    fac = self_affine_prefactor(nx, ny, sx, sy, Hurst, rms_height=rms_height,
                                rms_slope=rms_slope, short_cutoff=short_cutoff,
                                long_cutoff=long_cutoff)

    kny = ny//2+1

    # Memory mapped arrays
    rarr = np.memmap(rfn, np.float64, 'w+', shape=(nx, ny))
    karr = np.memmap(kfn, np.complex128, 'w+', shape=(nx, kny))

    print('Creating Fourier representation:')

    qy = 2*np.pi*np.arange(kny)/sy
    for x in range(nx):
        sys.stdout.write('{}\r'.format(progressbar(x, nx-1)))
        sys.stdout.flush()
        if x > nx//2:
            qx = 2*np.pi*(nx-x)/sx
        else:
            qx = 2*np.pi*x/sx
        q_sq = qx**2 + qy**2
        if x == 0:
            q_sq[0] = 1.
        phase = np.exp(2*np.pi*np.random.rand(kny)*1j)
        ran = fac * phase*np.random.normal(size=kny)
        karr[x, :] = ran * q_sq**(-(1+Hurst)/2)
        karr[x, q_sq > q_max**2] = 0.
        if q_min is not None:
            mask = q_sq < q_min**2
            karr[x, mask] = rolloff * ran[mask] * q_min**(-(1+Hurst))
    if nx % 2 == 0:
        karr[0, 0] = np.real(karr[0, 0])
        karr[1:nx//2, 0] = karr[-1:nx//2:-1, 0].conj()
    else:
        karr[0, 0] = np.real(karr[0, 0])
        karr[nx//2, 0] = np.real(karr[nx//2, 0])
        karr[1:nx//2, 0] = karr[-1:nx//2+1:-1, 0].conj()

    #rarr = np.fft.irfft2(karr)
    print('\nInverse FFT:')
    irfft2(karr, rarr)
    print()
    return rarr, karr

###

def tuple2(s, type=float):
    try:
        x, y = (type(x) for x in s.split(','))
        return x, y
    except:
        raise ArgumentTypeError('Size must be sx,sy')

###

if __name__ == '__main__':
    commandline = ' '.join(sys.argv[:])
    parser = ArgumentParser(description='Create a self-affine, randomly rough'
                                        'surface using a Fourier-filtering '
                                        'algorithm.')
    parser.add_argument('filename', metavar='FILENAME',
                        help='name of topography output file')
    parser.add_argument('--rms-height', dest='rms_height',
                        type=float,
                        help='create surface with root mean square height HEIGHT',
                        metavar='HEIGHT')
    parser.add_argument('--rms-slope', dest='rms_slope',
                        type=float,
                        help='create surface with root mean square slope SLOPE',
                        metavar='SLOPE')
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
    parser.add_argument('--resolution', dest='shape',
                        type=lambda x: tuple2(x, int), default=(128, 128),
                        help='output topography map has resolution RESOLUTION',
                        metavar='RESOLUTION')
    parser.add_argument('--size', dest='size',
                        type=tuple2, default=(1., 1.),
                        help='size of surface is SIZE',
                        metavar='SIZE')
    parser.add_argument('--unit', dest='unit',
                        type=str, default='m',
                        help='length unit is UNIT (this information is simply '
                             'dumped to the output file)',
                        metavar='UNIT')
    arguments = parser.parse_args()
    print('filename = {}'.format(arguments.filename))
    print('rms-height = {}'.format(arguments.rms_height))
    print('rms-slope = {}'.format(arguments.rms_slope))
    print('Hurst = {}'.format(arguments.Hurst))
    print('short-cutoff = {}'.format(arguments.short_cutoff))
    print('long-cutoff = {}'.format(arguments.long_cutoff))
    print('rolloff = {}'.format(arguments.rolloff))
    print('resolution = {}'.format(arguments.shape))
    print('size = {}'.format(arguments.size))
    print('unit = {}'.format(arguments.unit))

    nx, ny = arguments.shape
    sx, sy = arguments.size
    surface, surfaceq = Fourier_synthesis(nx, ny, sx, sy, arguments.Hurst,
                                          rms_height=arguments.rms_height,
                                          rms_slope=arguments.rms_slope,
                                          short_cutoff=arguments.short_cutoff,
                                          long_cutoff=arguments.long_cutoff,
                                          rolloff=arguments.rolloff)
    print('Saving outfile topography...')
    np.savetxt(arguments.filename, surface, header='{}\nWidth: {} {}\n'
               'Height: {} {}\nValue units: {}'.format(commandline,
                                                       sx, arguments.unit,
                                                       sy, arguments.unit,
                                                       arguments.unit))
    print('...done')
