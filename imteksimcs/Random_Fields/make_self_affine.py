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
from scipy.special import gamma

###

def progressbar(x, maxx, len=60):
    xx = x/maxx
    return '|'+int(xx*len)*'#'+((len-int(xx*len))*'-')+'| {:>3}% ' \
           '({}/{})'.format(int(xx*100), x, maxx)

###

def irfftn(karr):
    k = 0
    ndim = len(karr.shape)
    nk = np.sum([np.prod(karr.shape)//i for i in karr.shape])
    for dim in range(ndim):
        # FFT in direction dim
        shape = karr.shape[:dim] + karr.shape[dim+1:]
        for c in np.ndindex(shape):
            sys.stdout.write('{}\r'.format(progressbar(k, nk-1)))
            sys.stdout.flush()
            sl = c[:dim] + (slice(None),) + c[dim:]
            #if dim == ndim-1:
            #    karr[sl] = np.fft.irfft(karr[sl], n=rarr.shape[dim])
            #else:
            karr[sl] = np.fft.ifft(karr[sl])
            k += 1

###

def nsphere(dim):
    """surface area of the dim-dimensional sphere embedded in a space of dim+1"""
    return 2*np.pi**((dim+1)/2)/gamma((dim+1)/2)

def self_affine_prefactor(nb_grid_pts, physical_sizes, Hurst, rms_amplitude=None,
                          rms_derivative=None, short_cutoff=None, long_cutoff=None):
    nb_grid_pts = np.asarray(nb_grid_pts)
    physical_sizes = np.asarray(physical_sizes)

    if short_cutoff is not None:
        q_short = 2*np.pi/short_cutoff
    else:
        q_short = np.pi*np.min(nb_grid_pts/physical_sizes)

    if long_cutoff is not None:
        q_long = 2*np.pi/long_cutoff
    else:
        q_long = 2*np.pi*np.max(1/physical_sizes)

    dim = len(nb_grid_pts)
    volume = np.prod(physical_sizes)
    if rms_amplitude is not None:
        fac = 2*rms_amplitude/np.sqrt(q_long**(-2*Hurst)-q_short**(-2*Hurst))*np.sqrt(Hurst*np.pi)

        #fac = rms_amplitude*np.sqrt(nsphere(dim-1)/(2*Hurst*(2*np.pi)**dim*(q_long**(-2*Hurst)-q_short**(-2*Hurst))))
    elif rms_derivative is not None:
        raise NotImplemented
        #fac = 2*rms_derivative/np.sqrt(q_short**(2-2*Hurst)-q_long**(2-2*Hurst))*np.sqrt((1-Hurst)*np.pi)
    else:
        raise ValueError('Neither rms amplitude nor rms derivative is defined!')
    return fac * np.prod(nb_grid_pts)/np.sqrt(volume)

def Fourier_synthesis(nb_grid_pts, physical_sizes, Hurst, rms_amplitude=None, rms_derivative=None,
                      short_cutoff=None, long_cutoff=None, rolloff=1.0,
                      rfn='rarr.numpy', kfn='karr.numpy'):
    nb_grid_pts = np.asarray(nb_grid_pts)
    physical_sizes = np.asarray(physical_sizes)

    if short_cutoff is not None:
        q_short = 2*np.pi/short_cutoff
    else:
        q_short = np.pi*np.min(nb_grid_pts/physical_sizes)

    if long_cutoff is not None:
        q_long = 2*np.pi/long_cutoff
    else:
        q_long = None

    fac = self_affine_prefactor(nb_grid_pts, physical_sizes, Hurst, rms_amplitude=rms_amplitude,
                                rms_derivative=rms_derivative, short_cutoff=short_cutoff,
                                long_cutoff=long_cutoff)

    ndim = len(nb_grid_pts)
    #knlast = nb_grid_pts[-1]//2+1
    nb_fourier_pts = nb_grid_pts #np.append(nb_grid_pts[:-1], [knlast])

    # Memory mapped arrays
    karr = np.memmap(kfn, np.complex128, 'w+', shape=tuple(nb_fourier_pts))

    print('Creating Fourier representation:')

    wavevectors = [2*np.pi*np.fft.fftfreq(n)*(n/s) for n, s in zip(nb_fourier_pts, physical_sizes)]

    k = 0
    for coord in np.ndindex(tuple(nb_fourier_pts[:-1])):
        sys.stdout.write('{}\r'.format(progressbar(k, np.prod(nb_fourier_pts[:-1])-1)))
        sys.stdout.flush()
        q_sq = sum(wavevectors[i][c]**2 for i, c in enumerate(coord)) + wavevectors[-1]**2
        if all(wavevectors[i][c] == 0 for i, c in enumerate(coord)):
            q_sq[0] = 1
        phase = np.exp(2*np.pi*np.random.rand(nb_fourier_pts[-1])*1j)
        ran = fac * phase*np.random.normal(size=nb_fourier_pts[-1])
        karr[coord] = ran * q_sq**(-(ndim+2*Hurst)/4)
        karr[coord][q_sq > q_short**2] = 0
        if q_long is not None:
            karr[coord][q_sq < q_long**2] = rolloff * ran[mask] * q_long**(-(ndim+2*Hurst)/2)
        k += 1

    karr[(0,)*ndim] = 0 #np.real(karr[(0,)*ndim])

    print('\nInverse FFT:')
    irfftn(karr)
    print()
    return np.real(karr)
