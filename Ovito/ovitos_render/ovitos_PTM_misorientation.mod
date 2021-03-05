
# Ovito Script 3.x syntax
 
# Import OVITO modules.
from ovito.io import *
from ovito.modifiers import *
from ovito.modifiers import PythonScriptModifier
from ovito.data import *
from ovito.vis import *
from ovito import dataset

# Import PYTHON modules
import math
import numpy as np

from fractions import gcd

import os
import sys

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import NullFormatter 

import inspect

def lineno():
    """Returns the current line number in our program."""
    return inspect.currentframe().f_back.f_lineno

def progressbarstr(x, maxx, len=60):
    xx = x/maxx
    if xx != 1:
      return '|'+int(xx*len)*'#'+((len-int(xx*len))*'-')+'| {:>3}% ' \
             '({}/{})'.format(int(xx*100), x, maxx)
    else:
      return '|'+int(xx*len)*'#'+((len-int(xx*len))*'-')+'| {:>3}% ' \
             '({}/{})\n'.format(int(xx*100), x, maxx)

def progressbar(x, maxx, len=60):
    sys.stdout.write('{}\r'.format(progressbarstr(x, maxx, len=len)))
    sys.stdout.flush()

def quaternions_to_orientation_matrix(qs):
    """ Takes a list of quaternions (Nx4 array) and returns corresponding orientation matrices """
    
    if len(qs.shape) != 2:
        raise RuntimeError("qs must be a 2-dimensional array")
 
    if qs.shape[1] != 4:
        raise RuntimeError("qs must be a n x 4 dimensional array")
    
    # From PTM algo :
    # qs[:,0] = qs.X 
    # qs[:,1] = qs.Y 
    # qs[:,2] = qs.Z 
    # qs[:,3] = qs.W 
    
    # Set the quaternion as q = qr + i*qi + j*qj + k*qk
    q = np.zeros_like(qs)
    q[:,[1,2,3,0]] = qs[:,[0,1,2,3]]

    qr, qi, qj, qk = q[:,0], q[:,1], q[:,2], q[:,3]
    # q0, q1, q2, q3 = q[:,0], q[:,1], q[:,2], q[:,3]
    
    s = (np.sqrt(qr**2 + qi**2 + qj**2 + qk**2))**-2

    om = np.zeros([qs.shape[0],3,3], dtype=float)
    
    om[:,0,0] = 1 - 2 * s * (qj**2 + qk**2)
    om[:,1,1] = 1 - 2 * s * (qi**2 + qk**2)
    om[:,2,2] = 1 - 2 * s * (qi**2 + qj**2)

    om[:,0,1] = 2 * s * (qi * qj - qk * qr)
    om[:,0,2] = 2 * s * (qi * qk + qj * qr)

    om[:,1,0] = 2 * s * (qi * qj + qk * qr)
    om[:,1,2] = 2 * s * (qj * qk - qi * qr)

    om[:,2,0] = 2 * s * (qi * qk - qj * qr)
    om[:,2,1] = 2 * s * (qj * qk + qi * qr)
    
    # P = -1
    # q_ = q0**2 - (q1**2 + q2**2 + q3**2)
    
    # om[:,0,0] = q_ + 2*q1**2
    # om[:,1,1] = q_ + 2*q2**2
    # om[:,2,2] = q_ + 2*q3**2

    # om[:,0,1] = 2 * (q1*q2 - P*q0*q3)
    # om[:,0,2] = 2 * (q1*q3 + P*q0*q2)

    # om[:,1,0] = 2 * (q1*q2 + P*q0*q3)
    # om[:,1,2] = 2 * (q2*q3 - P*q0*q1)

    # om[:,2,0] = 2 * (q1*q3 - P*q0*q2)
    # om[:,2,1] = 2 * (q2*q3 - P*q0*q1)
    
    
    # print(om[0,::])
    # print()
    # print(np.linalg.inv(om)[0,:,:])
    
    return om
    
def apply_symmetry(position,space_group):
    if space_group == 195:
        # set of matrices
        R1 = np.array([[1,0,0],
                       [0,1,0],
                       [0,0,1]])
        R2 = np.array([[-1,0,0],
                       [0,-1,0],
                       [0,0,1]])
        R3 = np.array([[-1,0,0],
                       [0,1,0],
                       [0,0,-1]])
        R4 = np.array([[1,0,0],
                       [0,-1,0],
                       [0,0,-1]])

        R5 = np.array([[0,0,1],
                       [1,0,0],
                       [0,1,0]])
        R6 = np.array([[0,0,1],
                       [-1,0,0],
                       [0,-1,0]])
        R7 = np.array([[0,0,-1],
                       [-1,0,0],
                       [0,1,0]])
        R8 = np.array([[0,0,-1],
                       [1,0,0],
                       [0,-1,0]])
        
        R9 = np.array([[0,1,0],
                       [0,0,1],
                       [1,0,0]])
        R10 = np.array([[0,-1,0],
                       [0,0,1],
                       [-1,0,0]])
        R11 = np.array([[0,1,0],
                       [0,0,-1],
                       [-1,0,0]])
        R12 = np.array([[0,-1,0],
                       [0,0,-1],
                       [1,0,0]])
                       
        
        general_symmetry = np.array([R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11,R12])
        
        assert(general_symmetry.shape[0] == 12)
        
        # opp_general_symmetry = - general_symmetry
        # stack = np.vstack((general_symmetry,opp_general_symmetry))
        
    elif space_group == 221:
        # Space Group 221 (Pm-3m)
        R1 = np.array([[1,0,0],
                       [0,1,0],
                       [0,0,1]])
        R2 = np.array([[-1,0,0],
                       [0,-1,0],
                       [0,0,1]])
        R3 = np.array([[-1,0,0],
                       [0,1,0],
                       [0,0,-1]])
        R4 = np.array([[1,0,0],
                       [0,-1,0],
                       [0,0,-1]])
                       
        R5 = np.array([[0,0,1],
                       [1,0,0],
                       [0,1,0]])
        R6 = np.array([[0,0,1],
                       [-1,0,0],
                       [0,-1,0]])
        R7 = np.array([[0,0,-1],
                       [-1,0,0],
                       [0,1,0]])
        R8 = np.array([[0,0,-1],
                       [1,0,0],
                       [0,-1,0]])
                       
        R9 = np.array([[0,1,0],
                       [0,0,1],
                       [1,0,0]])
        R10 = np.array([[0,-1,0],
                       [0,0,1],
                       [-1,0,0]])
        R11 = np.array([[0,1,0],
                       [0,0,-1],
                       [-1,0,0]])
        R12 = np.array([[0,-1,0],
                       [0,0,-1],
                       [1,0,0]])
                       
        R13 = np.array([[0,1,0],
                       [1,0,0],
                       [0,0,-1]])
        R14 = np.array([[0,-1,0],
                       [-1,0,0],
                       [0,0,-1]])
        R15 = np.array([[0,1,0],
                       [-1,0,0],
                       [0,0,1]])
        R16 = np.array([[0,-1,0],
                       [1,0,0],
                       [0,0,1]])
                       
        R17 = np.array([[1,0,0],
                       [0,0,1],
                       [0,-1,0]])
        R18 = np.array([[-1,0,0],
                       [0,0,1],
                       [0,1,0]])                       
        R19 = np.array([[-1,0,0],
                       [0,0,-1],
                       [0,-1,0]])                       
        R20 = np.array([[1,0,0],
                       [0,0,-1],
                       [0,1,0]])                       
                       
        R21 = np.array([[0,0,1],
                       [0,1,0],
                       [-1,0,0]])                       
        R22 = np.array([[0,0,1],
                       [0,-1,0],
                       [1,0,0]])                       
        R23 = np.array([[0,0,-1],
                       [0,1,0],
                       [1,0,0]])                       
        R24 = np.array([[0,0,-1],
                       [0,-1,0],
                       [-1,0,0]])                       
                       
        R25 = np.array([[-1,0,0],
                       [0,-1,0],
                       [0,0,-1]])                       
        R26 = np.array([[1,0,0],
                       [0,1,0],
                       [0,0,-1]])                       
        R27 = np.array([[1,0,0],
                       [0,-1,0],
                       [0,0,1]])                      
        R28 = np.array([[-1,0,0],
                       [0,1,0],
                       [0,0,1]])                       
                       
        R29 = np.array([[0,0,-1],
                       [-1,0,0],
                       [0,-1,0]])                       
        R30 = np.array([[0,0,-1],
                       [1,0,0],
                       [0,1,0]])                       
        R31 = np.array([[0,0,1],
                       [1,0,0],
                       [0,-1,0]])                       
        R32 = np.array([[0,0,1],
                       [-1,0,0],
                       [0,1,0]])                       
                       
        R33 = np.array([[0,-1,0],
                       [0,0,-1],
                       [-1,0,0]])                       
        R34 = np.array([[0,1,0],
                       [0,0,-1],
                       [1,0,0]])                       
        R35 = np.array([[0,-1,0],
                       [0,0,1],
                       [1,0,0]])                       
        R36 = np.array([[0,1,0],
                       [0,0,1],
                       [-1,0,0]])                       
                       
        R37 = np.array([[0,-1,0],
                       [-1,0,0],
                       [0,0,1]])                       
        R38 = np.array([[0,1,0],
                       [1,0,0],
                       [0,0,1]])
        R39 = np.array([[0,-1,0],
                       [1,0,0],
                       [0,0,-1]])
        R40 = np.array([[0,1,0],
                       [-1,0,0],
                       [0,0,-1]])

        R41 = np.array([[-1,0,0],
                       [0,0,-1],
                       [0,1,0]])
        R42 = np.array([[1,0,0],
                       [0,0,-1],
                       [0,-1,0]])        
        R43 = np.array([[1,0,0],
                       [0,0,1],
                       [0,1,0]])
        R44 = np.array([[-1,0,0],
                       [0,0,1],
                       [0,-1,0]])

        R45 = np.array([[0,0,-1],
                       [0,-1,0],
                       [1,0,0]])
        R46 = np.array([[0,0,-1],
                       [0,1,0],
                       [-1,0,0]])
        R47 = np.array([[0,0,1],
                       [0,-1,0],
                       [-1,0,0]])
        R48 = np.array([[0,0,1],
                       [0,1,0],
                       [1,0,0]])

        general_symmetry = np.array([R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11,R12,
                                     R13,R14,R15,R16,R17,R18,R19,R20,R21,R22,R23,R24,
                                     R25,R26,R27,R28,R29,R30,R31,R32,R33,R34,R35,R36,
                                     R37,R38,R39,R40,R41,R42,R43,R44,R45,R46,R47,R48])

        assert(general_symmetry.shape[0] == 48)
                       
    dot_prod = np.zeros((position.shape[0],general_symmetry.shape[0],3))
    for i, row in enumerate(position):
        for j, R in enumerate(general_symmetry):
            dot_prod[i,j,:] = np.dot(R,row)
    
    return dot_prod.reshape(-1, dot_prod.shape[-1])

        
def stereographic_projection(poles, type = "equal_area", bin_and_reduce = False, SST = False):
    norm = np.linalg.norm(poles,axis=-1)
            
    poles =  (poles.T / norm.T).T
    
    # Transpose pole to sample ref. frame
    g = np.array([[0,0,1],
                  [0,1,0],
                  [1,0,0]])
    
    # Spherical angles from uvw
    u,v,w = poles[:,0],poles[:,1],poles[:,2]
        
    alpha = np.arccos(w)
    phi = np.arctan2(v,u) - np.pi/2             # -pi/2 require to orient the poles in the standard direction
    
    alpha[alpha > np.pi/2.] = np.pi - alpha[alpha > np.pi/2.]
        
    if type == "stereographic":
        X = np.tan(alpha/2) * np.cos(phi)
        Y = np.tan(alpha/2) * np.sin(phi)
        # do not consider lower hemisphere
        X[alpha > np.pi/2.] = np.inf
        Y[alpha > np.pi/2.] = np.inf
    elif type == "equal_area":
        r = np.sqrt(2*(1-np.cos(alpha)))
        X = r * np.cos(phi)
        Y = r * np.sin(phi)
        # do not consider lower hemisphere
        # X[alpha > np.pi/2.] = np.inf
        # Y[alpha > np.pi/2.] = np.inf
        
        
    if bin_and_reduce: # Return pole density to plot contour map 
        N = 200
               
        if SST:
            x_edges, d_x = np.linspace(0, stereographic_projection(poles_to_SST(np.array([[0,1,1]]),SG = 221))[0][0], N, retstep=True)
            y_edges, d_y = np.linspace(0, stereographic_projection(poles_to_SST(np.array([[1,1,1]]),SG = 221))[0][1], N, retstep=True)
        else:
            x_edges, d_x = np.linspace(np.sqrt(2)*np.cos(np.pi), np.sqrt(2)*np.cos(0), N, retstep=True)
            y_edges, d_y = np.linspace(np.sqrt(2)*np.sin(- np.pi/2.), np.sqrt(2)*np.sin(np.pi/2.), N, retstep=True)
               
        H, x_edges, y_edges = np.histogram2d(X, Y, bins=(x_edges, y_edges))
        x_centers = (x_edges[1:] + x_edges[:-1]) /2.
        y_centers = (y_edges[1:] + y_edges[:-1]) /2.
        
        X, Y = np.meshgrid(x_centers, y_centers)
        
        mask2d = np.sqrt(X[:,:]**2+Y[:,:]**2) >= np.sqrt(2)
        
        X, Y, H = np.ma.array(X,mask=mask2d), np.ma.array(Y,mask=mask2d), np.ma.array(H,mask=mask2d)
 
        return X,Y,H.T # Return pole density after stereographic projection in order to plot contour map 
    
    return np.array([X, Y]).T # Return pole Cartesian coordinates after stereographic projection
    
def equivalent_poles(poles, SG):
    """ Takes a list of poles (Nx3 array) and returns a list of ALL equivalent directions
    (Nx3 array) depending on the given symmetry of the space group SG
    """
    equivalent_directions = apply_symmetry(poles,SG)
    unique_poles =  np.vstack({tuple(row) for row in equivalent_directions}) 

    return unique_poles
        
def poles_to_SST(poles, SG):
    """ Takes a list of poles (Nx3 array) and returns the corresponding
    list of minimal poles (Nx3 array) lying in the standard stereographic triangle
    """

    if len(poles.shape) != 2:
        raise RuntimeError("poles must be a 2-dimensional array")
 
    if poles.shape[1] != 3:
        raise RuntimeError("poles must be a n x 3 dimensional array")
        
    # Get all equivalent directions from the symmetry of the space group 'SG'
    equivalent_directions = apply_symmetry(poles,SG)
            
    # Filter poles belonging to the standard stereographic triangle (SST), minimal symmetry zone for cubic symmetry.
    mask_SST = np.logical_and(equivalent_directions[:,2]-equivalent_directions[:,1]>=0,np.logical_and(equivalent_directions[:,1]+equivalent_directions[:,0]>=0,-equivalent_directions[:,0]>=0))
            
    min_poles =  np.vstack({tuple(row) for row in equivalent_directions[mask_SST]}) 
   
    return min_poles

    

def poles_to_RGB(poles):
    """ Takes a list of poles (Nx3 array) and returns a list of
        corresponding RGB colors (Nx3 array) """
        
    if len(poles.shape) != 2:
        raise RuntimeError("poles must be a 2-dimensional array")
 
    if poles.shape[1] != 3:
        raise RuntimeError("poles must be a n x 3 dimensional array")
    
    rgb = np.zeros_like(poles)
        
    h,k, l = poles.T 
    
    rgb[:,0] = l-k
    rgb[:,1] = k+h
    rgb[:,2] = -h
            
    # Normalize values.
    max = np.amax(rgb, axis = 1)
    rgb[:,0] = rgb[:,0] / max
    rgb[:,1] = rgb[:,1] / max
    rgb[:,2] = rgb[:,2] / max
     
    return rgb

    
def HSV_to_RGB(HSV):
    """
    Converts from *HSV* colourspace to *RGB* colourspace.

    Parameters
    ----------
    HSV : array_like
        *HSV* colourspace array.

    Returns
    -------
    ndarray
        *RGB* colourspace array.

    Notes
    -----
    -   Input *HSV* colourspace array is in domain [0, 1].
    -   Output *RGB* colourspace array is in range [0, 1].

    References
    ----------
    -   :cite:`EasyRGBn`
    -   :cite:`Smith1978b`
    -   :cite:`Wikipediacg`

    Examples
    --------
    >>> HSV = np.array([0.27867384, 0.74400000, 0.98039216])
    >>> HSV_to_RGB(HSV)  # doctest: +ELLIPSIS
    array([ 0.4901960...,  0.9803921...,  0.2509803...])
    """

    H, S, V = HSV[:,0]/(2*np.pi),HSV[:,1]/np.amax(HSV[:,1]),HSV[:,2]

    h = np.asarray(H * 6)
    h[np.asarray(h == 6)] = 0

    i = np.floor(h)
    j = V * (1 - S)
    k = V * (1 - S * (h - i))
    l = V * (1 - S * (1 - (h - i)))  # noqa

    i = np.array((i, i, i)).astype(np.uint8)

    RGB = np.choose(
        i, [
            np.array([V, l, j]),
            np.array([k, V, j]),
            np.array([j, V, l]),
            np.array([j, k, V]),
            np.array([l, j, V]),
            np.array([V, j, k]),
        ],
        mode='clip')

    return RGB.T
    

def plot_IPF(orientation_matrix):
    om = orientation_matrix
    
    # Main poles forming the standard stereographic triangle (reduced symmetry zone)
    familly_to_plot = np.array([[0,0,1],
                              [0,1,5],
                              [0,1,3],
                              [0,1,2],
                              [0,2,3],
                              [0,1,1],
                              [-1,3,3],
                              [-1,2,2],
                              [-1,1,1],
                              [-1,1,2],
                              [-1,1,3],
                              [-1,1,5],
                              [-1,2,3]])
                              
    ax=plt.subplot(221)

    ax.set_aspect(1)

    for at in range(om.shape[0]):
        progressbar(at,om.shape[0]) 
        spx = stereographic_projection(om[at,0,:])
        spy = stereographic_projection(om[at,1,:])
        spz = stereographic_projection(om[at,2,:])

        ax.scatter(spx[0],spx[1],color="r",s=50)
        ax.scatter(spy[0],spy[1],color="g",s=50)
        ax.scatter(spz[0],spz[1],color="b",s=50)
    
    legend_elements = [Line2D([0], [0], marker='o', color='r', label='X-projection',
                          markerfacecolor='r', markersize=5),
                      Line2D([0], [0], marker='o', color='g', label='Y-projection',
                          markerfacecolor='g', markersize=5),
                      Line2D([0], [0], marker='o', color='b', label='Z-projection',
                          markerfacecolor='b', markersize=5)]
    ax.legend(handles=legend_elements,scatterpoints = 1, loc='best')

    # Plot main directions and all the equivalent directions in cubic symmetry
    print_label = 0
    max_lim = 0
    for direction in familly_to_plot:
        equivalent_directions = apply_symmetry(direction,space_group = 221)
        pole_projections = stereographic_projection(direction, equiv = 1, sg= 221)
        if np.amax(pole_projections[:,0]) > max_lim : 
            max_lim = np.amax(pole_projections[:,0])
        for i, p in enumerate(pole_projections):
            ax.scatter(p[0],p[1],color="k",marker='x',s=10)
            if print_label :
                ax.annotate(np.array_str(equivalent_directions[i]), (p[0]+0.05,p[1]+0.05))
            
    flag_lines = False

    if flag_lines:
        # trace unit circle
        theta = np.linspace(-np.pi, np.pi, 200)
        ax.plot(np.sin(theta), np.cos(theta))
        # trace main x,y axis
        ax.plot((0,0),(-1,1),color="b")
        ax.plot((-1,1),(0,0),color="b")

    max_lim += 0.1
    plt.xlim([-max_lim,max_lim])
    plt.ylim([-max_lim,max_lim])

    ax.xaxis.set_major_formatter(NullFormatter())
    ax.yaxis.set_major_formatter(NullFormatter())

    ### Standard Stereographic Triangle (SST) plot ###

    # Projection along x-direction
    ax=plt.subplot(222)
    ax.set_title('Projection along X-direction')
    for at in range(om.shape[0]):
        spx_sst, min_pole = stereographic_projection(om[at,0,:], SST = 1, sg = 221)
        print(spx_sst)
        print(min_pole.shape)
        rgb = poles_to_colors(min_pole)
        ax.scatter(spx_sst[0],spx_sst[1],color=rgb,s=10,linewidths=.2,edgecolors='k')
    ax.xaxis.set_major_formatter( NullFormatter() )
    ax.yaxis.set_major_formatter( NullFormatter() )

    # IPF plot boundaries
    pole_boundaries = np.array([[0,0,1],
                              [0,1,1],
                              [-1,1,1],
                              [-1,3,3],
                              [-1,2,2]])
    pole_boundaries_txt = ['001',
                '011',
                '-111',
                '-133',
                '-122']

    pole_boundaries_projections = [stereographic_projection(i) for i in pole_boundaries]

    for i, p in enumerate(pole_boundaries_projections):
        ax.scatter(p[0],p[1],color="k",s=10)
        ax.annotate(pole_boundaries_txt[i], (p[0]+0.01,p[1]+0.01))
    ax.plot((pole_boundaries_projections[0][0],pole_boundaries_projections[1][0]),(pole_boundaries_projections[0][1],pole_boundaries_projections[1][1]),color='k')
    ax.plot((pole_boundaries_projections[0][0],pole_boundaries_projections[2][0]),(pole_boundaries_projections[0][1],pole_boundaries_projections[2][1]),color='k')
    
    # ----

    # Projection along y-direction
    ax = plt.subplot(223)
    ax.set_title('Projection along Y-direction')
    for at in range(om.shape[0]):
        spy_sst, min_pole = stereographic_projection(om[at,1,:], SST = 1, sg = 221)
        ax.scatter(spy_sst[0],spy_sst[1],color="g",s=10)
    ax.xaxis.set_major_formatter( NullFormatter() )
    ax.yaxis.set_major_formatter( NullFormatter() )


    for i, p in enumerate(pole_boundaries_projections):
        ax.scatter(p[0],p[1],color="k",s=10)
        ax.annotate(pole_boundaries_txt[i], (p[0]+0.01,p[1]+0.01))
    ax.plot((pole_boundaries_projections[0][0],pole_boundaries_projections[1][0]),(pole_boundaries_projections[0][1],pole_boundaries_projections[1][1]),color='k')
    ax.plot((pole_boundaries_projections[0][0],pole_boundaries_projections[2][0]),(pole_boundaries_projections[0][1],pole_boundaries_projections[2][1]),color='k')

    # ----

    # Projection along z-direction
    ax=plt.subplot(224)
    ax.set_title('Projection along Z-direction')
    for at in range(om.shape[0]):
        spz_sst, min_pole = stereographic_projection(om[at,2,:], SST = 1, sg = 221)
        ax.scatter(spz_sst[0],spz_sst[1],color="b",s=10)
    ax.xaxis.set_major_formatter( NullFormatter() )
    ax.yaxis.set_major_formatter( NullFormatter() )

    for i, p in enumerate(pole_boundaries_projections):
        ax.scatter(p[0],p[1],color="k",s=10)
        ax.annotate(pole_boundaries_txt[i], (p[0]+0.01,p[1]+0.01))
    ax.plot((pole_boundaries_projections[0][0],pole_boundaries_projections[1][0]),(pole_boundaries_projections[0][1],pole_boundaries_projections[1][1]),color='k')
    ax.plot((pole_boundaries_projections[0][0],pole_boundaries_projections[2][0]),(pole_boundaries_projections[0][1],pole_boundaries_projections[2][1]),color='k')

    plt.show()

# ---- Ovitos custom modifier

# unskew the cell, shift cell origine half-z up and wrap atoms back at PBC for clarity and get velocity = 0 at the median plane -- require creation of a new modifier
def unskew_xz(frame, input, output):
    global home, ID
    # Original simulation cell is passed through by default.
    # Output simulation cell is just a reference to the input cell.
    assert(output.cell is input.cell)

    # Make a copy of the simulation cell:
    cell = output.copy_if_needed(output.cell)

    # copy_if_needed() made a deep copy of the simulation cell object.
    # Now the the input and output each point to different objects.
    assert(cell is output.cell)
    assert(cell is not input.cell)

    # Now it's safe to modify the object copy:
    mat = cell.matrix.copy()
    
    # Using the marray to access a modifiable array containing the initial input information
    pos = input.particle_properties.position.marray
        
    pos[:,0] -= mat[0,2]/2 + flip[0]*mat[0,0]/2
    
    flip[1] = mat[0,2]
    mat[0,2] = 0.0
    cell.matrix = mat
    
def assign_color_ptm(frame, input, output):
    global proj
    orientation_property = input.particle_properties.orientation.array
    print(">>>Projecting {0} FCC atoms along {1}-directions for frame {2}".format(len(orientation_property),projection_axis, frame))
    om = quaternions_to_orientation_matrix(orientation_property)
    
    xy_full = stereographic_projection(om[:,proj,:])
    
    HSV = np.zeros((om[:,proj,:].shape[0],3))
    HSV[:,0] = np.arctan2(xy_full[:,1],xy_full[:,0])+np.pi
    HSV[:,1] = np.linalg.norm(xy_full,axis = 1)
    HSV[:,2] = 1
    RGB = HSV_to_RGB(HSV)

    Cprop = input.particles['Color']
    with Cprop:
        Cprop[...] = RGB

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------

############
### Main ###
######## 

''' Inverse (100) pole figure and unit standard stereographic triangle (SST) '''

usage = "usage: ovitos script.py fn.nc rendering_axis projection_axis frame flip"


# Settings 
production = 1      # [0,1] production run uses NetCDF traj file, non production file use a sample crystal from Atomsk
show_plot = 0

rendering_axis = "x"
projection_axis = "z" # [x,y,z]
projection_pole = np.array([[1,1,1]])

flag_main_axis = 0
flag_legend = 0

if production:
    try:
        fn = sys.argv[1] #  "../traj_0001-0002.nc"
        rendering_axis = str(sys.argv[2])
        projection_axis = str(sys.argv[3])
        frame_nb = int(sys.argv[4])
        flip = np.array([int(sys.argv[5]),1])
        slice = 200
    except:
        print(usage)
        quit()
else:
    flip = np.array([0,0])

if not production:
    theta_small = np.radians(-0.001)

    fn = "fcc_base.lmp"
    frame = 0

    try:
        os.system("rm {0}".format(fn)) # Removing Atomsk file
    except:
        pass
    
    # Various test cases from Atomsk configurations
    # os.system("$HOME/bin/atomsk --create fcc 4.08 Au orient [11-1] [2-11] [011]  -duplicate 1 1 1 -ig {0}".format(fn))
    os.system("$HOME/bin/atomsk --create fcc 4.08 Au orient [100] [010] [001] -duplicate 4 4 4 -ig {0}".format(fn))
    # Polycrystal test
    #os.system("$HOME/bin/atomsk --create fcc 4.046 Al -ig aluminium.xsf")
    #os.system("$HOME/bin/atomsk atomsk --polycrystal aluminium.xsf polycrystal.txt -ig {0} -ow".format(fn))

if projection_axis == "x":
    proj = 0
elif projection_axis == "y":
    proj = 1
elif projection_axis == "z":
    proj = 2
else:
    raise ValueError('Incorrect axis "{0}", please specify x, y or z.'.format(projection_axis))

# Obtaining a list of quaternions orientation for the system from the Polyhedral template matching for different orientation
# Import atomsk or trajectory file 'fn'
pipeline = import_file(fn,multiple_frames = True)
pipeline.add_to_scene()
cell_vis = pipeline.get_vis(SimulationCellVis)
cell_vis.render_cell = False
data = pipeline.compute(frame = frame_nb)
cell = data.expect(SimulationCell)
slice_dist = cell[1,3]

if not production:
    ''' If crystal too perfect, PTM needs introduction of small variation'''
    # Rotate cell around x-axis
    mod_x = AffineTransformationModifier(
            transform_particles = True,
            transform_box = True,
            transformation = [[1,0,0,0],
                            [0,np.cos(np.radians(theta_small)),-np.sin(np.radians(theta_small)),0],
                            [0,np.sin(np.radians(theta_small)),np.cos(np.radians(theta_small)),0]])
    node.modifiers.append(mod_x)

# Run PTM analysis
modptm = PolyhedralTemplateMatchingModifier()
modptm.output_orientation = True
pipeline.modifiers.append(modptm)

# Unskew the cell and set middle plan velocity to 0 for render
pipeline.modifiers.append(PythonScriptModifier(function = unskew_xz))

# Wrap PBC
modifier = WrapPeriodicImagesModifier()
pipeline.modifiers.append(modifier)


if production:
    # Select only a slice to reduce computational costs
    modslice = SliceModifier()
    if rendering_axis == 'x':
        modslice.normal = (1,0,0)
    if rendering_axis == 'y':
        modslice.normal = (0,1,0)
        modslice.distance = slice_dist
    if rendering_axis == 'z':
        modslice.normal = (0,0,1)
        modslice.distance = cell[2,2]/2.
    modslice.slice_width = slice
    pipeline.modifiers.append(modslice)

#select only FCC type, invert selection, and delete
SelectMod = SelectTypeModifier(property = "Structure Type")
SelectMod.types = { CommonNeighborAnalysisModifier.Type.FCC }
pipeline.modifiers.append(SelectMod)
InvMod = InvertSelectionModifier()
pipeline.modifiers.append(InvMod)
DelMod = DeleteSelectedModifier()
pipeline.modifiers.append(DelMod)

# Assign new RGB color depending on PTM projection
pipeline.modifiers.append(PythonScriptModifier(function = assign_color_ptm))

# Extract orientation from PTM
output = pipeline.compute(frame = frame_nb)
orientation_property = output.particle_properties.orientation.array
print(">>>Projecting {0} FCC atoms along {1}-directions for plots".format(len(orientation_property),projection_axis))
om = quaternions_to_orientation_matrix(orientation_property)

# Projection along axis "proj"
min_poles = poles_to_SST(om[:,proj,:], SG = 221)
RGB_min = poles_to_RGB(min_poles)
xy = stereographic_projection(min_poles)
xy_full = stereographic_projection(om[:,proj,:])
xy_max = stereographic_projection(np.array([[0,1,0]]))[:,0]

#---- Inverse Pole Figures ----#

#--- Plot 1 ---#
ax=plt.subplot(221)

# Plot limits and main pole directions
# main pole direction
familly_to_plot = np.array([[0,0,1],
                              [0,1,5],
                              [0,1,3],
                              [0,1,2],
                              [0,2,3],
                              [0,1,1],
                              [-1,3,3],
                              [-1,2,2],
                              [-1,-1,-1],
                              [-1,1,2],
                              [-1,1,3],
                              [-1,1,5],
                              [-1,2,3]])
main_poles = poles_to_SST(familly_to_plot, SG = 221)
xy_sst_poles = stereographic_projection(main_poles)
plt.scatter(xy[:,0],xy[:,1],s=30,c='k')

ax.scatter(xy[:,0],xy[:,1],s=10,c=RGB_min)
ax.scatter(xy_sst_poles[:,0],xy_sst_poles[:,1],s=30,c='k',marker='x')
ax.xaxis.set_major_formatter(NullFormatter())
ax.yaxis.set_major_formatter(NullFormatter())
ax.set_aspect(1)

# Plot 2
ax=plt.subplot(222)

HSV = np.zeros((om[:,proj,:].shape[0],3))
HSV[:,0] = np.arctan2(xy_full[:,1],xy_full[:,0])+np.pi
HSV[:,1] = np.linalg.norm(xy_full,axis = 1)
HSV[:,2] = 1
RGB = HSV_to_RGB(HSV)

ax.scatter(xy_full[:,0],xy_full[:,1],s=10,c=RGB)

ax.set_xlim([-np.sqrt(2)-.1,np.sqrt(2)+.1])
ax.set_ylim([-np.sqrt(2)-.1,np.sqrt(2)+.1])


# Plot limits and main pole directions
# trace unit circle
theta = np.linspace(-np.pi, np.pi, 200)
ax.plot(xy_max*np.sin(theta), xy_max*np.cos(theta),color='k')
# trace main x,y axis
ax.plot((0,0),(-xy_max,xy_max),color="k")
ax.plot((-xy_max,xy_max),(0,0),color="k")
ax.xaxis.set_major_formatter(NullFormatter())
ax.yaxis.set_major_formatter(NullFormatter())
ax.set_aspect(1)
# main pole direction
main_poles = equivalent_poles(familly_to_plot, SG = 221)
xy = stereographic_projection(main_poles)
plt.scatter(xy[:,0],xy[:,1],s=30,c='k',marker='x')

#--- Plot 3 ---#
ax=plt.subplot(224)

X,Y,Z = stereographic_projection(om[:,proj,:], bin_and_reduce = True)

if production:
    xlim = np.linspace(-np.sqrt(2),np.sqrt(2),50)
    ylim = (2.0-xlim**2)
    ylim[ylim < 1e-5] = 0.0  # avoid -1e-16 error with sqrt
    ylim = np.sqrt(ylim) 
    xlim = np.hstack((xlim,xlim[::-1]))
    ylim = np.hstack((ylim,-ylim[::-1]))
    merged_lim = np.vstack((xlim,ylim))
    merged = np.vstack((np.vstack((X.flatten(),Y.flatten())),Z.flatten()))
    np.savetxt('./pole_001_{0}_{1}_slice_{2}_histogram2d.txt'.format(projection_axis,frame_nb,slice),merged.T)
    np.savetxt('./pole_lim.txt',merged_lim.T)


# print(np.mean(Z))
cf = ax.contourf(X, Y, Z,cmap=plt.cm.jet, levels = np.linspace(0,5*np.mean(Z),200), extend="both")
ax.set_xlim([-np.sqrt(2)-.1,np.sqrt(2)+.1])
ax.set_ylim([-np.sqrt(2)-.1,np.sqrt(2)+.1])
ax.set_aspect(1)

# Plot limits and main pole directions
# trace unit circle
theta = np.linspace(-np.pi, np.pi, 200)
ax.plot(xy_max*np.sin(theta), xy_max*np.cos(theta),color='k')
# trace main x,y axis
ax.plot((0,0),(-xy_max,xy_max),color="k")
ax.plot((-xy_max,xy_max),(0,0),color="k")
ax.xaxis.set_major_formatter(NullFormatter())
ax.yaxis.set_major_formatter(NullFormatter())
ax.set_aspect(1)
# main pole direction
familly_to_plot = np.array([[0,0,1],
                              [0,1,5],
                              [0,1,3],
                              [0,1,2],
                              [0,2,3],
                              [0,1,1],
                              [-1,3,3],
                              [-1,2,2],
                              [-1,-1,-1],
                              [-1,1,2],
                              [-1,1,3],
                              [-1,1,5],
                              [-1,2,3]])
                              
main_poles = equivalent_poles(familly_to_plot, SG = 221)
xy = stereographic_projection(main_poles)
plt.scatter(xy[:,0],xy[:,1],s=30,c='k',marker='x')

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)

plt.colorbar(cf, cax=cax, ticks=[])

#--- Plot 4 ---#
ax=plt.subplot(223)

X,Y,Z = stereographic_projection(min_poles, bin_and_reduce = True, SST = True)

if production:
    SST_boundary = np.array([[0,0,1],
                              [1,1,1],
                              [5,5,4],
                              [4,4,3],
                              [7,7,5],
                              [3,3,2],
                              [5,5,3],
                              [2,2,1],
                              [7,7,3],
                              [5,5,2],
                              [3,3,1],
                              [4,4,1],
                              [5,4,3],
                              [5,5,1],
                              [6,6,1],
                              [7,7,1],
                              [1,1,0]])
    pole_boundary = poles_to_SST(SST_boundary, SG = 221)
    xy_boundary = stereographic_projection(pole_boundary)
    
    merged = np.vstack((np.vstack((X.flatten(),Y.flatten())),Z.flatten()))
    np.savetxt('./pole_SST_001_{0}_{1}_slice_{2}_histogram2d.txt'.format(projection_axis,frame_nb,slice),merged.T)
    np.savetxt('./pole_SST_lim.txt',xy_boundary)


cf = ax.contourf(X, Y, Z,cmap=plt.cm.jet, levels = np.linspace(0,5*np.mean(Z),200), extend="both") 
ax.set_aspect(1)

ax.scatter(xy_sst_poles[:,0],xy_sst_poles[:,1],s=30,c='k',marker='x')

ax.xaxis.set_major_formatter(NullFormatter())
ax.yaxis.set_major_formatter(NullFormatter())
ax.set_aspect(1)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)

plt.colorbar(cf, cax=cax, ticks=[])

#---- Pole Figures ----#

ax2=plt.subplot(111)

# Plot limits and main pole directions
# main pole direction
familly_to_plot = np.array([[0,0,1],
                              [0,1,5],
                              [0,1,3],
                              [0,1,2],
                              [0,2,3],
                              [0,1,1],
                              [-1,3,3],
                              [-1,2,2],
                              [-1,-1,-1],
                              [-1,1,2],
                              [-1,1,3],
                              [-1,1,5],
                              [-1,2,3]])
main_poles = poles_to_SST(familly_to_plot, SG = 221)
xy_sst_poles = stereographic_projection(main_poles)
plt.scatter(xy[:,0],xy[:,1],s=30,c='k')

ax2.scatter(xy[:,0],xy[:,1],s=10,c=RGB_min)
ax2.scatter(xy_sst_poles[:,0],xy_sst_poles[:,1],s=30,c='k',marker='x')
ax2.xaxis.set_major_formatter(NullFormatter())
ax2.yaxis.set_major_formatter(NullFormatter())
ax2.set_aspect(1)

if production:
    # Rendering
    # pipeline.add_to_scene()
    # cell_vis = pipeline.get_vis(SimulationCellVis)
    # cell_vis.render_cell = False

    vp = Viewport(type = Viewport.Type.Ortho)

    if rendering_axis == 'x':
        vp.camera_dir = (1,0,0)
    if rendering_axis == 'y':
        vp.camera_dir = (0,1,0)
    if rendering_axis == 'z':
        vp.camera_dir = (0,0,1)
    vp.zoom_all()
        
    vp.render_image(filename='./projection_001_view-{0}_proj-{1}_frame-{2}_slice-{3}.png'.format(rendering_axis,projection_axis,frame_nb,slice),
                frame = frame_nb,
                size=(1600,1200),
                background=(0,0,0), 
                renderer=TachyonRenderer(ambient_occlusion=False, shadows=False))
    # vp.render_anim(filename='./test_.png',
                # range = (0,5),
                # size=(1600,1200),
                # background=(0,0,0), 
                # renderer=TachyonRenderer(ambient_occlusion=False, shadows=False))

    plt.savefig('./pole_SST_001_view-{0}_proj-{1}_frame-{2}_slice-{3}.png'.format(rendering_axis,projection_axis,frame_nb,slice), dpi=1200)
    pipeline.remove_from_scene()  # Required to avoid a final pipeline evaluation at the end of the script
    
if show_plot:
    plt.show()

# Misc. #
if flag_main_axis:
    ''' Plot main axis of the SST '''
    familly_to_plot = np.array([[0,0,1],
                                  [0,1,5],
                                  [0,1,3],
                                  [0,1,2],
                                  [0,2,3],
                                  [0,1,1],
                                  [-1,3,3],
                                  [-1,2,2],
                                  [-1,-1,-1],
                                  [-1,1,2],
                                  [-1,1,3],
                                  [-1,1,5],
                                  [-1,2,3]])
                                  
    min_poles = poles_to_SST(familly_to_plot, SG = 221)
    RGB = poles_to_RGB(min_poles)
    xy = stereographic_projection(min_poles)

    plt.scatter(xy[:,0],xy[:,1],s=100,c=RGB)

    plt.show()

if flag_legend:
    ''' Plot SST colormap '''
    
    familly_to_plot = np.array([[0,0,1],
                                  #[0,1,5],
                                  #[0,1,3],
                                  #[0,1,2],
                                  #[0,2,3],
                                  [0,1,1],
                                  #[-1,3,3],
                                  #[-1,2,2],
                                  [-1,-1,-1],
                                  [-1,1,2],
                                  #[-1,1,3],
                                  #[-1,1,5],
                                  [-1,2,3]])
    
    # Create a (Nx3) array of random pole direction to fill up the SST and get entire color map
    N = 75

    color_bar_array = np.zeros((N**3,3))

    i = 0
    for a in np.linspace(0,1,N):
        for b in np.linspace(0,1,N):
            for c in np.linspace(0,1,N):
                progressbar(i,N**3) 
                if c == 0:
                    color_bar_array[i] = np.array([a,b,1.])
                else:
                    color_bar_array[i] = np.array([a,b,c])
                i+=1

                
    min_poles = poles_to_SST(color_bar_array, SG = 221)
    RGB = poles_to_RGB(min_poles)
    xy = stereographic_projection(min_poles)

    plt.scatter(xy[:,0],xy[:,1],s=60,c=RGB,linewidths=0)
    
        # main pole direction
    main_poles = poles_to_SST(familly_to_plot, SG = 221)
    xy_sst_poles = stereographic_projection(main_poles)
    plt.scatter(xy_sst_poles[:,0],xy_sst_poles[:,1],s=30,c='k',marker='x')
    
    plt.axes().set_aspect('equal', 'datalim')
    
    plt.savefig('./SST_color_map_N_{0}.png'.format(N), dpi=1200)
    plt.show()
    
    # Create a (Nx3) array of random pole direction to fill up the (100) inverse pole figure 
    
    i = 0
    for a in np.linspace(-1,1,N):
        for b in np.linspace(-1,1,N):
            for c in np.linspace(-1,1,N):
                progressbar(i,N**3) 
                if c == 0:
                    color_bar_array[i] = np.array([a,b,1.])
                else:
                    color_bar_array[i] = np.array([a,b,c])
                i+=1
    
    xy = stereographic_projection(color_bar_array)
    HSV = np.zeros((color_bar_array.shape[0],3))
    HSV[:,0] = np.arctan2(xy[:,1],xy[:,0])+np.pi
    HSV[:,1] = np.linalg.norm(xy,axis = 1)
    HSV[:,2] = 1
    RGB = HSV_to_RGB(HSV)
    plt.scatter(xy[:,0],xy[:,1],s=60,c=RGB,linewidths=0)
    
        # main poles direction
    main_poles = equivalent_poles(familly_to_plot, SG = 195)
    xy = stereographic_projection(main_poles)
    np.savetxt('./main_poles_xy_coordinates.txt',xy)
    plt.scatter(xy[:,0],xy[:,1],s=30,c='k',marker='x')
    
    plt.axes().set_aspect('equal', 'datalim')
    
    plt.savefig('./stereo_color_map_N_{0}.png'.format(N), dpi=1200)
    plt.show()


plot = 0
if plot:
    plot_IPF(om)
