# Ovito Script 3.x syntax


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------

############
### Main ###
############

# Import OVITO modules.
from ovito.io import *
#from ovito.pipeline import StaticSource, Pipeline
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
from matplotlib.ticker import NullFormatter, NullLocator
from matplotlib.colors import LogNorm
from matplotlib.font_manager import FontProperties

class PTM_SST():
    def __init__(self, fn, rendering_axis, projection_axis, frame_nb, flip, slice, plot_data, show_plot = 0):
        self.fn = fn 
        self.rendering_axis = rendering_axis 
        self.projection_axis = projection_axis 
        self.frame = frame_nb 
        self.flip = flip 
        self.slice = slice
        self.show_plot = show_plot
        self.plot_data = plot_data
        
        if __name__ == "__main__":
            # Printing parameters
            print("trajectory file:", self.fn)
            print("rendering axis:", self.rendering_axis)
            print("rendering slice:", self.slice)
            print("projection axis:", self.projection_axis)
            print("processing frame:", self.frame)
            print("accounting for shear flip:",self.flip)
            print("plotting data:",self.plot_data)
            print("showing_plot:",self.show_plot)

        #### Important poles directions ####
        # Discrete poles forming the boundary of the standard stereographic triangle for cubic system
        self.SST_boundaries = np.array([[0,0,1],
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

        self.projection_pole = np.array([[1,1,1]])

        if self.projection_axis == "x":
            self.proj = 0
        elif self.projection_axis == "y":
            self.proj = 1
        elif self.projection_axis == "z":
            self.proj = 2
        else:
            raise ValueError('Incorrect axis "{0}", please specify x, y or z.'.format(self.projection_axis))
            
        # Obtaining a list of quaternions orientation for the system from the Polyhedral template matching for different orientation
        # Import trajectory file 'fn'
        self.pipeline = import_file(self.fn,multiple_frames = True)
        self.pipeline.add_to_scene()
        cell_vis = self.pipeline.get_vis(SimulationCellVis)
        cell_vis.render_cell = False
        data = self.pipeline.compute(frame = self.frame)
        cell = data.expect(SimulationCell)
        
        if __name__ == "__main__":
            print('cell')
            print(cell)
            print('distances')
            print('cell_1_1',cell[1,1])
            print('cell_1_3',cell[1,3])
            print('slice distance', cell[1,1]+cell[1,3])
            print('slice',self.slice)

        # Run PTM analysis
        modptm = PolyhedralTemplateMatchingModifier()
        modptm.output_orientation = True
        self.pipeline.modifiers.append(modptm)

        # Unskew the cell and set middle plan velocity to 0 for render
        self.pipeline.modifiers.append(PythonScriptModifier(function = self.unskew_xz))

        # Wrap PBC
        modifier = WrapPeriodicImagesModifier()
        self.pipeline.modifiers.append(modifier)


        # Select only a slice to reduce computational costs
        modslice = SliceModifier()
        if self.rendering_axis == 'x':
            modslice.normal = (1,0,0)
            modslice.distance = cell[0,0]
        if self.rendering_axis == 'y':
            modslice.normal = (0,1,0)
            # modslice.distance = cell[1,1]+cell[1,3]
            modslice.distance = cell[1,3]   # cell origin
        if self.rendering_axis == 'z':
            modslice.normal = (0,0,1)
            modslice.distance = cell[2,2]/2.
        modslice.slice_width = self.slice
        self.pipeline.modifiers.append(modslice)

        #select only FCC type, invert selection, and delete
        SelectMod = SelectTypeModifier(property = "Structure Type")
        SelectMod.types = { CommonNeighborAnalysisModifier.Type.FCC }
        self.pipeline.modifiers.append(SelectMod)
        InvMod = InvertSelectionModifier()
        self.pipeline.modifiers.append(InvMod)
        DelMod = DeleteSelectedModifier()
        self.pipeline.modifiers.append(DelMod)

        # Assign new RGB color depending on PTM projection, 
        # SELECT option -1- or -2-
        
        # - 1 - Assign color to atoms from the HSV mapping of the stereographic projection area (unwrapped symmetry)
        # print("Coloring according to full-PTM stereographic projection, i.e. 4 quadrants")
        # pipeline.modifiers.append(PythonScriptModifier(function = assign_color_ptm_full))
        # - 2 - Assign color to atoms from RGB mapping of the standard stereographic triangle area (wrap all symmetry)
        print("Coloring according to SST mapping")
        self.pipeline.modifiers.append(PythonScriptModifier(function = self.assign_color_ptm_SST))


        if self.plot_data == 1 or self.plot_data == True:
            # Extract orientation from PTM
            output = self.pipeline.compute(frame = self.frame)
            orientation_property = output.particle_properties.orientation.array
            print(">>>Projecting {0} FCC atoms along {1}-directions for plots".format(len(orientation_property),self.projection_axis))
            om = self.quaternions_to_orientation_matrix(orientation_property)
            self.plot_IPF(om)

        # Rendering
        vp = Viewport(type = Viewport.Type.Ortho)

        if rendering_axis == 'x':
            vp.camera_dir = (1,0,0)
        if rendering_axis == 'y':
            vp.camera_dir = (0,1,0)
        if rendering_axis == 'z':
            vp.camera_dir = (0,0,1)
        # if __name__ != "__main__":
        vp.zoom_all()
            
        vp.render_image(filename='./projection_001_view_{0}_proj_{1}_slice_{2}_frame_{3:04d}.png'.format(rendering_axis,self.projection_axis,slice,self.frame),
                    frame = self.frame,
                    size=(1600,1200),
                    background=(0,0,0), 
                    renderer=TachyonRenderer(ambient_occlusion=False, shadows=False))
                    
        self.pipeline.remove_from_scene()  # Required to avoid a final pipeline evaluation at the end of the script
        del self.pipeline
        
        # if __name__ != "__main__":
            # self.get_histogram_data(self.frame)
        
        
            
    def get_histogram_data(self, frame):
        # Get histogram data
        output = self.pipeline.compute(frame = frame)
        orientation_property = output.particle_properties.orientation.array
        print(">>>Projecting {0} FCC atoms along {1}-directions for histogram data".format(len(orientation_property),self.projection_axis))
        om = self.quaternions_to_orientation_matrix(orientation_property)
        min_poles = self.poles_to_SST_fast(om[:,self.proj,:], SG = 221)
        X,Y,Z = self.stereographic_projection(min_poles, bin_and_reduce = True, SST = True)        
        merged = np.vstack((np.vstack((X.flatten(),Y.flatten())),Z.flatten()))
        np.savetxt('./pole_SST_001_{0}_{1}_slice_{2}_histogram2d.txt'.format(self.projection_axis,self.frame,self.slice),merged.T)
        
    def assign_color_ptm_full(self, frame, input, output):
        orientation_property = input.particle_properties.orientation.array
        print(">>>Projecting {0} FCC atoms along {1}-directions for frame {2}".format(len(orientation_property),self.projection_axis, self.frame))
        om = self.quaternions_to_orientation_matrix(orientation_property)
        
        xy_full = self.stereographic_projection(om[:,self.proj,:])
        
        HSV = np.zeros((om[:,proj,:].shape[0],3))
        HSV[:,0] = np.arctan2(xy_full[:,1],xy_full[:,0])+np.pi
        HSV[:,1] = np.linalg.norm(xy_full,axis = 1)
        HSV[:,2] = 1
        RGB = self.HSV_to_RGB(HSV)

        Cprop = input.particles['Color']
        with Cprop:
            Cprop[...] = RGB
        
    def assign_color_ptm_SST(self, frame, input, output):
        orientation_property = input.particle_properties.orientation.array
        print(">>>Projecting {0} FCC atoms along {1}-directions for frame {2}".format(len(orientation_property),self.projection_axis, self.frame))
        om = self.quaternions_to_orientation_matrix(orientation_property) 
        
        min_poles = self.poles_to_SST_sorted(om[:,self.proj,:], SG = 221)
        RGB_min = self.poles_to_RGB(min_poles)
        
        Cprop = input.particles['Color']
        with Cprop:
            Cprop[...] = RGB_min
            
        if __name__ != "__main__":
            print(">>>Projecting {0} FCC atoms along {1}-directions for histogram data".format(len(orientation_property),self.projection_axis))
            X,Y,Z = self.stereographic_projection(min_poles, bin_and_reduce = True, SST = True)        
            merged = np.vstack((np.vstack((X.flatten(),Y.flatten())),Z.flatten()))
            np.savetxt('./pole_SST_001_{0}_{1}_slice_{2}_histogram2d.txt'.format(self.projection_axis,self.frame,self.slice),merged.T)
            del merged, X,Y,Z, min_poles, om, orientation_property, RGB_min
        
    def unskew_xz(self, frame, input, output):
        ''' Ovitos custom modifier to unskew cell and shift cell origin half-z up under simple shear deformation'''
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
        
        flip_is_number = self.is_number(self.flip)
        if flip_is_number:
            pos[:,0] -= mat[0,2]/2 + self.flip*mat[0,0]/2
        else:
            flip_file = open(self.flip,"r").readlines()
            len_file = len(flip_file)
            # smaller than first flip registered
            if frame < int(flip_file[0].split(' ')[0]):
                F = 0
            # larger than last flip registered
            elif frame >= int(flip_file[len_file-1].split(' ')[0]):
                F = int(flip_file[len_file-1].split(' ')[1])
            else:
                for index,line in enumerate(flip_file):
                    if frame == int(line.split(' ')[0]):
                        F = int(flip_file[index].split(' ')[1])
                        break
                    elif frame < int(line.split(' ')[0]):
                        F = int(flip_file[index-1].split(' ')[1])
                        break
            pos[:,0] -= mat[0,2]/2 + F*mat[0,0]/2

        mat[0,2] = 0.0
        cell.matrix = mat
        
    def is_number(self, s):
        try:
            float(s)
            return True
        except ValueError:
            pass
        try:
            import unicodedata
            unicodedata.numeric(s)
            return True
        except (TypeError, ValueError):
            pass
        return False

    def quaternions_to_orientation_matrix(self, qs):
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
            
        return om
    
    def apply_symmetry(self, position,space_group):
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

        
    def stereographic_projection(self, poles, type = "equal_area", bin_and_reduce = False, SST = False):
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

            
            
        if bin_and_reduce: # Return pole density to plot contour map 
            N = 200
                   
            if SST:
                x_edges, d_x = np.linspace(0, self.stereographic_projection(self.poles_to_SST_fast(np.array([[0,1,1]]),SG = 221))[0][0], N, retstep=True)
                y_edges, d_y = np.linspace(0, self.stereographic_projection(self.poles_to_SST_fast(np.array([[1,1,1]]),SG = 221))[0][1], N, retstep=True)
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
    
    
    def equivalent_poles(self, poles, SG):
        """ Takes a list of poles (Nx3 array) and returns a list of ALL equivalent directions
        (Nx3 array) depending on the given symmetry of the space group SG
        """
        equivalent_directions = self.apply_symmetry(poles,SG)
        unique_poles =  np.vstack({tuple(row) for row in equivalent_directions}) 

        return unique_poles
        
    def poles_to_SST_fast(self, poles, SG):
        """ Takes a list of poles (Nx3 array) and returns the corresponding
        list of minimal poles (Nx3 array) lying in the standard stereographic triangle
        """
        
        # Return a shuffled list of poles, i.e. issue to assign correct color to correct atom ID

        if len(poles.shape) != 2:
            raise RuntimeError("poles must be a 2-dimensional array")
     
        if poles.shape[1] != 3:
            raise RuntimeError("poles must be a n x 3 dimensional array")

        # Get all equivalent directions from the symmetry of the space group 'SG'
        equivalent_directions = self.apply_symmetry(poles,SG)
                
        # Filter poles belonging to the standard stereographic triangle (SST), minimal symmetry zone for cubic symmetry.
        mask_SST = np.logical_and(equivalent_directions[:,2]-equivalent_directions[:,1]>=0,np.logical_and(equivalent_directions[:,1]+equivalent_directions[:,0]>=0,-equivalent_directions[:,0]>=0))

        min_poles =  np.vstack({tuple(row) for row in equivalent_directions[mask_SST]}) 
       
        return min_poles

    def poles_to_SST_sorted(self, poles, SG):
        """ Takes a list of poles (Nx3 array) and returns the corresponding
        list of minimal poles (Nx3 array) lying in the standard stereographic triangle
        """

        if len(poles.shape) != 2:
            raise RuntimeError("poles must be a 2-dimensional array")
     
        if poles.shape[1] != 3:
            raise RuntimeError("poles must be a n x 3 dimensional array")

        # Assigning each pole one by one to avoid any shuffling in the ID
        for i, p in enumerate(poles):
            p = p.reshape((1,3))
            equivalent_directions = self.apply_symmetry(p,SG)
            mask_SST = np.logical_and(equivalent_directions[:,2]-equivalent_directions[:,1]>=0,np.logical_and(equivalent_directions[:,1]+equivalent_directions[:,0]>=0,-equivalent_directions[:,0]>=0))
            min_p = equivalent_directions[mask_SST]
            poles[i] = min_p
       
        return poles

    def poles_to_RGB(self, poles):
        """ Takes a list of poles (Nx3 array) and returns a list of
            corresponding RGB colors (Nx3 array) """
            
        if len(poles.shape) != 2:
            raise RuntimeError("poles must be a 2-dimensional array")
     
        if poles.shape[1] != 3:
            raise RuntimeError("poles must be a n x 3 dimensional array")
        
        rgb = np.zeros_like(poles)
            
        h,k,l = poles.T 
        
        rgb[:,0] = l-k
        rgb[:,1] = k+h
        rgb[:,2] = -h
                
        # Normalize values.
        max = np.amax(rgb, axis = 1)
        rgb[:,0] = rgb[:,0] / max
        rgb[:,1] = rgb[:,1] / max
        rgb[:,2] = rgb[:,2] / max
         
        return rgb

    
    def HSV_to_RGB(self, HSV):
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
        
        
    def plot_IPF(self,om):
        # Projection along axis "proj"
        min_poles = self.poles_to_SST_fast(om[:,self.proj,:], SG = 221)
        RGB_min = self.poles_to_RGB(min_poles)
        xy = self.stereographic_projection(min_poles)
        xy_full = self.stereographic_projection(om[:,self.proj,:])
        xy_max = self.stereographic_projection(np.array([[0,1,0]]))[:,0]

        #---- Inverse Pole Figures ----#
        #--- Plot 1 - SST ---#
        ax=plt.subplot(221)
        # Plot limits and main pole directions
        # Main poles directions
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
        main_poles = self.poles_to_SST_fast(familly_to_plot, SG = 221)
        xy_sst_poles = self.stereographic_projection(main_poles)

        ax.scatter(xy[:,0],xy[:,1],s=10,c=RGB_min)
        ax.scatter(xy_sst_poles[:,0],xy_sst_poles[:,1],s=30,c='k',marker='x')
        ax.xaxis.set_major_formatter(NullFormatter())
        ax.yaxis.set_major_formatter(NullFormatter())
        ax.xaxis.set_major_locator(NullLocator())
        ax.yaxis.set_major_locator(NullLocator())
        ax.set_aspect(1)

        #--- Plot 2 - Stereo projection ---#
        ax=plt.subplot(222)

        HSV = np.zeros((om[:,self.proj,:].shape[0],3))
        HSV[:,0] = np.arctan2(xy_full[:,1],xy_full[:,0])+np.pi
        HSV[:,1] = np.linalg.norm(xy_full,axis = 1)
        HSV[:,2] = 1
        RGB = self.HSV_to_RGB(HSV)

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
        ax.xaxis.set_major_locator(NullLocator())
        ax.yaxis.set_major_locator(NullLocator())
        ax.set_aspect(1)
        # main pole direction
        main_poles = self.equivalent_poles(familly_to_plot, SG = 221)
        xy = self.stereographic_projection(main_poles)
        plt.scatter(xy[:,0],xy[:,1],s=20,c='k',marker='x')

        #--- Plot 3 - SST density ---#
        ax=plt.subplot(223)

        X,Y,Z = self.stereographic_projection(min_poles, bin_and_reduce = True, SST = True)
        max_level = len(str(int(np.amax(Z))))-1
        
        levls = np.logspace(0,max_level,400)*5
        font = FontProperties()
        font.set_size("xx-large")
        font.set_weight("bold")
        
        ax.patch.set_facecolor("#050080") # fill background color to avoid white spaces for values == 0
        cf = ax.contourf(X, Y, Z, norm = LogNorm(),levels = levls) 
                
        ax.text(.04,.01,"[001]",color="w",fontproperties=font)
        ax.text(.66,.01,"[011]",color="w",fontproperties=font)
        ax.text(.52,.62,"[111]",color="k",fontproperties=font)

        ax.xaxis.set_major_formatter(NullFormatter())
        ax.yaxis.set_major_formatter(NullFormatter())
        ax.xaxis.set_major_locator(NullLocator())
        ax.yaxis.set_major_locator(NullLocator())
        ax.set_aspect(1)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)

        plt.colorbar(cf, cax=cax, ticks=[])
            

        pole_boundary = self.poles_to_SST_fast(self.SST_boundaries, SG = 221)
        xy_boundary = self.stereographic_projection(pole_boundary)
        
        merged = np.vstack((np.vstack((X.flatten(),Y.flatten())),Z.flatten()))
        np.savetxt('./pole_SST_001_{0}_{1}_slice_{2}_histogram2d.txt'.format(self.projection_axis,self.frame,self.slice),merged.T)
        np.savetxt('./pole_SST_lim.txt',xy_boundary)

        #--- Plot 4 - stereographic projection density ---#
        ax=plt.subplot(224)

        X,Y,Z = self.stereographic_projection(om[:,self.proj,:], bin_and_reduce = True)

        xlim = np.linspace(-np.sqrt(2),np.sqrt(2),50)
        ylim = (2.0-xlim**2)
        ylim[ylim < 1e-5] = 0.0  # avoid -1e-16 error with sqrt
        ylim = np.sqrt(ylim) 
        xlim = np.hstack((xlim,xlim[::-1]))
        ylim = np.hstack((ylim,-ylim[::-1]))
        merged_lim = np.vstack((xlim,ylim))
        merged = np.vstack((np.vstack((X.flatten(),Y.flatten())),Z.flatten()))
        np.savetxt('./pole_001_{0}_{1}_slice_{2}_histogram2d.txt'.format(self.projection_axis,self.frame,self.slice),merged.T)
        np.savetxt('./pole_lim.txt',merged_lim.T)

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
        ax.xaxis.set_major_locator(NullLocator())
        ax.yaxis.set_major_locator(NullLocator())
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
                                      
        main_poles = self.equivalent_poles(familly_to_plot, SG = 221)
        xy = self.stereographic_projection(main_poles)
        plt.scatter(xy[:,0],xy[:,1],s=20,c='k',marker='x')

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)

        plt.colorbar(cf, cax=cax, ticks=[])

        plt.tight_layout()
        if self.show_plot:
            plt.show()
        plt.savefig('./pole_SST_001_view-{0}_proj-{1}_frame-{2}_slice-{3}.png'.format(rendering_axis,projection_axis,frame_nb,slice), dpi=1200)
        plt.clf()


if __name__ == "__main__":
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
    from matplotlib.ticker import NullFormatter, NullLocator

    import time
    
    start_time = time.time()
    ''' Inverse (100) pole figure and unit standard stereographic triangle (SST) '''

    usage = "usage: ovitos script.py fn.nc rendering_axis projection_axis frame flip slice"
    
    try:
        fn = sys.argv[1]
        rendering_axis = str(sys.argv[2])
        projection_axis = str(sys.argv[3])
        frame_nb = int(sys.argv[4])
        flip = int(sys.argv[5])
        slice = int(sys.argv[6])
    except:
        print(usage)
        quit()
        
    show_plot = 0
    plot_data = 1
    
    # Printing parameters
    # print("trajectory file:", fn)
    # print("rendering axis:", rendering_axis)
    # print("rendering slice:", slice)
    # print("projection axis:", projection_axis)
    # print("processing frame:", frame_nb)
    # print("accounting for shear flip:",flip)
    # print("plotting data:",plot_data)
    # print("showing_plot:",show_plot)
    
    PTM_SST(fn, rendering_axis, projection_axis, frame_nb, flip, slice,plot_data, show_plot)
    
    print("--- ran for %s seconds ---" % (time.time() - start_time))
