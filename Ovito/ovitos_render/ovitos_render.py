#!/usr/bin/python3
""" Usage : ovitos ovitos_render.py <parameter file>

    This script creates graphical render from a NetCDF trajectory file
    
    Parameters:
        infile
            Path to input NetCDF trajectory file, str
        outfile
            Path to the output png file(s), str
        modfile -- optional
            Path to the file containing extra ovitos modifiers to apply before rendering, str
        width 
            width in pixel of the rendered png, int
        height
            height in pixel of the rendered png, int
        type
            type of projection used in the viewport [ORTHO/PERSPECTIVE], str
        pos
            Position of the camera in the viewport, 3 floats -- separator = spaces
        dir
            Direction where the camera is pointing, 3 floats -- separator = spaces
        fov
            Field of view of the camera, float
        range
            range of the animation [full/single_frame/custom], str
        frame -- if single_frame
            frame number to export, int
        start -- if custom 
            starting frame for the custom range, int
        end -- if custom
            ending frame for the custom range, int
                
    Here is a minimal example for a parameter file:
    
        [files]
        infile = traj_indent.nc
        outfile = ./rendered_png/dxa_only_.png
        modfile = ./ovitos_dxa.mod

        [settings]
        width = 1600
        height = 1200
        type = ORTHO
        pos = 167.17 159.382 299.442
        dir = -0.867618 0.49721 0.00468931
        fov = 299.741
        range = single_frame
        frame = 100
    
"""

# Import OVITO modules.
from ovito.io import *
from ovito.modifiers import *
from ovito.modifiers import PythonScriptModifier
from ovito.vis import *

# Import NumPy module.
import numpy as np

import configparser
import sys
import os.path

################################################################################################

def main():
    # Parse configuration file
    configfile = sys.argv[1]
    if not os.path.isfile(configfile):
        raise ValueError('invalid path to config file %s', configfile)
    config = configparser.RawConfigParser()
    config.read(configfile)

    # Get file names
    infile = config.get('files','infile')
    outfile = config.get('files','outfile') 
    
    # Read animation type
    rr = config.get('settings','range')
    if rr == "full":
        rs = RenderSettings.Range.ANIMATION
    elif rr == "custom":
        rs = RenderSettings.Range.CUSTOM_INTERVAL
    elif rr == "single_frame":
        rf = config.get('settings','frame')
        rs = RenderSettings.Range.CURRENT_FRAME
    else:
        raise ValueError('invalid range setting "{0}" in config file', rr)

    # Load file
    node = import_file(infile, multiple_frames = True)
    node.add_to_scene()
    cell = node.source.cell
    cell.display.enabled = False         # This hides the simulation cell 
    
    
    
    # Insert progress bar
    def progress(frame, input, output):
        n_frames = int(node.source.num_frames)
        sys.stdout.write('{0}/{1}\r'.format(frame, n_frames, len=len))
        sys.stdout.flush() 
    node.modifiers.append(PythonScriptModifier(function = progress))
    
    
    
    
    if "modfile" in config['files']:
        exec(open(config.get('files','modfile')).read())

    if rr == "single_frame":
        ovito.dataset.anim.current_frame = int(rf)    
        node.compute()                                      # Ensure pipeline output exists
    elif rr == "custom":
        ovito.dataset.anim.current_frame = int(config.get('settings','start'))
        node.compute()
        settings = RenderSettings(
            range = rs,
            custom_range = (int(config.get('settings','start')),int(config.get('settings','end'))),
            filename = outfile,
            size = (int(config.get('settings','width')), int(config.get('settings','height')))
        )
    else:
        settings = RenderSettings(
            range = rs,
            filename = outfile,
            size = (int(config.get('settings','width')), int(config.get('settings','height')))
        )
    
    settings.renderer = TachyonRenderer() # Activate the TachyonRenderer backend
    settings.renderer.shadows = False     # Turn off cast shadows

    vp = Viewport()
    if config.get('settings','type') == "ORTHO":
        vp.type = Viewport.Type.ORTHO 
    elif config.get('settings','type') == "PERSPECTIVE":
         vp.type = Viewport.Type.PERSPECTIVE
    else:
        raise ValueError('invalid viewport type "{0}" in config file', config.get('settings','type'))
    vp.camera_pos = tuple([float(x) for x in config.get('settings','pos').split()])
    vp.camera_dir = tuple([float(x) for x in config.get('settings','dir').split()])
    vp.fov = float(config.get('settings','fov'))

    ##vp.zoom_all()

    vp.render(settings)

if __name__=='__main__':
    main()
