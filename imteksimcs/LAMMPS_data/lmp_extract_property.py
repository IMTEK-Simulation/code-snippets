#!/usr/bin/env ovitos
""" Extracts properties from LAMMPS data files and NetCDF trajectories """
# 8 Apr 2019, extracts properties from LAMMPS data files and NetCDF
# Johannes Hörmann, johannes.hoermann@imtek.uni-freiburg.de
# tested with Ovito 3.0.0-dev234 

import sys, logging
import numpy as np

from ovito.io import import_file, export_file
from ovito.modifiers import CommonNeighborAnalysisModifier
from ovito.modifiers import SelectParticleTypeModifier, InvertSelectionModifier
from ovito.modifiers import DeleteSelectedParticlesModifier
from ovito.modifiers import LoadTrajectoryModifier
from ovito.modifiers import CoordinationNumberModifier #CoordinationAnalysisModifier

logger = logging.getLogger(__name__)

# adapted from
# https://ovito.org/manual/python/modules/ovito_modifiers.html#ovito.modifiers.LoadTrajectoryModifier


def select_fcc_only(pipeline):
    """Deletes all particles not in FCC environment"""
    # Ovito 3.0.0-dev362 standard:
    # from ovito.modifiers import CommonNeighborAnalysisModifier
    # from ovito.modifiers import SelectTypeModifier, InvertSelectionModifier
    # from ovito.modifiers import DeleteSelectedModifier
  
    # neigh_mod = CommonNeighborAnalysisModifier(
    #   mode = CommonNeighborAnalysisModifier.Mode.AdaptiveCutoff )
    # sel_mod = SelectTypeModifier( 
    #  operate_on = "particles",
    #  property   = "Structure Type",
    #  types     =  { CommonNeighborAnalysisModifier.Type.FCC } )
    # inv_mod = InvertSelectionModifier( operate_on = "particles" )
    # del_mod = DeleteSelectedModifier( operate_on = "particles" )
  
    # Ovito 2.9.0 standard:
    # from ovito.modifiers import CommonNeighborAnalysisModifier
    # from ovito.modifiers import SelectParticleTypeModifier, InvertSelectionModifier
    # from ovito.modifiers import DeleteSelectedParticlesModifier
  
    neigh_mod = CommonNeighborAnalysisModifier(
      mode = CommonNeighborAnalysisModifier.Mode.AdaptiveCutoff )
    sel_mod = SelectParticleTypeModifier( 
      property   = "Structure Type",
      types     =  { CommonNeighborAnalysisModifier.Type.FCC } )
    inv_mod = InvertSelectionModifier()
    del_mod = DeleteSelectedParticlesModifier()
  
    pipeline.modifiers.append(neigh_mod)
    pipeline.modifiers.append(sel_mod)
    pipeline.modifiers.append(inv_mod)
    pipeline.modifiers.append(del_mod)
  
    return pipeline
   
def extract_nothing(pipeline, frames=[0]):
    """Dummy for testing purposes"""
    logger.info("Nothing extracted.")
    return None

def extract_rdf(pipeline, frames=[0], cutoff=8.0, number_of_bins=20000):
    """Returns M x 1 array of RDF averaged over N frames, where N = len(frames)
       and M = number_of_bins"""
    # Insert the modifier into the pipeline:
    # rdf_mod = CoordinationAnalysisModifier(cutoff = 8.0, number_of_bins = 20000)
    rdf_mod = CoordinationNumberModifier(cutoff=8.0, number_of_bins=20000)
    pipeline.modifiers.append(rdf_mod)

    # Initialize array for accumulated RDF histogram to zero:
    total_rdf    = np.zeros((rdf_mod.number_of_bins, 2))      
    # Iterate over all frames of the sequence.
    for n, frame in enumerate(frames):
        # Evaluate pipeline to let the modifier compute the RDF of the current frame:
        data = pipeline.compute(frame)
        # Accumulate RDF histograms:
        # total_rdf += data.series['coordination-rdf'].as_table()
        total_rdf += rdf_mod.rdf

    # Averaging:
    total_rdf /= len(frames)

    return total_rdf


def extract_box_measures(pipeline, frames=[0]):
    """Returns N x 3 array with box dimensions, where N = len(frames)"""
    cell_measures = np.zeros((len(frames),3))
    for n, frame in enumerate(frames):
        data = pipeline.compute(frame)
        cell = data.cell
        cell_measures[n,:] = [ cell[i,i] for i in range(3) ]

    average_cell_measures = np.mean(cell_measures, axis=0)
    logger.info("Initial and final cell measures: {}, {} ".format(cell_measures[0,:],cell_measures[-1,:]))
    logger.info("Average cell measures over {:d} frames: '{}'".format(len(frames),average_cell_measures))
    return cell_measures

def extract_property(
    topology_file   = 'datafile.lammps', # lammps data file
    trajectory_file = None, # netcdf trajectory
    frames          = None, # list of frame numbers
    property_func   = extract_box_measures,
    modifier_func   = None,
    outfile         = 'outfile.txt',
    atom_style      = 'full'):
    """Extracts property from single frame in LAMMPS data file or from
       slected or all frames in a NetCDF trajectory.
    Parameters
    ----------
    topology_file : str, optional
        LAMMPS data file to define topology (and possibly first frame).
    trajectory_file: str, optional
        Trajectory file (atom positions) in NetCDF format.
    frames: list of int, optional
        Arbitrary selection of frames. Per default, all frames if 
        trajectory_file not None, otherwise single frame.
    property_func: function(pipeline, frames)
        A function handle returning a single numpy array derived from
        the passed ovitos pipeline and list of int frame selection.
    modifier_func: (list of) function(pipeline)
        All functions in this list are applied to the pipeline 
        in the specified order before evaulation (but after loading the
        trajectory, if any has been specified).
    outfile: str, optional
        The property in form ao a numpy arrey is written to this text file.
    atom_style: str, optional
        LAMMPS atom sytle as understood by ovitos.
    """

    #logger = logging.getLogger(__name__)
    logger.info("Reading topology from '{:s}'".format(topology_file))

    # Load static topology data from a LAMMPS data file.
    pipeline = import_file(topology_file, atom_style = atom_style)
    
    if trajectory_file:
      logger.info("Reading trajectory from '{:s}'".format(trajectory_file))
      # Load atom trajectories from separate LAMMPS dump file.
      traj_mod = LoadTrajectoryModifier()
      traj_mod.source.load(trajectory_file, multiple_frames = True)

      # Insert modifier into modification pipeline.
      pipeline.modifiers.append(traj_mod)
      
      # if a trajectory is loaded, but no frames selected, then select all
      # frames per default:
      if frames is None: 
        frames = range(traj_mod.source.num_frames)
        logger.info("Selected all {:d} frames from trajectory '{:s}'.".format(
          len(frames), trajectory_file))

    # apply modifications if desired:
    if type(modifier_func) is not list:
      modifier_func = [modifier_func]
    for mf in modifier_func:
      if mf is not None:
        logger.info("Applying modifier function {:s}()...'.".format(mf.__name__))
        mf(pipeline)
      
    # if NO trajectory is loaded and no frames selected, then select '0'th frame:   
    if frames is None:
      frames = [0]

    logger.info("Selected {:d} frames.".format(len(frames)))

    property = property_func(pipeline, frames)  
    # otherwise, just pass None to evaluate "single frame" LAMMPS data file property
    if property is not None:
      np.savetxt(outfile, property)   
      logger.info("Property extracted by function {:s}() written to '{:s}'.".format(
        property_func.__name__, outfile))
    else:
      logger.error("Property function {:s}() returned nothing!".format(property_func.__name__))

    return
