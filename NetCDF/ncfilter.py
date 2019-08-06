#! /usr/bin/env python

# ncfilter.py
#
# Copyright (C) 2019 IMTEK Simulation
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

# TODO:
# 2019/07/18: Hangs when more tasks than discrete entries in paralelization dim.
#             Porbably issue with MPI output.
'''
Filters a LAMMPS NetCDF trajectory by atom indices from GROMACS sytle .ndx file.

LAMMPS command "group2ndx" creates such a file from groups defined in LAMMPS if
LAMMPS has been compilesd with module "USER-COLVARS".
(refer to https://lammps.sandia.gov/doc/group2ndx.html)

Example: 

  mpirun -n 20 ncfilter.py --log ncfilter.log \\ 
      default.nc filtered.nc groups.ndx indenter

will select all atoms in the 'indenter' group as defined in 'groups.ndx'
by atom ids from 'default.nc' and write a copy of this subset to 'filtered.nc'.
While only error are logged to the terminal, specifying a log file will hold
messages on all verbosity levels.

Requires GromacsWrapper, tested with version 0.8.0 (April 30, 2019)
* https://gromacswrapper.readthedocs.io
* https://github.com/Becksteinlab/GromacsWrapper

Requires MPIFileHandler for logging,
* https://gist.github.com/chengdi123000/42ec8ed2cbef09ee050766c2f25498cb

Source: https://github.com/IMTEK-Simulation/code-snippets
'''
import os
import logging
import numpy as np
from mpi4py import MPI
from netCDF4 import Dataset
# tested with GromacsWrapper 0.8.0 (April 30, 2019)
from gromacs.fileformats.ndx import NDX

# https://gist.github.com/chengdi123000/42ec8ed2cbef09ee050766c2f25498cb
from MPIFileHandler import MPIFileHandler

# global settings
standard_loglevel   = logging.ERROR
standard_logformat  = ''.join(("%(asctime)s:%(name)s:%(levelname)s",
  "[ %(filename)s:%(lineno)s - %(funcName)s() ]: %(message)s"))

def ncfilter(
    infile      = 'default.nc',
    outfile      = 'filtered.nc',
    ndx_file    = "groups.ndx",
    g_sel       = "indenter",
    selvarname  = 'id',    # select subset based on NetCDF variable 'id'
    seldimname  = "atom",  # dimension to resize
    pardimname  = "frame", # dimension for chunk-wise parallelization
    format      = 'NETCDF4',
    collective  = True, # NetCDF MPIIO mode: collective or independent
    loglvl      = standard_loglevel,
    logfmt      = standard_logformat,
    logout      = None
  ):
  """Filters NetCDF trajectory by atom indices from GROMACS sytle .ndx file.

  Args:
      infile (str, optional):     NetCDF trajectory, defaults to 'default.nc'.
      outfile (str, optional):    NetCDF trajectory, defaults to 'filtered.nc'.
      ndx_file (str, optional):   GROMACS style .ndx file,
                                  defaults to 'groups.ndx'.
      g_sel (str, optional):      Name of group to select, as found in .ndx
                                  file, defaults to 'indenter'.
      selvarname (str, optional): Name of NetCDF variable to compare against
                                  indices in .ndx file, defaults to 'id'.
      seldimname (str, optional): Name of NetCDF dimension to reduce. Defaults
                                  to 'atom'.
      pardimname (str, optional): Name of NetCDF dimension to chop into chunks
                                  for parallel processing via MPI. Defaults to
                                  'frame'.
      format (str, optional):     NetCDF format, defaults to 'NETCDF4'.
      collective (bool, optional):Parallel IO mode, defaults to 'True'.
                                  Independent mode is used and unlimited are
                                  written as finite dimensions if 'False'.
                                  Please refert to
          http://unidata.github.io/netcdf4-python/netCDF4/index.html#section13
      loglvl (int, optional):     Log level, defaults to 'logging.ERROR'.
      logfmt (str, optional):     Override default log format.
      logout (str, optional):     Name of log file. Defaults to 'None',
                                  in which case the log streams to the terminal.

  Returns:
      Nothing.
  """


  # MPI communicator
  comm = MPI.COMM_WORLD

  size = comm.Get_size()
  rank = comm.Get_rank()

  # set up logging to terminal (and to file)
  # use same format for file and stream
  logger        = logging.getLogger("rank[{:03d}]:{}".format(rank,__name__))
  logger.setLevel(loglvl)

  logformatter  = logging.Formatter(logfmt)

  ch = logging.StreamHandler()
  ch.setFormatter(logformatter)
  ch.setLevel(loglvl)
  logger.addHandler(ch)

  if isinstance(logout,str):
    logger.setLevel(logging.DEBUG)
    fh = MPIFileHandler(logout)
    fh.setFormatter(logformatter)
    fh.setLevel(logging.DEBUG) # always verbose to log file
    logger.addHandler(fh)

  logger.debug('Hello from rank {}/{}.'.format(rank, size))

  # read groups from .ndx file:
  ndx = NDX()

  if rank == 0:
    logger.info(
      "Looking for '{:s}' in current working directory '{:s}'...".format(
        ndx_file, os.getcwd()) )

  ndx.read(ndx_file)

  if rank == 0:
    logger.info("Read {:d} groups from '{:s}':".format(len(ndx),ndx_file))
    for g in ndx.keys():
      logger.info("{:48s} - {: 24d} atoms".format(g,len(ndx[g])))

  # https://unidata.github.io/netcdf4-python/netCDF4/index.html#section15 You
  # can create netCDF Datasets whose content is held in memory instead of in a
  # disk file. There are two ways to do this. If you don't need access to the
  # memory buffer containing the Dataset from within python, the best way is to
  # use the diskless=True keyword argument when creating the Dataset. If you
  # want to save the Dataset to disk when you close it, also set persist=True.
  with Dataset(infile, 'r', parallel=True, comm=comm, info=MPI.Info(),
      format=format) as ncin, \
    Dataset(outfile, 'w', parallel=True, comm=comm, info=MPI.Info(),
      format=format, diskless=True, persist=True) as ncout:

    # get atom ids from netcdf, assume unsorted
    selvar = ncin.variables[selvarname]

    # look for parallelization dimension in selection variable:
    selvardnames = np.array(selvar.dimensions,dtype=str)
    selvarndims = len(selvardnames) # number of dimensions in selvar

    # select everything in selection variable to begin with
    selvarsel = [slice(None)]*selvarndims

    # check for parallelization along pardimname
    if rank == 0:
      logger.debug(
        "Parallelization dimension '{}' in selection variable {}{} ?".format(
          pardimname,selvarname,selvardnames))

    selvarpardimpos = np.argwhere(selvardnames == pardimname).flatten()
    if rank == 0:
      logger.debug("np.argwhere returned {}".format(selvarpardimpos))

    # in the standard case, selection bases on variable 'id' with dimensions
    # ('frame','atom'), where processing ist parallelized chunk-wise over
    # 'frame', while only a subset is processed over 'atom'.
    if len(selvarpardimpos) != 1:
      if rank == 0:
        logger.error(
          "Selection variable {}{} has no parallelization dimension {}!".format(
            selvarname, selvar.shape, pardimname) )
        raise ValueError()
    else:
      selvarpardimpos = selvarpardimpos[0]
      if rank == 0:
        logger.info(' '.join(("Selection variable {}{} has parallelization",
          "dimension {} at position {}")).format(
            selvarname, selvar.shape, pardimname, selvarpardimpos))

    # create a reference selection based on the first index in pardim
    # to assure the total numbert of particles not to change along pardim
    selvarsel[selvarpardimpos] = 0

    # TODO: split ndx-specific functionality from this function
    # with the ndx object here, selvar is expected to be only 1d after reduction
    # by parallelization dimension:
    refsel    = np.isin(selvar[selvarpardimpos], ndx[g_sel])
    refnumsel = np.count_nonzero(refsel)
    # number of selected entries along selection dimension

    if rank == 0:
      logger.info(' '.join(("Selection variable {}{} has {} entries",
        "along selection dimension {} in selection '{}'")).format(
          selvarname, selvarsel, refnumsel, seldimname, g_sel ) )

    # discretize trajectory along pardim and process chunk-wise
    pardimlen = len(ncin.dimensions[pardimname])
    if rank == 0:
      logger.info('Parallelization dimension: {}[{}]'.format(
        pardimname, pardimlen))

    if pardimlen < size: # more ranks than pardimlen needed!
        logger.error(' '.join(('Number of ranks {} exceeds size of',
           'parallelization dimension: {}[{}]. Reduce!')).format(
             size, pardimname, pardimlen))
        raise ValueError()

    n1 = rank*(pardimlen//size)
    n2 = (rank+1)*(pardimlen//size)
    logger.info('Parallelization dimension: {}[{}]'.format(
      pardimname, pardimlen))

    # TODO: this intended treatment of pardimlen < size does not work
    # treatment for special case where pardimlen < size
    # if n1 == 0 and n2 == 0:
    #   if rank == 0:
    #     logger.warn(' '.join(('Number of ranks {} exceeds size of',
    #       'parallelization dimension: {}[{}]')).format(
    #         size, pardimname, pardimlen))
    #   # treatment for special case where rank >= size
    #   if rank < pardimlen:
    #     n1 = rank  
    #     n2 = rank+1
    #   else:
    #     logger.warn(' '.join(('Apparently, rank {}/{} exceeds parallelization',
    #       'dimension {}[{}] and will thus be idle.')).format(
    #         rank, size, pardimname, pardimlen))
    #     n1 = 0
    #     n2 = 0
    # elif rank == size-1: # treatment for last rank if pardimlen >= size
    #   n2 = pardimlen
    
    if rank == size-1: # treatment for last rank if pardimlen >= size
      n2 = pardimlen


    logger.info(
      'Rank {}/{} treats parallelization dimension slice {}[{}:{}]'.format(
        rank, size, pardimname, n1, n2))

    # create global selection along all dimensions in selection variable
    # selection_shape = [ len(dim) for dim in selvar.get_dims() ]
    selection = np.zeros(selvar.shape, dtype=bool)

    # TODO: split ndx-specific functionality from this function
    for i in range(n1, n2):
      # select only selvar values that are within the specified index group
      selvarsel[selvarpardimpos] = i
      selvari = selvar[selvarsel]
      selection[selvarsel]  = np.isin(selvari, ndx[g_sel])
      numsel = np.count_nonzero(selection[selvarsel])
      logger.debug(' '.join(("Selection variable {}{} has {} entries along",
        "selection dimension {} in selection '{}'")).format(
          selvarname, selvarsel, numsel, seldimname, g_sel ))
      if numsel != refnumsel:
        logger.error(' '.join(("Selection variable {}{} has {} entries along",
          "selection dimension {} in selection '{}', differing from {}",
          "in reference selection!")).format(
          selvarname, selvarsel, numsel, seldimname, g_sel, refnumsel ))
        raise ValueError()
      elif np.any( np.not_equal( selection[i], refsel ) ):
        logger.warn(' '.join(("Selection variable {}{} and reference selection",
          "both have {} entries along selection dimension {}",
          "in selection '{}', but in differing order.")).format(
            selvarname, selvarsel, numsel, seldimname, g_sel ))

    # copy attributes, src:
    # https://stackoverflow.com/questions/15141563/
    # python-netcdf-making-a-copy-of-all-variables-and-attributes-but-one
    if rank == 0:
      logger.info("Copying global attributes...")

    for aname in ncin.ncattrs():
      if rank == 0:
        logger.info("Copy attibute {}".format(aname))
      ncout.setncattr(aname, ncin.getncattr(aname))

    # copy dimensions, src: https://gist.github.com/guziy/8543562
    if rank == 0:
      logger.info("Copying dimensions...")
    for dname, dim in ncin.dimensions.items():
      if dname == seldimname:
        lendim = refnumsel # reduce size along filter dimension
        if rank == 0:
          logger.info("Shrink dimension {}[{}] to [{}]".format(
            dname, len(dim), lendim))
      else:
        if collective:
          lendim = len(dim) if not dim.isunlimited() else None
        else:
          lendim = len(dim) # ignore unlimited dimension in this case

        if rank == 0:
          logger.info("Copy dimension {}[{}] to {}[{}]".format(
            dname, len(dim), dname, lendim))

      ncout.createDimension(dname, lendim)

    # copy variables
    if rank == 0:
      logger.info("Copying variables...")

    for vname, inVar in ncin.variables.items():
      logger.debug(
        "Going to copy innput variable {}{}: {}. Rank {}/{} still alive".format(
          vname, inVar.shape, inVar.datatype, rank, size))
      outVar = ncout.createVariable(vname, inVar.datatype, inVar.dimensions)
      if rank == 0:
        logger.info(
          "Copy input variable {}{}: {} to output variable {}{}: {}".format(
            vname, inVar.shape, inVar.datatype,
            vname, outVar.shape, outVar.datatype))
      if collective:
        outVar.set_collective(True) # output in collective mode
      else:
        outVar.set_collective(False) # output in independent mode

      # Copy variable attributes
      outVar.setncatts({k: inVar.getncattr(k) for k in inVar.ncattrs()})

      # look at dimensions in variable
      dnames = np.array(inVar.dimensions,dtype=str)
      ndims = len(dnames) # number of dimensions in current variable

      # select everything in variable to begin with
      in_selection = [slice(None)]*ndims
      out_selection = [slice(None)]*ndims

      # check for selection dimension:
      if rank == 0:
        logger.debug(' '.join(("Is the selection dimension '{}'",
          "in dimensions {} of variable {}{}?")).format(
          seldimname, dnames, vname, inVar.shape))
      seldimpos = np.argwhere(dnames == seldimname).flatten()
      if rank == 0:
        logger.debug("np.argwhere(...) returned '{}'".format(seldimpos))
      if len(seldimpos) != 1:
        seldimpos = None
        if rank == 0:
          logger.info(
            "Input variable {}{} has no selection dimension {}.".format(
              vname, inVar.shape, seldimname))
      else:
        seldimpos = seldimpos[0]
        if rank == 0:
          logger.info(' '.join(("Input variable {}{} has selection",
            "dimension {} at position {}")).format(
            vname, inVar.shape, seldimname, seldimpos))

      # check for parallelization dimension:
      if rank == 0:
        logger.debug(' '.join(("Is the parallelization dimension '{}'",
          "in dimensions {} of variable {}{}?")).format(
          pardimname, dnames, vname, inVar.shape))

      pardimpos = np.argwhere(dnames == pardimname).flatten()
      if rank == 0:
        logger.debug("np.argwhere(...) returned {}".format(pardimpos))

      # no parallelization dimension in current inVar:
      if len(pardimpos) != 1:
        if rank == 0:
          logger.info(
            "Input variable {}{} has no parallelization dimension {}.".format(
              vname, inVar.shape, pardimname))

        # in the current implementation, the selection dimension is expected
        # only to occur in variables that have parallelization dimension as well
        if (seldimpos is not None) and (rank == 0):
          logger.error(' '.joint((
            "Input variable {}{} has selection dimension {},"
            "but no parallelization dimension {}!")).format(
            vname, inVar.shape, seldimname, pardimname))

        if rank == 0:
          logger.info("Copy input variable {}{} as is at once.".format(
            vname, inVar.shape))
        outVar[:] = inVar[:]

      # parallelization dimension found in current inVar:
      else:
        pardimpos = pardimpos[0]
        if rank == 0:
          logger.info(' '.join((
            "Input variable {}{} has parallelization dimension {}",
            "at position {}")).format(
              vname, inVar.shape, pardimname, pardimpos))
        
        if collective:
           # issue: outVar.shape has zero entry for unlimited dim
          if rank == 0:
            logger.info("outVar {}{} has {} unlimited dims".format( 
              vname, outVar.shape, np.equal( outVar.shape, 0) ))
          tmpVar_shape = np.where( np.equal( outVar.shape, 0 ), inVar.shape, outVar.shape)
          tmpVar = np.zeros( tmpVar_shape, dtype=outVar.dtype )
        else:
          tmpVar = np.zeros( outVar.shape, dtype=outVar.dtype )

        if rank == 0:
          logger.debug("np.array tmpVar{}: {} created".format(
            tmpVar.shape, tmpVar.dtype))

        # https://unidata.github.io/netcdf4-python/netCDF4/index.html#section6
        # Boolean array and integer sequence indexing behaves differently for
        # netCDF variables than for numpy arrays. Only 1-d boolean arrays and
        # integer sequences are allowed, and these indices work independently
        # along each dimension (similar to the way vector subscripts work in
        # fortran).

        # Here, we iterate subsequently over temporal axis as selection may
        # vary for every frame and multidimensional selection arrays are no
        # option.
        for i in range(n1, n2):
          if seldimpos is not None:
            selvarsel[selvarpardimpos] = i
            in_selection[seldimpos] = selection[selvarsel]

          in_selection[pardimpos]  = i
          out_selection[pardimpos] = i

          logger.debug("Assigning input variable {}{} to tmpVar{}".format(
              vname, in_selection, out_selection))

          tmpVar[out_selection] = inVar[in_selection]

          logger.debug("Assigned input variable {}{} to tmpVar{}".format(
              vname, in_selection, out_selection))

        # TODO:
        # this commented conditional clause intends to alleviate issue when 
        # number of ranks exceeds number of prallelizable entried in par dim.
        # if n2 > n1: # only if there is a finite range of elements to assign:
        out_selection[pardimpos] = slice(n1,n2)
        logger.info(' '.join(("Assigning filtered input variable {}{}",
          "to output var {}{}")).format(
            vname, out_selection, vname, out_selection))

        # apparently, accumulating the desired subset index by index in an numpy
        # array tmpVar and then assigning the accumulated chunk to the NetCDF
        # variable outVar works better than assigning subset to NetCDF variable
        # index by index directly
        outVar[out_selection] = tmpVar[out_selection]
        logger.info(' '.join(("Assigned filtered input variable {}{}",
          "to output var {}{}")).format(
            vname, out_selection, vname, out_selection))
        # else:
        #  logger.info("Skipped input variable {} on rank {}/{}.".format(
        #    vname, rank, size))
          
    # nc.close() not necessary due to 'with' statement
  logger.debug('Goodbye from rank {}/{}.'.format(rank, size))

def main():
  import argparse

  # in order to have both:
  # * preformatted help text and ...
  # * automatic display of defaults
  class ArgumentDefaultsAndRawDescriptionHelpFormatter(
      argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass

  parser = argparse.ArgumentParser(description=__doc__,
      formatter_class = ArgumentDefaultsAndRawDescriptionHelpFormatter)

  parser.add_argument(
    'infile', default='default.nc', help='NetCDF trajectory')
  parser.add_argument(
    'outfile', default='filtered.nc', help='NetCDF trajectory')
  parser.add_argument(
    'ndx', default='group.ndx', help='GROMACS style .ndx file')
  parser.add_argument(
    '-n', '--ndx', default='group.ndx', dest='ndx',
    help='GROMACS style .ndx file')
  parser.add_argument(
    'group', default='indenter', help='Group selection as named in .ndx')
  parser.add_argument(
    '-g','--group', default='indenter', dest='group',
    help='Group selection as named in .ndx')
  parser.add_argument(
    '-sv','--selection-variable', default='id', dest='selvarname',
    help='Name of NetCDF variable to compare against indices in .ndx file')
  parser.add_argument(
    '-sd','--selection-dimension', default='atom', dest='seldimname',
    help='Name of NetCDF dimension to reduce')
  parser.add_argument(
    '-pd','--parallelization-dimension', default='frame', dest='pardimname',
    help='Name of NetCDF dimension for chunk-wise parallel processing.')
  parser.add_argument(
    '-f','--fmt', '--format', default='NETCDF4', dest='format',
    help=' '.join(('NetCDF format. Refer to',
    'http://unidata.github.io/netcdf4-python/netCDF4/index.html#section1'))
    )
  parser.add_argument(
    '-i','--independent', action='store_false', dest='collective',
    help=' '.join(('Parallel IO mode. Independent mode is used',
      'and unlimited are written as finite dimensions if specified.',
      'Otherwise, collective mode is used by default. Please refert to',
      'http://unidata.github.io/netcdf4-python/netCDF4/index.html#section13'))
    )
  parser.add_argument(
    '-l','--log', default=None, const='log.out', dest='logfile', nargs='?',
    help=' '.join(("Write output to 'log.out' or any other log file whose name",
      "is specified optionally after this flag, instead of stream output to",
      "the terminal"))
    )
  parser.add_argument('--verbose', '-v', action='count', dest='verbose',
      default=0, help='Make this tool more verbose')
  parser.add_argument('--debug','-d', action='store_true',
      help='Make this tool print debug info')
  args = parser.parse_args()

  loglevel = logging.ERROR
  
  if args.verbose > 0:
    loglevel = logging.WARN
  if args.verbose > 1:
    loglevel = logging.INFO
  if args.debug or (args.verbose > 2):
    loglevel = logging.DEBUG

  ncfilter(
    infile      = args.infile,
    outfile      = args.outfile,
    ndx_file    = args.ndx,
    g_sel       = args.group,
    selvarname  = args.selvarname,
    seldimname  = args.seldimname,
    pardimname  = args.pardimname,
    format      = args.format,
    collective  = args.collective,
    loglvl      = loglevel,
    logout      = args.logfile
  )

if __name__ == '__main__':
  main()
