# SelectionFromNdx.ovito-3.0.0-dev349.py
#
# Copyright (C) 2019 IMTEK Simulation
# Author: Johannes Hoermann, johannes.hoermann@imtek.uni-freiburg.de
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
Creates selection from GROMACS-type .ndx file.

LAMMPS command "group2ndx" creates such a file from groups defined in LAMMPS if
LAMMPS has been compilesd with module "USER-COLVARS".
(refer to https://lammps.sandia.gov/doc/group2ndx.html)

ndx file input functionality adopted from "GromacsWrapper" for convenience
(2019-06-06 git master, https://github.com/Becksteinlab/GromacsWrapper)

GPLv3 adopted from "GromacsWrapper" as well.
"""
#TODO: pass ndx_file and group name somewhere as argument
ndx_file    = "groups.ndx"
g_sel       = "indenter"

from ovito.data import *

import os
import re
import numpy
from collections import OrderedDict as odict

# source: https://github.com/Becksteinlab/GromacsWrapper/blob/master/gromacs/utilities.py
class FileUtils(object):
    """Mixin class to provide additional file-related capabilities."""

    #: Default extension for files read/written by this class.
    default_extension = None

    def _init_filename(self, filename=None, ext=None):
        """Initialize the current filename :attr:`FileUtils.real_filename` of the object.
        Bit of a hack.
        - The first invocation must have ``filename != None``; this will set a
          default filename with suffix :attr:`FileUtils.default_extension`
          unless another one was supplied.
        - Subsequent invocations either change the filename accordingly or
          ensure that the default filename is set with the proper suffix.
        """

        extension = ext or self.default_extension
        filename = self.filename(filename, ext=extension, use_my_ext=True, set_default=True)
        #: Current full path of the object for reading and writing I/O.
        self.real_filename = os.path.realpath(filename)

# source: https://github.com/Becksteinlab/GromacsWrapper/blob/master/gromacs/fileformats/ndx.py
class NDX(odict, FileUtils):
    """Gromacs index file.
    Represented as a ordered dict where the keys are index group names and
    values are numpy arrays of atom numbers.
    Use the :meth:`NDX.read` and :meth:`NDX.write` methods for
    I/O. Access groups by name via the :meth:`NDX.get` and
    :meth:`NDX.set` methods.
    Alternatively, simply treat the :class:`NDX` instance as a
    dictionary. Setting a key automatically transforms the new value
    into a integer 1D numpy array (*not* a set, as would be the
    :program:`make_ndx` behaviour).
    .. Note::
       The index entries themselves are ordered and can contain
       duplicates so that output from NDX can be easily used for
       :program:`g_dih` and friends. If you need set-like behaviour
       you will have do use :class:`gromacs.formats.uniqueNDX` or
       :class:`gromacs.cbook.IndexBuilder` (which uses
       :program:`make_ndx` throughout).
    **Example**
      Read index file, make new group and write to disk::
        ndx = NDX()
        ndx.read('system.ndx')
        print ndx['Protein']
        ndx['my_group'] = [2, 4, 1, 5]   # add new group
        ndx.write('new.ndx')
      Or quicker (replacing the input file ``system.ndx``)::
        ndx = NDX('system')          # suffix .ndx is automatically added
        ndx['chi1'] = [2, 7, 8, 10]
        ndx.write()
    """
    default_extension = "ndx"

    # match:  [ index_groupname ]
    SECTION = re.compile("""\s*\[\s*(?P<name>\S.*\S)\s*\]\s*""")

    #: standard ndx file format: 15 columns
    ncol = 15
    #: standard ndx file format: '%6d'
    format = '%6d'

    def __init__(self, filename=None, **kwargs):
        super(NDX, self).__init__(**kwargs)  # can use kwargs to set dict! (but no sanity checks!)

        if filename is not None:
            self._init_filename(filename)
            self.read(filename)

    def read(self, filename=None):
        """Read and parse index file *filename*."""
        self._init_filename(filename)

        data = odict()
        with open(self.real_filename) as ndx:
            current_section = None
            for line in ndx:
                line = line.strip()
                if len(line) == 0:
                    continue
                m = self.SECTION.match(line)
                if m:
                    current_section = m.group('name')
                    data[current_section] = []  # can fail if name not legal python key
                    continue
                if current_section is not None:
                    data[current_section].extend(map(int, line.split()))

        super(NDX,self).update(odict([(name, self._transform(atomnumbers))
                                     for name, atomnumbers in data.items()]))

    def write(self, filename=None, ncol=ncol, format=format):
        """Write index file to *filename* (or overwrite the file that the index was read from)"""
        with open(self.filename(filename, ext='ndx'), 'w') as ndx:
            for name in self:
                atomnumbers = self._getarray(name)  # allows overriding
                ndx.write('[ {0!s} ]\n'.format(name))
                for k in range(0, len(atomnumbers), ncol):
                    line = atomnumbers[k:k+ncol].astype(int)   # nice formatting in ncol-blocks
                    n = len(line)
                    ndx.write((" ".join(n*[format])+'\n') % tuple(line))
                ndx.write('\n')

    def get(self, name):
        """Return index array for index group *name*."""
        return self[name]

    def set(self, name, value):
        """Set or add group *name* as a 1D numpy array."""
        self[name] = value

    def size(self, name):
        """Return number of entries for group *name*."""
        return len(self[name])

    @property
    def groups(self):
        """Return a list of all groups."""
        return self.keys()

    @property
    def sizes(self):
        """Return a dict with group names and number of entries,"""
        return {name: len(atomnumbers) for name, atomnumbers in self.items()}

    @property
    def ndxlist(self):
        """Return a list of groups in the same format as  :func:`gromacs.cbook.get_ndx_groups`.
        Format:
           [ {'name': group_name, 'natoms': number_atoms, 'nr':  # group_number}, ....]
        """
        return [{'name': name, 'natoms': len(atomnumbers), 'nr': nr+1} for
                nr,(name,atomnumbers) in enumerate(self.items())]

    def _getarray(self, name):
        """Helper getter that is used in write().
        Override when using a _transform that stores something that
        cannot be indexed, e.g. when using set()s.
        """
        return self[name]

    def _transform(self, v):
        """Transform input to the stored representation.
        Override eg with ``return set(v)`` for index lists as sets.
        """
        return numpy.ravel(v).astype(int)

    def __setitem__(self, k, v):
        super(NDX, self).__setitem__(k, self._transform(v))

    def setdefault(*args,**kwargs):
        raise NotImplementedError


    def filename(self,filename=None,ext=None,set_default=False,use_my_ext=False):
        """Supply a file name for the class object.
        Typical uses::
           fn = filename()             ---> <default_filename>
           fn = filename('name.ext')   ---> 'name'
           fn = filename(ext='pickle') ---> <default_filename>'.pickle'
           fn = filename('name.inp','pdf') --> 'name.pdf'
           fn = filename('foo.pdf',ext='png',use_my_ext=True) --> 'foo.pdf'
        The returned filename is stripped of the extension
        (``use_my_ext=False``) and if provided, another extension is
        appended. Chooses a default if no filename is given.
        Raises a ``ValueError`` exception if no default file name is known.
        If ``set_default=True`` then the default filename is also set.
        ``use_my_ext=True`` lets the suffix of a provided filename take
        priority over a default ``ext`` tension.
        .. versionchanged:: 0.3.1
           An empty string as *ext* = "" will suppress appending an extension.
        """
        if filename is None:
            if not hasattr(self,'_filename'):
                self._filename = None        # add attribute to class
            if self._filename:
                filename = self._filename
            else:
                raise ValueError("A file name is required because no default file name was defined.")
            my_ext = None
        else:
            filename, my_ext = os.path.splitext(filename)
            if set_default:                  # replaces existing default file name
                self._filename = filename
        if my_ext and use_my_ext:
            ext = my_ext
        if ext is not None:
            if ext.startswith(os.extsep):
                ext = ext[1:]  # strip a dot to avoid annoying mistakes
            if ext != "":
                filename = filename + os.extsep + ext
        return filename

    def check_file_exists(self, filename, resolve='exception', force=None):
        """If a file exists then continue with the action specified in ``resolve``.
        ``resolve`` must be one of
        "ignore"
              always return ``False``
        "indicate"
              return ``True`` if it exists
        "warn"
              indicate and issue a :exc:`UserWarning`
        "exception"
              raise :exc:`IOError` if it exists
        Alternatively, set *force* for the following behaviour (which
        ignores *resolve*):
        ``True``
              same as *resolve* = "ignore" (will allow overwriting of files)
        ``False``
              same as *resolve* = "exception" (will prevent overwriting of files)
        ``None``
              ignored, do whatever *resolve* says
        """
        def _warn(x):
            msg = "File {0!r} already exists.".format(x)
            logger.warn(msg)
            warnings.warn(msg)
            return True
        def _raise(x):
            msg = "File {0!r} already exists.".format(x)
            logger.error(msg)
            raise IOError(errno.EEXIST, x, msg)
        solutions = {'ignore': lambda x: False,      # file exists, but we pretend that it doesn't
                     'indicate': lambda x: True,     # yes, file exists
                     'warn': _warn,
                     'warning': _warn,
                     'exception': _raise,
                     'raise': _raise,
                     }

        if force is True:
            resolve = 'ignore'
        elif force is False:
            resolve = 'exception'

        if not os.path.isfile(filename):
            return False
        else:
            return solutions[resolve](filename)

    def infix_filename(self, name, default, infix, ext=None):
        """Unless *name* is provided, insert *infix* before the extension *ext* of *default*."""
        if name is None:
            p, oldext = os.path.splitext(default)
            if ext is None:
                ext = oldext
            if ext.startswith(os.extsep):
                ext = ext[1:]
            name = self.filename(p+infix, ext=ext)
        return name

    def __repr__(self):
        fmt = "{0!s}(filename=%r)".format(self.__class__.__name__)
        try:
            fn =  self.filename()
        except ValueError:
            fn = None
        return fmt % fn

def modify(frame, data):
    ndx = NDX()
    print("Looking for '{:s}' in current working directory '{:s}'...".format(
        ndx_file,os.getcwd()))

    ndx.read(ndx_file)
    print("Read {:d} groups from '{:s}':".format(len(ndx),ndx_file))

    for g in ndx.keys():
        print("{:48s} - {: 24d} atoms".format(g,len(ndx[g])))

    selection = data.particles_.create_property('Selection')

    yield "Creating selection from group '{:s}' in '{:s}'...".format(
        g_sel, ndx_file)

    # source: http://www.ovito.org/manual_testing/python/introduction/examples.html#example-select-overlapping-particles
    with selection:
        for index, particle_id in enumerate(ndx[g_sel]):
            # Update progress display in the status bar.
            yield (index / ndx.size(g_sel))

            selection[
                data.particles['Particle Identifier'] == particle_id ] = 1

    print("Selected {:d} atoms.".format(numpy.count_nonzero(selection)))
