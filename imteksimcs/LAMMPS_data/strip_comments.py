#!/usr/bin/env python
"""
Strips all comments from LAMMPS data file to be processed by pizza.py
"""

def strip_comments(infile, outfile):
    """Removes all trailing comments from a LAMMPS data file.
       Necessary to make them pizza.py-processible"""
    import re
    regex = re.compile(r"\s*#.*$")
    with open(infile) as i:
        with open(outfile, 'w') as o:
            for line in i:
                line = regex.sub('',line)
                o.write(line)
