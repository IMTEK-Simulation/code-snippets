#!/usr/bin/env python
"""
Converts LAMMPS pair_coeff lines to 'pair_style hybrid' format.
"""
# 29 Mar 2019, converts LAMMPS force field pair_coeff entries to hybrid style


def convert(infile, outfile, pair_style, regular_expression):
  import re
  regex = re.compile(regular_expression)
  with open(infile,'r') as f, open(outfile,'w') as g:
    for line in f:
      match = re.match(regex, line)
      if match is not None:
        line = '{:s} {:s} {:s}'.format(match.group(1), pair_style, match.group(2))
      g.write(line.rstrip() + '\n')
