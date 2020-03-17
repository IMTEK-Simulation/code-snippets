#!/bin/bash
set -e

# This is a sample on how to use selectByAtomType.nco
# perf stat -B prints performance statistics for each step

# The -D5 option makes NetCDF operators print level five debug information

# delete unwanted quantities
echo "01 -- DELETE"
perf stat -B ncks -D5 -O -x -v c_peratom_stress,f_peratom_stress_ave,mass,mol,velocities $1 _thinned_$2

# Note that renaming a dimension to the name of a dependent variable can be used to invert the relationship between an independent coordinate variable and a dependent variable.  In  this  case,
#       the named dependent variable must be one-dimensional and should have no missing values.  Such a variable will become a coordinate variable.

# permute dimensions (apparently, this is the most demanding step)
echo "02 -- PERMUTE"
perf stat -B ncpdq -D5 -O -a atom,frame,spatial _thinned_$2 _thinned_permuted_sorted_$2

# the following steps are not necessary if NetCDF is already sorted by atom id

# rename dimension
# ncrename -d atom,id _thinned_permuted_$2 _thinned_permuted_renamed_$2

# sort NetCDF by atom ID
# ncap2 -O -v -S sortByAtomID.nco _thinned_permuted_$2 _thinned_permuted_sorted_$2

# select based on variable-dependent criterion
# You must sort and copy all variables of interest explicitly within
# selectByAtomType.nco! Otherwise, they will not appear in the output.
echo "03 -- SELECT"
perf stat -B ncap2 -D5 -O -v -S selectByAtomType.nco _thinned_permuted_sorted_$2 _thinned_permuted_sorted_selected_$2

# drop old atom coordinate
echo "04 -- DROP"
perf stat -B ncks -D5 -O -x -v atom _thinned_permuted_sorted_selected_$2 _thinned_permuted_sorted_selected_deleted_$2

# rename new atom coordinate to
echo "05 -- RENAME"
perf stat -B ncrename -D5 -O -d atom_out,atom _thinned_permuted_sorted_selected_deleted_$2 _thinned_permuted_sorted_selected_deleted_renamed_$2

# permute back
echo "06 -- PERMUTE"
perf stat -B ncpdq -D5 -O -a frame,atom,spatial _thinned_permuted_sorted_selected_deleted_renamed_$2 $2

# now, you have a NetCDF trajectory of your designated atom subset only
