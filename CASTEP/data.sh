#!/bin/bash
# 
sed -n '/Smax/p' $1 | awk '{print $4}' > SRelax.dat
sed -n '/|F|/p' $1 | awk '{print $4}' > FRelax.dat
sed -n '/dE\/ion/p' $1 | awk '{print $4}' > ERelax.dat
#
sed -n '/alpha = /p' $1 | awk '{print $3}' > a.dat
sed -n '/beta  = /p' $1 | awk '{print $3}' > b.dat
sed -n '/gamma = /p' $1 | awk '{print $3}' > c.dat