# Snippets for visualization software PyMol

https://github.com/schrodinger/pymol-open-source

PyMOL supports parallel rendering via multithreading. However, this involves
reading serially through the whole trajectory initially, which may take hours
for long ones. Furthermore, that does not allow for parallelization across
multiple nodes.

Since rendering each frame is a completely independent operation, splitting a
trajectory into chunks, i.e. with methods in `imteksimcs.GROMACS.gmx_split_traj`
and then rendering those chunks via independent processes is trivial and
implemented within `pymol_mpi`.
