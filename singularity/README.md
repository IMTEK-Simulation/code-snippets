Singularity recipes
===================

As a general rule, the singularity recipe needs to contains the same version
of OpenMPI that is running on the target cluster. This required because process
distribution is orchestrated by the OpenMPI installation that resides on the
machine, while the code is run via the installation in the container.
Communication between these instances only works reliably if they have the same
version.

NEMO
----

NEMO uses OpenMPI 4.0.2. Please make sure to load the module prior to running
the container via `ml mpi/openmpi/4.0-gnu-9.2`.

Use the `openmpi-4.0.2_ucx_1.9.0_gcc-9-ubuntu_18.def` base container for
running on NEMO.