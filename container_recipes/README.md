Singularity recipes
===================

As a general rule, the singularity recipe needs to contains the same version
of OpenMPI that is running on the target cluster. This required because process
distribution is orchestrated by the OpenMPI installation that resides on the
machine, while the code is run via the installation in the container.
Communication between these instances only works reliably if they have the same
version.

_All_ manually compiled code is installed in `/usr/local`. All base container
contain a simple benchmark in `/opt/mpiBench`.

NEMO
----

NEMO uses OpenMPI 4.0.2. Please make sure to load the module prior to running
the container via `ml mpi/openmpi/4.0-gnu-9.2`.

Use the `openmpi-4.0.2_psm2-11.2.185.def` base container for running on NEMO.

_Notes_: NEMO uses OmniPath. While OpenMPI can support OmniPath through UCX,
PSM2 gives much better throughput in simple benchmarks.

An example script for NEMO looks as follows:
```bash
#MSUB -l nodes=1:ppn=2
#MSUB -l walltime=00:15:00
#MSUB -l pmem=1000mb
#MSUB -q express
#MSUB -m bea

set -e

ml system/modules/testing
ml mpi/openmpi/4.0-gnu-9.2
ml tools/singularity/3.5

cd $MOAB_SUBMITDIR

# Run in container
mpirun -n $PBS_NP singularity exec mpibase.sif /opt/mpiBench/mpiBench
```

JUWELS
------

JUWELS automatically import the host MPI to the container. This leads to
several problems, in particular with version of some libraries inside the
container. (JUWELS MPI uses some libraries from the JUWELS software stack
and some - such aus PMIX - from inside the container. Even worse, some
libraries have newer version than in the CentOS base system.).

Use the `openmpi-4.1.0rc1_ucx-1.9.0.def` base container for running on JUWELS.
You will need to copy the proprietary GPFS library and headers to your
computer before building the container. Execute
```bash
tar -cvzf gpfs.tar.gz -C /usr/lpp/mmfs lib/libgpfs.so include
```
on JUWELS and copy the file `gpfs.tar.gz` to the directory where you build the
base container.

You need to specify `--mpi=pspmix` with `srun`. Make sure that you _do not_
load MPI through the Lmod system. An example script on JUWELS looks as
follows:
```bash
#!/bin/bash -x
#SBATCH --account=hka18
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=00:30:00
#SBATCH --partition=batch

PROJECT_PATH=/p/project/chka18

module load GCC/9.3.0 Singularity-Tools

srun --mpi=pspmix singularity run ${PROJECT_PATH}/muspectre.sif ${PROJECT_PATH}/projects/SurfaceRoughness/Continuum/biaxial_free_surface.py -d 0.01 -n 10 -g 128,128,128
```
