# Build all base images for NEMO

sudo singularity build mpibase.sif openmpi-4.0.2_psm2-11.2.185_fedora-33.def
sudo singularity build io.sif ../pnetcdf-1.12.2_fftw-3.3.9_openmpi-4.0.2_gcc-7_ubuntu-18.def
sudo singularity build mpibench.sif ../mpibench.def
