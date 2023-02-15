# Build all base images for NEMO

#singularity build -F mpibase.sif openmpi-4.0.2_psm2-11.2.185.def
singularity build -F nemo_mpi_fftw_and_io.sif ../../STACK/pnetcdf-1.12.2_fftw-3.3.9.def

