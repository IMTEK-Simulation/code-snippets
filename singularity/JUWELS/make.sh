# Build all base images for NEMO

singularity build -F mpibase.sif openmpi-4.1.0rc1_ucx-1.8.1.def
singularity build -F fftw_and_io.si ../pnetcdf-1.12.2_fftw-3.3.9.def
singularity build -F muspectre.sif ../muspectre.def
