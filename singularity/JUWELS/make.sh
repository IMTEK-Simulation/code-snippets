# Build all base images for NEMO

singularity build -F mpibase.sif openmpi-4.1.0_fedora-32.def
singularity build -F io.sif ../pnetcdf-1.12.2_fftw-3.3.9.def
singularity build -F mpibench.sif ../mpibench.def
sudo singularity build -F muspectre.sif ../muspectre.def
