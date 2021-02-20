# Build all base images for NEMO

sudo singularity build -F mpibase.sif openmpi-4.0.2_psm2-11.2.185.def
sudo singularity build -F io.sif ../pnetcdf-1.12.2_fftw-3.3.9.def
sudo singularity build -F mpibench.sif ../mpibench.def
sudo singularity build -F muspectre.sif ../muspectre.def
