# Build all base images for NEMO

singularity build -F mpibase.sif openmpi-4.0.2_psm2-11.2.185.def
singularity build -F fftw_and_io.sif ../pnetcdf-1.12.2_fftw-3.3.9_pfft-master.def
singularity build -F muspectre.sif ../muspectre.def
