# GPAW

This directory contains a Docker recipe that builds [GPAW](https://wiki.fysik.dtu.dk/gpaw/) with the Intel HPC (legacy) compiler suite. The Docker recipe uses a two stage build to minimize the image size. It compiles both the Python interface and the standalone code.

Build the Docker image with:
```bash
docker build -t gpaw .
```

# Singularity

The Docker images can be converted into a Singularity image for use on HPC systems. You need to export the image first and then convert it to Singularity.

Convert the image with:
```bash
docker save gpaw -o gpaw.tar
singularity build gpaw.sif docker-archive://gpaw.tar
```

# Troubleshooting

If you get an error
```
OMP: Error #179: Function Can't open SHM2 failed:
OMP: System error #2: No such file or directory
```
you need to bind `/run/shm` into the container. Try executing:
```
OMP_NUM_THREADS=4 singularity run --bind /run/shm:/run/shm dftb.sif mdcore-1.0.1
```