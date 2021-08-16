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
