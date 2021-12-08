# DFTB

This directory contains a Docker recipe that builds [Atomistica](https://github.com/Atomistica/atomistica) and [DFTB+](https://github.com/dftbplus/dftbplus) with the Intel HPC (legacy) compiler suite. The Docker recipe uses a two stage build to minimize the image size. It compiles both the Python interface and the standalone code.

Build the Docker image with:
```bash
docker build -t dftb .
```

# Singularity

The Docker images can be converted into a Singularity image for use on HPC systems. You need to export the image first and then convert it to Singularity.

Convert the image with:
```bash
docker save dftb -o dftb.tar
singularity build dftb.sif docker-archive://dftb.tar
```

## Running the standalone code

To run the standalone code, execute:

```bash
OMP_NUM_THREADS=4 singularity run dftb.sif mdcore-1.0.1
```

If you get an error
```
OMP: Error #179: Function Can't open SHM2 failed:
OMP: System error #2: No such file or directory
```
you need to bind `/run/shm` into the container. Try executing:
```
OMP_NUM_THREADS=4 singularity run --bind /run/shm:/run/shm dftb.sif mdcore-1.0.1
```

## Running the Python interface

The Python interface can run along the same lines:
```bash
OMP_NUM_THREADS=4 singularity run --bind /run/shm:/run/shm dftb.sif python3 my_python_script.py
```
