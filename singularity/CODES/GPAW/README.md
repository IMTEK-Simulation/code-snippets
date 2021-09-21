# GPAW

## Docker

This directory contains a Docker recipe that builds [GPAW](https://wiki.fysik.dtu.dk/gpaw/) with the Intel HPC (legacy) compiler suite. The Docker recipe uses a two stage build to minimize the image size. It compiles both the Python interface and the standalone code.

Build the Docker image with:
```bash
docker build -t gpaw .
```

The image is on Docker Hub [here](https://hub.docker.com/repository/docker/pastewka/gpaw).

## Singularity

The Docker images can be converted into a Singularity image for use on HPC systems. Convert the image with:
```bash
sudo singularity build gpaw.sif docker-daemon://gpaw:latest
```

## Docker Hub

The image is also available on [Docker Hub](https://hub.docker.com). To get it from there, execute:
```bash
docker pull pastewka/gpaw:210817
sudo singularity build gpaw.sif docker-daemon://pastewka/gpaw:210817
```

## Troubleshooting

### Performance issues

The container has a simple MPI benchmark installed. The exectuable is located at `/usr/local/bin/mpiBench`.

### Debugging libfabric issues

Set `FI_LOG_LEVEL=debug` before running the code.

### OMP: Error #179

If you get an error
```
OMP: Error #179: Function Can't open SHM2 failed:
OMP: System error #2: No such file or directory
```
you need to bind `/run/shm` into the container. Try executing:
```
OMP_NUM_THREADS=4 singularity run --bind /run/shm:/run/shm dftb.sif mdcore-1.0.1
```
