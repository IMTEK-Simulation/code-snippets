# muSpectre

## Docker

This directory contains a Docker recipe that builds [muSpectre](https://gitlab.com/muspectre/muspectre) with the Intel HPC (legacy) compiler suite. The Docker recipe uses a two stage build to minimize the image size. It compiles both the Python interface and the standalone code.

Build the Docker image with:
```bash
docker build -t muspectre .
```

You can run a shell inside the container with:
```bash
docker run -it --mount type=bind,source="${HOME}",target="/home" muspectre /bin/bash
```
Note that this binds your home directory to `/home`.

The image is on Docker Hub [here](https://hub.docker.com/repository/docker/pastewka/muspectre).

## Singularity

The Docker images can be converted into a Singularity image for use on HPC systems. Convert the image with:
```bash
sudo singularity build muspectre.sif docker-daemon://muspectre:latest
```

## Docker Hub

The image is also available on [Docker Hub](https://hub.docker.com). To get it from there, execute:
```bash
docker pull pastewka/muspectre:210817
sudo singularity build muspectre.sif docker-daemon://pastewka/muspectre:210817
```

## Troubleshooting

### Performance issues

The container has a simple MPI benchmark installed. The exectuable is located at `/usr/local/bin/mpiBench`.

### Debugging libfabric issues

Set `FI_LOG_LEVEL=debug` before running the code.
