# To generate a singularity container: singularity build image.sif docker-daemon://<image tag>:latest
FROM intel/oneapi:os-tools-ubuntu18.04 AS builder

ARG DEBIAN_FRONTEND=noninteractive

# Latest release of Atomistica; run docker build --build-arg ATOMISTICA_VERSION=master to fetch master
ARG ATOMISTICA_VERSION="1.0.4"

# Latest release of DFTB+; run docker build --build-arg DFTBPLUS_VERSION=main to fetch main
ARG DFTBPLUS_VERSION="22.1"

# Latest reslease of ASE; run docuer build --build-arg ASE_VERSION=master to fetch master
ARG ASE_VERSION="3.22.1"

# WARNING: if you change this version number, also change it below
ENV INTEL_VERSION="2022.1.0"

RUN apt-get update -y \
	&& apt-get install -y \
	software-properties-common \
	lsb-release \
	wget \
	&& wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | gpg --dearmor - | tee /etc/apt/trusted.gpg.d/kitware.gpg >/dev/null \
	&& apt-add-repository "deb https://apt.kitware.com/ubuntu/ $(lsb_release -cs) main" \
	&& apt-get update -y \
	&& apt-get install -y kitware-archive-keyring \
	&& rm /etc/apt/trusted.gpg.d/kitware.gpg \
	&& apt-get update -y \
	&& apt-get install -y \
	cmake \
	curl \
	git \
	cmake \
	gfortran \
	python3 \
	python3-dev \
	python3-pip \
	intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic \
	intel-oneapi-compiler-fortran \
	intel-oneapi-mkl-devel \
	zlib1g \
	zlib1g-dev \
	libjpeg-dev

# Setting environment variables
ENV PATH="/opt/intel/oneapi/mkl/${INTEL_VERSION}/bin/intel64:/opt/intel/oneapi/compiler/${INTEL_VERSION}/linux/bin/intel64:/opt/intel/oneapi/compiler/${INTEL_VERSION}/linux/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
ENV CPATH="/opt/intel/oneapi/mkl/${INTEL_VERSION}/include:/opt/intel/oneapi/compiler/${INTEL_VERSION}/linux/include"
ENV LD_LIBRARY_PATH="/opt/intel/oneapi/mkl/${INTEL_VERSION}/lib/intel64:/opt/intel/oneapi/compiler/${INTEL_VERSION}/linux/lib:/opt/intel/oneapi/compiler/${INTEL_VERSION}/linux/lib/x64:/opt/intel/oneapi/compiler/${INTEL_VERSION}/linux/compiler/lib/intel64_lin"
ENV LIBRARY_PATH="/opt/intel/oneapi/mkl/${INTEL_VERSION}/lib/intel64:/opt/intel/oneapi/compiler/${INTEL_VERSION}/linux/compiler/lib/intel64_lin:/opt/intel/oneapi/compiler/${INTEL_VERSION}/linux/lib"
ENV MKLROOT="/opt/intel/oneapi/mkl/${INTEL_VERSION}"

# Sanity check to see if PATH, etc. are correct
RUN icc --version

# Working around Ubuntu's 18.04's Python 3.6
# Numpy 1.19.5 and scipy 1.5.4 are last usable versions
RUN python3 -m pip install "setuptools<49.2.0" "cython==0.29.21"
ENV NUMPY_VERSION="1.19.5" SCIPY_VERSION="1.5.4"

# Fetching sources
WORKDIR /tmp
RUN ( curl -L https://github.com/numpy/numpy/archive/refs/tags/v${NUMPY_VERSION}.tar.gz | tar xz ) \
    && ( curl -L https://github.com/scipy/scipy/archive/refs/tags/v${SCIPY_VERSION}.tar.gz | tar xz )

# Installing Numpy
WORKDIR /tmp/numpy-${NUMPY_VERSION}
RUN python3 setup.py config --compiler=intelem build_clib --compiler=intelem build_ext --compiler=intelem install

# Installing Scipy
WORKDIR /tmp/scipy-${SCIPY_VERSION}
RUN python3 setup.py config --compiler=intelem --fcompiler=intelem build_clib --compiler=intelem --fcompiler=intelem build_ext --compiler=intelem --fcompiler=intelem install

WORKDIR /tmp/atomistica
RUN curl -L https://github.com/Atomistica/atomistica/tarball/${ATOMISTICA_VERSION} | tar -xz --strip-components=1

# Install Atomistica - Python module and standalone code
RUN python3 -m pip install git+https://gitlab.com/ase/ase.git@${ASE_VERSION} \
	&& cp setup.cfg.intel_omp setup.cfg \
	&& python3 setup.py install --prefix=/usr/local --force \
	&& cd build_standalone \
	&& cp Makefile.intel Makefile \
	&& DEBUG=0 OPENMP=1 make factories \
	&& DEBUG=0 OPENMP=1 make mdcore \
	&& mkdir -p /usr/local/bin \
	&& cp mdcore-* /usr/local/bin

WORKDIR /tmp/dftbplus
RUN curl -L https://github.com/dftbplus/dftbplus/tarball/${DFTBPLUS_VERSION} | tar xz --strip-components=1

# Build and install DFTB+
WORKDIR build
RUN cmake -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_INSTALL_PREFIX=/usr/local \
	  -DCMAKE_C_COMPILER=icc \
	  -DCMAKE_Fortran_COMPILER=ifort \
	  -DWITH_SOCKETS=ON .. \
    && cmake --build . \
    && cmake --install .

#-------------------------------------------------------------------------------

# Stripping the full installation into a fresh ubuntu image
FROM ubuntu:18.04

RUN apt-get -qq update -y\
	&& apt-get install -y -qq \
	python3 \
	python3-pip \
	strace \
	less

# Need to redefine correct version
ENV INTEL_VERSION="2022.1.0"

# Environment variables for MKL
ENV LD_LIBRARY_PATH="/opt/intel/oneapi/compiler/${INTEL_VERSION}/linux/compiler/lib/intel64_lin:/opt/intel/oneapi/mkl/${INTEL_VERSION}/lib/intel64"

# Copy numpy, scipy, Atomistica and DFTB+
COPY --from=builder /usr/local/ /usr/local/

# Copy Fortran runtime - but only the libraries we actually need
COPY --from=builder /opt/intel/oneapi/compiler/${INTEL_VERSION}/linux/compiler/lib/intel64_lin/*.so* /opt/intel/oneapi/compiler/${INTEL_VERSION}/linux/compiler/lib/intel64_lin/

# Copy MKL
COPY --from=builder /opt/intel/oneapi/mkl/${INTEL_VERSION}/lib/intel64/*.so* /opt/intel/oneapi/mkl/${INTEL_VERSION}/lib/intel64/

CMD ["python3"]
