#
# This is the Docker image used for building this...
#

FROM intel/oneapi:os-tools-ubuntu18.04 AS builder

ARG DEBIAN_FRONTEND=noninteractive

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
	intel-oneapi-mpi-devel

ENV PATH='/opt/intel/oneapi/mkl/2021.3.0/bin/intel64:/opt/intel/oneapi/compiler/2021.3.0/linux/bin/intel64:/opt/intel/oneapi/compiler/2021.3.0/linux/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin'
ENV CPATH='/opt/intel/oneapi/mkl/2021.3.0/include:/opt/intel/oneapi/compiler/2021.3.0/linux/include'
ENV LD_LIBRARY_PATH='/opt/intel/oneapi/mkl/2021.3.0/lib/intel64:/opt/intel/oneapi/compiler/2021.3.0/linux/lib:/opt/intel/oneapi/compiler/2021.3.0/linux/lib/x64:/opt/intel/oneapi/compiler/2021.3.0/linux/compiler/lib/intel64_lin'
ENV LIBRARY_PATH='/opt/intel/oneapi/mkl/2021.3.0/lib/intel64:/opt/intel/oneapi/compiler/2021.3.0/linux/compiler/lib/intel64_lin:/opt/intel/oneapi/compiler/2021.3.0/linux/lib'
ENV MKLROOT='/opt/intel/oneapi/mkl/2021.3.0'

#
# Numpy 1.19.5 and scipy 1.5.4 are last versions that support Ubuntu 18.04's Python 3.6
#

ENV NUMPY_VERSION="1.19.5"
WORKDIR /tmp
RUN python3 -m pip install "setuptools<49.2.0" "cython==0.29.21" \
	&& curl -L https://github.com/numpy/numpy/archive/refs/tags/v${NUMPY_VERSION}.tar.gz | tar xz
WORKDIR /tmp/numpy-${NUMPY_VERSION}
RUN python3 setup.py config --compiler=intelem build_clib --compiler=intelem build_ext --compiler=intelem install

ENV SCIPY_VERSION="1.5.4"
WORKDIR /tmp
RUN curl -L https://github.com/scipy/scipy/archive/refs/tags/v${SCIPY_VERSION}.tar.gz | tar xz
WORKDIR /tmp/scipy-${SCIPY_VERSION}
RUN python3 setup.py config --compiler=intelem --fcompiler=intelem build_clib --compiler=intelem --fcompiler=intelem build_ext --compiler=intelem --fcompiler=intelem install

#
# Install libxc; GPAW currently does not support libxc 5
#

ENV LIBXC_VERSION="4.3.4"
WORKDIR /tmp
RUN curl -L http://www.tddft.org/programs/libxc/down.php?file=${LIBXC_VERSION}/libxc-${LIBXC_VERSION}.tar.gz | tar xz
WORKDIR /tmp/libxc-${LIBXC_VERSION}
RUN CC=icc CXX=icpc ./configure --prefix=/usr/local --enable-shared \
    && make \
    && make check \
    && make install

# Install GPAW - Python module and standalone code
RUN python3 -m pip install --upgrade pip \
    && python3 -m pip install gpaw

#
# ...which we then copy to the Docker image for running things
#

FROM ubuntu:18.04

RUN apt-get update -y \
	&& apt-get install -y \
	python3 \
	python3-pip \
	strace \
	less

# Environment variables for MPI and MKL
ENV PATH="/opt/intel/oneapi/mpi/2021.3.0/libfabric/bin:/opt/intel/oneapi/mpi/2021.3.0/bin:/usr/local/bin:$PATH"
ENV LD_LIBRARY_PATH="/opt/intel/oneapi/compiler/2021.3.0/linux/compiler/lib/intel64_lin:/opt/intel/oneapi/mkl/2021.3.0/lib/intel64:/opt/intel/oneapi/mpi/2021.3.0/libfabric/lib:/opt/intel/oneapi/mpi/2021.3.0/lib:/usr/local/lib:$LD_LIBRARY_PATH"

# Copy numpy, scipy, Atomistica and DFTB+
COPY --from=builder /usr/local/ /usr/local/

# Copy Fortran runtime - but only the libraries we actually need
COPY --from=builder /opt/intel/oneapi/compiler/2021.3.0/linux/compiler/lib/intel64_lin/*.so* /opt/intel/oneapi/compiler/2021.3.0/linux/compiler/lib/intel64_lin/

# Copy MKL - but only the libaries we actually need
COPY --from=builder /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/libgcc_s.so* /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/
COPY --from=builder /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/libmkl_sequential.so* /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/
COPY --from=builder /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/libmkl_rt.so* /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/
COPY --from=builder /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/libmkl_core.so* /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/
COPY --from=builder /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/libmkl_intel_thread.so* /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/
COPY --from=builder /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/libmkl_intel_lp64.so* /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/
COPY --from=builder /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/libmkl_def.so* /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/
COPY --from=builder /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/libmkl_avx*.so* /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/
COPY --from=builder /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/libmemkind.so* /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/
COPY --from=builder /opt/intel/oneapi/mpi/2021.3.0/bin /opt/intel/oneapi/mpi/2021.3.0/
COPY --from=builder /opt/intel/oneapi/mpi/2021.3.0/lib /opt/intel/oneapi/mpi/2021.3.0/
COPY --from=builder /opt/intel/oneapi/mpi/2021.3.0/etc /opt/intel/oneapi/mpi/2021.3.0/

# ASE 3.18.0 belongs to GPAW 19.8.1
RUN pip install --upgrade pip \
	&& pip install ase==3.18.0

CMD ["python3"]
