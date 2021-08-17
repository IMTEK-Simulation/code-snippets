#
# This is the Docker image used for building this...
#

FROM intel/oneapi:os-tools-ubuntu18.04 AS builder

ENV ONEAPI_VERSION="2021.3.0"

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
	intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic-${ONEAPI_VERSION} \
	intel-oneapi-compiler-fortran-${ONEAPI_VERSION} \
	intel-oneapi-mkl-devel-${ONEAPI_VERSION} \
	intel-oneapi-mpi-devel-${ONEAPI_VERSION}

ENV PATH="/opt/intel/oneapi/mkl/${ONEAPI_VERSION}/bin/intel64:/opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/bin/intel64:/opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/bin:/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
ENV CPATH="/opt/intel/oneapi/mkl/${ONEAPI_VERSION}/include:/opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/include:/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/include"
ENV LD_LIBRARY_PATH="/opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64:/opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/lib:/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/lib/release:/opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/lib/x64:/opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/compiler/lib/intel64_lin:/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/lib:/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/libfabric/lib"
ENV LIBRARY_PATH="/opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64:/opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/compiler/lib/intel64_lin:/opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/lib:/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/lib/release:/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/lib:/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/libfabric/lib"
ENV MKLROOT="/opt/intel/oneapi/mkl/${ONEAPI_VERSION}"

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
# ASE 3.22.0 belongs to GPAW 21.6.0 - see https://wiki.fysik.dtu.dk/gpaw/releasenotes.html
ENV ASE_VERSION=3.22.0
RUN python3 -m pip install --upgrade pip \
	&& python3 -m pip install ase==${ASE_VERSION} \
	&& python3 -m pip install mpi4py matscipy
ENV GPAW_VERSION=21.6.0
WORKDIR /tmp
RUN echo 'mpicompiler="mpiicc"' > gpaw_siteconfig.py \
	&& echo 'libraries =["xc", "mkl_rt", "pthread", "m"]' >> gpaw_siteconfig.py \
	&& LANG=C.UTF-8 GPAW_CONFIG=/tmp/gpaw_siteconfig.py python3 -m pip -v install gpaw==${GPAW_VERSION}

# Install simple MPI benchmark
WORKDIR /tmp
RUN git clone https://github.com/LLNL/mpiBench.git \
	&& cd mpiBench \
	&& make \
	&& cp mpiBench /usr/local/bin

#
# ...which we then copy to the Docker image for running things
#

FROM ubuntu:18.04

ENV ONEAPI_VERSION="2021.3.0"

RUN apt-get update -y \
	&& apt-get install -y \
	python3 \
	python3-pip \
	strace \
	less \
	curl

# Environment variables for MPI and MKL
ENV PATH="/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/libfabric/bin:/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/bin:/usr/local/bin:$PATH"
ENV LD_LIBRARY_PATH="/opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/compiler/lib/intel64_lin:/opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64:/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/libfabric/lib:/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/lib:/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/lib/release:/usr/local/lib:$LD_LIBRARY_PATH"
ENV FI_PROVIDER_PATH=/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/libfabric/lib/prov

# Copy numpy, scipy, Atomistica and DFTB+
COPY --from=builder /usr/local/ /usr/local/

# Copy Fortran runtime - but only the libraries we actually need
COPY --from=builder /opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/compiler/lib/intel64_lin/*.so* /opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/compiler/lib/intel64_lin/

# Copy MKL - but only the libaries we actually need
COPY --from=builder /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/libgcc_s.so* /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/
COPY --from=builder /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/libmkl_sequential.so* /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/
COPY --from=builder /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/libmkl_rt.so* /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/
COPY --from=builder /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/libmkl_core.so* /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/
COPY --from=builder /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/libmkl_intel_thread.so* /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/
COPY --from=builder /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/libmkl_intel_lp64.so* /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/
COPY --from=builder /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/libmkl_def.so* /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/
COPY --from=builder /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/libmkl_avx*.so* /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/
COPY --from=builder /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/libmemkind.so* /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/
COPY --from=builder /opt/intel/oneapi/mpi/${ONEAPI_VERSION}/bin /opt/intel/oneapi/mpi/${ONEAPI_VERSION}/bin/
COPY --from=builder /opt/intel/oneapi/mpi/${ONEAPI_VERSION}/lib /opt/intel/oneapi/mpi/${ONEAPI_VERSION}/lib/
COPY --from=builder /opt/intel/oneapi/mpi/${ONEAPI_VERSION}/libfabric /opt/intel/oneapi/mpi/${ONEAPI_VERSION}/libfabric/
COPY --from=builder /opt/intel/oneapi/mpi/${ONEAPI_VERSION}/etc /opt/intel/oneapi/mpi/${ONEAPI_VERSION}/etc/

# Get GPAW setups
ENV GPAW_SETUPS_VERSION="0.9.20000"
RUN mkdir -p /usr/local/share \
	&& curl -L https://wiki.fysik.dtu.dk/gpaw-files/gpaw-setups-${GPAW_SETUPS_VERSION}.tar.gz | tar -xzC /usr/local/share

# Set environment variable for GPAW setups
ENV GPAW_SETUP_PATH=/usr/local/share/gpaw-setups-${GPAW_SETUPS_VERSION}

# Just make sure there is no OpenMP madness goingon
ENV OMP_NUM_THREADS=1
ENV MKL_NUM_THREADS=1

CMD ["python3"]
