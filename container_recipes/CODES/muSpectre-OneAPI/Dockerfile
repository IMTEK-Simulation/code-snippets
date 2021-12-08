#
# This is a Docker file to build muSpectre with FFTW, PFFT and PNetCDF using Intel oneAPI
#


FROM intel/oneapi:os-tools-ubuntu18.04 AS builder

# this is the intelmpi version of Galileo100
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
	libcurl4-openssl-dev \
	m4 \
	libtool \
	automake \
	python3 \
	python3-dev \
	python3-pip \
	libopenblas-dev \
	intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic-${ONEAPI_VERSION} \
	intel-oneapi-compiler-fortran-${ONEAPI_VERSION} \
	intel-oneapi-mpi-devel-${ONEAPI_VERSION}

ENV PATH="/opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/bin/intel64:/opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/bin:/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
ENV CPATH="/opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/include:/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/include"
ENV LD_LIBRARY_PATH="/opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/lib:/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/lib/release:/opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/lib/x64:/opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/compiler/lib/intel64_lin:/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/lib:/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/libfabric/lib"
ENV LIBRARY_PATH="/opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/compiler/lib/intel64_lin:/opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/lib:/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/lib/release:/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/lib:/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/libfabric/lib"

#
# Install simple MPI benchmark
#

WORKDIR /tmp
RUN git clone https://github.com/LLNL/mpiBench.git \
	&& cd mpiBench \
	&& make \
	&& cp mpiBench /usr/local/bin

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
# Manual install of MPI-parallel FFTW
#
ENV FFTW_VERSION=3.3.9
WORKDIR /tmp
RUN curl -L http://www.fftw.org/fftw-${FFTW_VERSION}.tar.gz | tar -xz
WORKDIR /tmp/fftw-${FFTW_VERSION}
RUN ./configure --enable-shared --enable-mpi --disable-fortran --enable-sse2 --enable-avx --enable-avx2 --prefix=/usr/local \
	&& make -j 4 \
	&& make install

#
# Install current master of PFFT
#
WORKDIR /tmp
RUN git clone https://github.com/mpip/pfft.git pfft
WORKDIR /tmp/pfft
RUN ./bootstrap.sh \
	&& CFLAGS="-O3 -march=sandybridge -mno-avx512f" ./configure --disable-fortran --enable-shared --prefix=/usr/local \
	&& make -j 4 \
	&& make install 

#
# Install parallel version of the NetCDF library from the sources.
# This is necessary because parallel compiles (if existing) are
# broken on most distributions.
#
ENV PNETCDF_VERSION=1.12.2
WORKDIR /tmp
RUN curl https://parallel-netcdf.github.io/Release/pnetcdf-${PNETCDF_VERSION}.tar.gz | tar -xz
WORKDIR /tmp/pnetcdf-${PNETCDF_VERSION}
RUN ./configure --disable-fortran --disable-cxx --enable-shared --prefix=/usr/local \
	&& make -j 4 \
	&& make install

#
# Install parallel version of HDF5
#
ENV HDF5_VERSION=1.10.6 
WORKDIR /tmp
RUN curl https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-${HDF5_VERSION}/src/hdf5-${HDF5_VERSION}.tar.gz | tar -xz
WORKDIR /tmp/hdf5-${HDF5_VERSION}
RUN ./configure --enable-parallel --prefix=/usr/local \
	&& make -j 4 \
	&& make install

#
# We need to compile NetCDF ourselves because there is no package that has
# parallel PnetCDF and HDF5 enabled.
#
ENV NETCDF4_VERSION=4.7.3
WORKDIR /tmp
RUN curl https://www.unidata.ucar.edu/downloads/netcdf/ftp/netcdf-c-${NETCDF4_VERSION}.tar.gz | tar -xz && mkdir /tmp/netcdf-c-build
WORKDIR /tmp/netcdf-c-build
RUN cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local -DUSE_PARALLEL=ON -DENABLE_PARALLEL4=ON -DENABLE_PNETCDF=ON -DENABLE_TESTS=OFF /tmp/netcdf-c-${NETCDF4_VERSION} \
	&& make -j 4 \
	&& make install

#
# Install netcdf4-python and make sure that it is compiled (no-binary),
# otherwise it will not have parallel support.
#
ENV NETCDF4_PYTHON_VERSION=1.5.5.1
RUN NETCDF4_DIR=/usr/local HDF5_DIR=/usr/local python3 -m pip install --no-binary netCDF4 netCDF4==${NETCDF4_PYTHON_VERSION}

ENV MPI4PY_VERSION=3.0.3 
RUN python3 -m pip install --no-binary mpi4py mpi4py==${MPI4PY_VERSION}

#
# Get muSpectre branch, compile and install
#
RUN apt-get install -y zlib1g-dev libjpeg-dev libboost-test-dev
RUN python3 -m pip install NuMPI meshio[all] matplotlib uvw
ENV MUSPECTRE_BRANCH=master
WORKDIR /tmp
RUN git clone -b ${MUSPECTRE_BRANCH} https://gitlab.com/muspectre/muspectre.git
# If you want to try compiling with the Intel compilers: Remove the -Weffc++
# option as this messes with the Intel Compiler (It complains that effective
# C++ does not allow C-style comments, which are used in muSpectre.)
# --- && sed -i 's/-Weffc++//g' muspectre/cmake/muspectreTools.cmake
WORKDIR /tmp/muspectre/build
# NEMO is Sandy Bridge, JUWELS is Skylake
# NEMO: -march=sandybridge -mno-avx512f
# export OPT_FLAGS="-O2 -march=sandybridge -mno-avx512f"
# JUWELS: -march=skylake
# export OPT_FLAGS="-O2 -march=skylake"
# Galileo100 is CPU Intel CascadeLake
#ENV OPT_FLAGS="-O2 -march=cascadelake"
ENV OPT_FLAGS="-O2"
RUN cmake -DCMAKE_C_COMPILER=mpigcc -DCMAKE_CXX_COMPILER=mpigxx -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local -DCMAKE_C_FLAGS_RELEASE="${OPT_FLAGS}" -DCMAKE_CXX_FLAGS_RELEASE="${OPT_FLAGS}" -DMUSPECTRE_MPI_PARALLEL=ON -DMUFFT_USE_FFTWMPI=ON -DMUFFT_USE_PFFT=ON -DNETCDF_IO=ON -DMUSPECTRE_MAKE_TESTS=OFF .. \
	&& make -j 4 \
	&& make install

ENV LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/usr/local/lib"

#
# Copy muSpectre to smaller container
#
FROM ubuntu:18.04

ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

ENV ONEAPI_VERSION="2021.3.0"

RUN apt-get update -y \
	&& apt-get install -y \
	python3 \
	python3-pip \
	strace \
	less \
	curl \
	libopenblas-base

# Environment variables for MPI and MKL
ENV PATH="/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/libfabric/bin:/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/bin:/usr/local/bin:$PATH"
ENV LD_LIBRARY_PATH="/opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/compiler/lib/intel64_lin:/opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64:/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/libfabric/lib:/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/lib:/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/lib/release:/usr/local/lib:$LD_LIBRARY_PATH"
ENV FI_PROVIDER_PATH=/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/libfabric/lib/prov

# Copy everything we installed above
COPY --from=builder /usr/local/ /usr/local/

# Copy compiler runtime - but only the libraries we actually need
COPY --from=builder /opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/compiler/lib/intel64_lin/*.so* /opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/compiler/lib/intel64_lin/

# Copy MKL - but only the libaries we actually need
#COPY --from=builder /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/libmkl_sequential.so* /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/
#COPY --from=builder /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/libmkl_rt.so* /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/
#COPY --from=builder /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/libmkl_core.so* /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/
#COPY --from=builder /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/libmkl_intel_thread.so* /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/
#COPY --from=builder /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/libmkl_intel_lp64.so* /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/
#COPY --from=builder /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/libmkl_def.so* /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/
#COPY --from=builder /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/libmkl_avx*.so* /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/

# Copy MPI - but only the libraries we actually need
COPY --from=builder /opt/intel/oneapi/mpi/${ONEAPI_VERSION}/bin /opt/intel/oneapi/mpi/${ONEAPI_VERSION}/bin/
COPY --from=builder /opt/intel/oneapi/mpi/${ONEAPI_VERSION}/lib /opt/intel/oneapi/mpi/${ONEAPI_VERSION}/lib/
COPY --from=builder /opt/intel/oneapi/mpi/${ONEAPI_VERSION}/libfabric /opt/intel/oneapi/mpi/${ONEAPI_VERSION}/libfabric/
COPY --from=builder /opt/intel/oneapi/mpi/${ONEAPI_VERSION}/etc /opt/intel/oneapi/mpi/${ONEAPI_VERSION}/etc/

# Just make sure there is no OpenMP madness goingon
ENV OMP_NUM_THREADS=1
ENV MKL_NUM_THREADS=1

CMD ["python3"]
