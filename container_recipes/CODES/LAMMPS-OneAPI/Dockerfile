FROM ubuntu:20.04 AS builder

ARG DEBIAN_FRONTEND=noninteractive
ARG APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE=1

RUN apt-get update -y && apt-get install -y apt-utils

RUN apt-get install -y --no-install-recommends \
    build-essential \
    pkg-config \
    cmake \
    ca-certificates \
    gnupg

ARG url=https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.PUB
ADD $url /
RUN file=$(basename "$url") && \
    apt-key add "$file" && \
    rm "$file"

ARG apt_repo=https://apt.repos.intel.com/oneapi
RUN echo "deb $apt_repo all main" > /etc/apt/sources.list.d/oneAPI.list

# Attention: There is a second definition in the second stage
# This is the OneAPI version on NEMO. Load it with:
#   source /opt/bwhpc/intel/oneapi/setvars.sh
# (It is not available via module)
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
	curl \
	git \
	cmake \
	gfortran \
	python3 \
	python3-dev \
	python3-pip \
	pybind11-dev \
	intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic-${ONEAPI_VERSION} \
	intel-oneapi-compiler-fortran-${ONEAPI_VERSION} \
	intel-oneapi-mkl-devel-${ONEAPI_VERSION} \
	intel-oneapi-mpi-devel-${ONEAPI_VERSION}

ENV PATH=/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/libfabric/bin:/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/bin:/opt/intel/oneapi/mkl/${ONEAPI_VERSION}/bin/intel64:/opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/bin/intel64:/opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
ENV CPATH=/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/include:/opt/intel/oneapi/mkl/${ONEAPI_VERSION}/include:/opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/include
ENV LD_LIBRARY_PATH=/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/libfabric/lib:/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/lib/release:/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/lib:/opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64:/opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/lib:/opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/lib/x64:/opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/compiler/lib/intel64_lin
ENV LIBRARY_PATH=/opt/intel/oneapi/mpi/${ONEAPI_VERSION}//libfabric/lib:/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/lib/release:/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/lib:/opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64:/opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/compiler/lib/intel64_lin:/opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/lib
ENV MKLROOT=/opt/intel/oneapi/mkl/${ONEAPI_VERSION}
ENV I_MPI_ROOT=/opt/intel/oneapi/mpi/${ONEAPI_VERSION}

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
# Install parallel version of the NetCDF library from the sources.
# This is necessary because parallel compiles (if existing) are
# broken on most distributions.
#
ENV PNETCDF_VERSION=1.12.2
WORKDIR /tmp
RUN apt-get install -y m4 \
	&& curl https://parallel-netcdf.github.io/Release/pnetcdf-${PNETCDF_VERSION}.tar.gz | tar -xz
WORKDIR /tmp/pnetcdf-${PNETCDF_VERSION}
RUN MPICC=mpiicc MPICXX=mpiicpc ./configure --disable-fortran --disable-cxx --enable-shared --prefix=/usr/local \
	&& make -j 4 \
	&& make install

#
# Install parallel version of HDF5
#
ENV HDF5_VERSION=1.10.6 
WORKDIR /tmp
RUN apt-get install -y libcurl4-openssl-dev \
	&& curl https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-${HDF5_VERSION}/src/hdf5-${HDF5_VERSION}.tar.gz | tar -xz
WORKDIR /tmp/hdf5-${HDF5_VERSION}
RUN MPICC=mpiicc MPICXX=mpiicpc ./configure --enable-parallel --prefix=/usr/local \
	&& make -j 4 \
	&& make install

#
# We need to compile NetCDF ourselves because there is no package that has
# parallel PnetCDF and HDF5 enabled.
#
ENV NETCDF4_VERSION=4.9.0
WORKDIR /tmp
RUN curl https://downloads.unidata.ucar.edu/netcdf-c/${NETCDF4_VERSION}/netcdf-c-${NETCDF4_VERSION}.tar.gz | tar -xz && mkdir /tmp/netcdf-c-build
WORKDIR /tmp/netcdf-c-build
RUN cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local -DMPI_CC_COMPILER=mpiicc -DMPI_CXX_COMPILER=mpiicpc -DUSE_PARALLEL=ON -DENABLE_PARALLEL4=ON -DENABLE_PNETCDF=ON -DENABLE_TESTS=OFF /tmp/netcdf-c-${NETCDF4_VERSION} \
	&& make -j 4 \
	&& make install

ENV LAMMPS_VERSION="stable_23Jun2022_update1"
WORKDIR /tmp
RUN apt-get install -y libjpeg62-dev libpng-dev \
	&& curl -L https://github.com/lammps/lammps/archive/refs/tags/${LAMMPS_VERSION}.tar.gz | tar xz
WORKDIR /tmp/lammps-${LAMMPS_VERSION}
RUN mkdir build \
	&& cd build \
	&& cmake \
	-DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_INSTALL_PREFIX=/usr/local \
	-DBUILD_MPI=YES \
	-DMPI_CXX_COMPILER=mpiicpc \
	-DBUILD_OMP=NO \
	-DBUILD_SHARED_LIBS=YES \
	-DFFT=MKL \
	-DFFT_MKL_THREADS=OFF \
	-DWITH_JPEG=YES \
	-DWITH_PNG=YES \
	-DPKG_MANYBODY=ON \
	-DPKG_KSPACE=ON \
	-DPKG_MOLECULE=ON \
	-DPKG_MPIIO=ON \
	-DPKG_NETCDF=ON \
	-DPKG_H5MD=ON \
	-DPKG_MOLFILE=ON \
	-DPKG_OPT=ON \
	-DPKG_PYTHON=ON \
	-DPKG_QEQ=ON \
	-DPKG_RIGID=ON \
	-DPKG_UEF=ON \
	-DPKG_ML-PACE=ON \
	../cmake \
	&& cmake --build . \
	&& make install

#
# Install simple MPI benchmark
#

WORKDIR /tmp
RUN git clone https://github.com/LLNL/mpiBench.git \
	&& cd mpiBench \
	&& make \
	&& cp mpiBench /usr/local/bin

FROM ubuntu:20.04

RUN apt-get update -y \
	&& apt-get install -y \
 	python3 \
 	python3-pip \
 	strace \
 	less \
 	libcurl4 \
 	libjpeg62 \
 	libpng16-16

ENV ONEAPI_VERSION="2021.3.0"

# Environment variables for MKL
ENV PATH=/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
ENV LD_LIBRARY_PATH=/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/libfabric/lib:/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/lib/release:/opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/compiler/lib/intel64_lin:/opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64:/usr/local/lib

# Copy everything we build ourselves
COPY --from=builder /usr/local/ /usr/local/

# Copy compiler runtime - but only the libraries we actually need
COPY --from=builder /opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/compiler/lib/intel64_lin/*.so* /opt/intel/oneapi/compiler/${ONEAPI_VERSION}/linux/compiler/lib/intel64_lin/

# Copy MPI - but only the libraries we actually need
COPY --from=builder /opt/intel/oneapi/mpi/${ONEAPI_VERSION}/bin /opt/intel/oneapi/mpi/${ONEAPI_VERSION}/bin/
COPY --from=builder /opt/intel/oneapi/mpi/${ONEAPI_VERSION}/lib /opt/intel/oneapi/mpi/${ONEAPI_VERSION}/lib/
COPY --from=builder /opt/intel/oneapi/mpi/${ONEAPI_VERSION}/libfabric /opt/intel/oneapi/mpi/${ONEAPI_VERSION}/libfabric/
COPY --from=builder /opt/intel/oneapi/mpi/${ONEAPI_VERSION}/etc /opt/intel/oneapi/mpi/${ONEAPI_VERSION}/etc/

# Copy MKL
COPY --from=builder /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/libmkl_rt.so* /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/
COPY --from=builder /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/libmkl_core.so* /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/
COPY --from=builder /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/libmkl_intel_thread.so* /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/
COPY --from=builder /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/libmkl_intel_lp64.so* /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/
COPY --from=builder /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/libmkl_def.so* /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/
COPY --from=builder /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/libmkl_avx*.so* /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/
# COPY --from=builder /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/libmemkind.so* /opt/intel/oneapi/mkl/${ONEAPI_VERSION}/lib/intel64/

# Just make sure there is no OpenMP madness goingon
ENV OMP_NUM_THREADS=1
ENV MKL_NUM_THREADS=1

# Fix FI_PROVIDER_PATH
ENV FI_PROVIDER_PATH=/opt/intel/oneapi/mpi/${ONEAPI_VERSION}/libfabric/lib/prov

CMD ["/usr/local/bin/lmp"]
