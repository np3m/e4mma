# # Stage 1: Build stage
FROM awsteiner/o2scl:v0.930a4_u24.04_py AS builder
MAINTAINER Andrew W. Steiner (awsteiner0@protonmail.com)
# --------------------------------------------------------------------------
# Additional packages

RUN apt-get update -y && apt-get install -y libopenmpi-dev xclip sudo && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# --------------------------------------------------------------------------
# UTK EOS
# --------------------------------------------------------------------------
# Set makefile variables
ENV CXX=g++
ENV MPI_CXX=mpicxx
ENV CFLAGS="-I/usr/lib/x86_64-linux-gnu/openmpi/include/ \
  -I/usr/include/python3.12 -DO2SCL_PYTHON \
  -I/usr/local/lib/python3.12/dist-packages/numpy/core/include \
  -I/usr/local/hdf5/include -DO2SCL_OPENMP \
  -DO2SCL_PLAIN_HDF5_HEADER -DO2SCL_HDF5_COMP \
  -I/usr/local/include -std=c++14"
ENV LDFLAGS="-L/usr/local/hdf5/lib -lpython3.12"
#ENV LD_LIBRARY_PATH="/usr/local/hdf5/lib"

RUN git clone https://github.com/np3m/e4mma && \
    cd e4mma && \
    git checkout muses 
 
RUN cd e4mma/src && make -j 4 eos_nuclei eos
# --------------------------------------------------------------------------
# Stage 2: Dependencies MUSES
FROM python:3.11-slim AS deps
# Install git for pip install
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -yq \
        git \
    && rm -rf /var/lib/apt/lists/*
## Install dependencies to /usr/local/lib/python3.11/site-packages/
COPY requirements.txt ./
RUN pip --no-cache-dir install -r requirements.txt --no-compile
# --------------------------------------------------------------------------
# Stage 3: Runtime environment
FROM python:3.11-slim
# Install curl to download eos tables
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -yq \
        curl \
    && rm -rf /var/lib/apt/lists/*

ARG USERNAME=e4mma
ARG UID=1000

## Create user with desired UID and set "home" to /opt
RUN useradd --uid $UID --no-log-init --no-create-home --home-dir /opt --shell /bin/bash $USERNAME
RUN chown -R $UID:$UID /opt
USER $USERNAME

## Install Python dependencies
COPY --from=deps --chown=$UID:$UID /usr/local/lib/python3.11/site-packages/ /usr/local/lib/python3.11/site-packages/
# --------------------------------------------------------------------------
# Copy built files from the previous stage
COPY --from=builder --chown=$UID:$UID /usr/bin/make \
    /usr/bin/mpi* \
    /usr/bin/orte* \
    /usr/bin/rsh* \
    /usr/bin/ssh* \
    /usr/bin/
# --------------------------------------------------------------------------
COPY --from=builder --chown=$UID:$UID /usr/local/lib/libo2scl.so.0 /usr/local/lib/
COPY --from=builder --chown=$UID:$UID /usr/local/hdf5/lib/libhdf5.so.310 \
    /usr/local/hdf5/lib/libhdf5_hl.so.310 /usr/local/hdf5/lib/
COPY --from=builder --chown=$UID:$UID /usr/lib64/ld-linux-x86-64.so.2 /usr/lib64/
COPY --from=builder --chown=$UID:$UID \
    /lib/x86_64-linux-gnu/libgsl.so.27 \
    /lib/x86_64-linux-gnu/libpython3.12.so.1.0 \
    /lib/x86_64-linux-gnu/libmpi_cxx.so.40 \
    /lib/x86_64-linux-gnu/libmpi.so.40 \
    /lib/x86_64-linux-gnu/libstdc++.so.6 \
    /lib/x86_64-linux-gnu/libm.so.6 \
    /lib/x86_64-linux-gnu/libgomp.so.1 \
    /lib/x86_64-linux-gnu/libgcc_s.so.1 \
    /lib/x86_64-linux-gnu/libc.so.6 \
    /lib/x86_64-linux-gnu/libfftw3.so.3 \
    /lib/x86_64-linux-gnu/libquadmath.so.0 \
    /lib/x86_64-linux-gnu/libz.so.1 \
    /lib/x86_64-linux-gnu/libgslcblas.so.0 \
    /lib/x86_64-linux-gnu/libexpat.so.1 \
    /lib/x86_64-linux-gnu/libopen-pal.so.40 \
    /lib/x86_64-linux-gnu/libopen-rte.so.40 \
    /lib/x86_64-linux-gnu/libhwloc.so.15 \
    /lib/x86_64-linux-gnu/libevent_core-2.1.so.7 \
    /lib/x86_64-linux-gnu/libevent_pthreads-2.1.so.7 \
    /lib/x86_64-linux-gnu/libudev.so.1 \
    /lib/x86_64-linux-gnu/libcap.so.2 \
    #-----------------------
    /lib/x86_64-linux-gnu/libibverbs.so.1 \
    /lib/x86_64-linux-gnu/libucp.so.0 \
    /lib/x86_64-linux-gnu/libucs.so.0 \
    /lib/x86_64-linux-gnu/libuct.so.0 \
    /lib/x86_64-linux-gnu/libfabric.so.1 \
    /lib/x86_64-linux-gnu/libmca_common_ofi.so.10 \
    /lib/x86_64-linux-gnu/libmca_common_sm.so.40 \
    /lib/x86_64-linux-gnu/libmca_common_monitoring.so.50 \
    /lib/x86_64-linux-gnu/libmca_common_ofi.so.10 \
    /lib/x86_64-linux-gnu/libmca_common_verbs.so.40 \
    /lib/x86_64-linux-gnu/libinfinipath.so.4 \
    /lib/x86_64-linux-gnu/libpsm2.so.2 \
    /lib/x86_64-linux-gnu/libpsm_infinipath.so.1 \
    /lib/x86_64-linux-gnu/libnuma.so.1 \
    /lib/x86_64-linux-gnu/libfabric.so.1 \
    /lib/x86_64-linux-gnu/librdmacm.so.1 \
    /lib/x86_64-linux-gnu/libmca_common_ucx.so.40 \
    /lib/x86_64-linux-gnu/libefa.so.1 \
    /lib/x86_64-linux-gnu/libatomic.so.1 \
    /lib/x86_64-linux-gnu/libucm.so.0 \
    /lib/x86_64-linux-gnu/libmunge.so.2 \
    /lib/x86_64-linux-gnu/libnl-3.so.200 \
    /lib/x86_64-linux-gnu/libnl-route-3.so.200 \
    /lib/x86_64-linux-gnu/libgomp.so.1 \
    /lib/x86_64-linux-gnu/libpmix.so.2 \
    /lib/x86_64-linux-gnu/libopen-rte.so.40 \
    /lib/x86_64-linux-gnu/
# --------------------------------------------------------------------------
COPY --from=builder --chown=$UID:$UID /usr/lib/x86_64-linux-gnu/openmpi /usr/lib/x86_64-linux-gnu/openmpi
COPY --from=builder --chown=$UID:$UID /usr/lib/x86_64-linux-gnu/pmix2 /usr/lib/x86_64-linux-gnu/pmix2
COPY --from=builder --chown=$UID:$UID /usr/share/openmpi /usr/share/openmpi

COPY --from=builder --chown=$UID:$UID /usr/local/share/o2scl /usr/local/share/o2scl
# --------------------------------------------------------------------------
# Copy utk-eos files to the runtime stage
COPY --from=builder --chown=$UID:$UID e4mma/data /opt/e4mma/data
COPY --from=builder --chown=$UID:$UID e4mma/src/eos_nuclei e4mma/src/eos /opt/e4mma/src/ 
COPY --chown=$UID:$UID src/makefile \
    src/yaml_validator.py \
    src/yaml_generator.py \
    src/postprocess.py \
    src/Status.py \
    /opt/e4mma/src/

COPY --chown=$UID:$UID manifest.yaml /opt/e4mma/
COPY --chown=$UID:$UID test /opt/e4mma/test
COPY --chown=$UID:$UID api /opt/e4mma/api
## Create input and output directories
RUN mkdir /opt/e4mma/input/ && chown $UID:$UID /opt/e4mma/input/
RUN mkdir /opt/e4mma/output/ && chown $UID:$UID /opt/e4mma/output/
# --------------------------------------------------------------------------
# Set environment variables
ENV O2SCL_ADDL_LIBS="/usr/lib/gcc/x86_64-linux-gnu/12/libgomp.so"
ENV LD_LIBRARY_PATH="/usr/local/lib:/usr/local/hdf5/lib" 
# --------------------------------------------------------------------------
# Set working directory
WORKDIR /opt/e4mma/src/