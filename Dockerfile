# Stage 1: Build dependencies
FROM ubuntu:23.04 AS builder

# Set up build dependencies
RUN apt-get -y update && apt-get \
  -y install g++ make libgsl-dev autoconf automake libtool git \
  libhdf5-dev libreadline-dev libboost-all-dev libeigen3-dev cmake \
  libopenblas-dev liblapack-dev libarpack2-dev libsuperlu-dev \
  libarmadillo-dev libfftw3-dev curl texlive dvipng texlive-latex-extra \
  cm-super libcairo2-dev libyaml-cpp-dev

WORKDIR /opt

# Clone and install o2scl
RUN git clone https://github.com/awsteiner/o2scl \
    && cd o2scl \
    && git checkout eee9fd83 && autoreconf -i \
# We disable static and make blank-doc to keep the image small
    && LDFLAGS="-larmadillo -llapack -lblas" CXXFLAGS="-O3 -DO2SCL_UBUNTU_HDF5 -DO2SCL_HDF5_PRE_1_12 -DO2SCL_REGEX -DO2SCL_HDF5_COMP -DO2SCL_NO_BOOST_MULTIPRECISION -I/usr/include" ./configure --enable-eigen --enable-armadillo --enable-openmp --enable-fftw --disable-static \
    && make blank-doc && make -j 3 && make install \ 
    && cd ../

RUN git clone https://github.com/np3m/e4mma && \
    cd e4mma && \
    git checkout muses 

COPY src/eos_nuclei.h /opt/e4mma/src/
COPY src/eos_nuclei.cpp /opt/e4mma/src/
RUN  cd e4mma/src/ && make -j 4 eos_nuclei

FROM python:3.11-slim AS deps
# Install git for pip install
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -yq \
        git \
    && rm -rf /var/lib/apt/lists/*
## Install dependencies to /usr/local/lib/python3.11/site-packages/
COPY requirements.txt ./
RUN pip --no-cache-dir install -r requirements.txt --no-compile

# Stage 2: Runtime environment
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

# Copy built files from the previous stage
COPY --from=builder --chown=$UID:$UID /usr/bin/make \
    /usr/bin/mpi* \
    /usr/bin/orte* \
    /usr/bin/rsh* \
    /usr/bin/ssh* \
    /usr/bin/
COPY --from=builder --chown=$UID:$UID /usr/lib/libarmadillo.so.* /usr/lib/
COPY --from=builder --chown=$UID:$UID /usr/lib/x86_64-linux-gnu/libfftw.so.* \
    /usr/lib/x86_64-linux-gnu/libhdf5_serial.so.103 \
    /usr/lib/x86_64-linux-gnu/libhdf5_serial_hl.so.100 \
    /usr/lib/x86_64-linux-gnu/libreadline.so.8 \
    /usr/lib/x86_64-linux-gnu/libgsl.so.27 \
    /usr/lib/x86_64-linux-gnu/libm.so.6 \
    /usr/lib/x86_64-linux-gnu/libblas.so.3 \
    /usr/lib/x86_64-linux-gnu/liblapack.so.3 \
    /usr/lib/x86_64-linux-gnu/libarpack.so.2 \
    /usr/lib/x86_64-linux-gnu/libsuperlu.so.5 \
    /usr/lib/x86_64-linux-gnu/libcurl.so.4 \
    /usr/lib/x86_64-linux-gnu/libsz.so.2 \
    /usr/lib/x86_64-linux-gnu/libgslcblas.so.0 \
    /usr/lib/x86_64-linux-gnu/libopenblas.so.0 \
    /usr/lib/x86_64-linux-gnu/libgfortran.so.5 \
    /usr/lib/x86_64-linux-gnu/libnghttp2.so.14 \
    /usr/lib/x86_64-linux-gnu/librtmp.so.1 \
    /usr/lib/x86_64-linux-gnu/libssh.so.4 \
    /usr/lib/x86_64-linux-gnu/libpsl.so.5 \
    /usr/lib/x86_64-linux-gnu/libldap.so.2 \
    /usr/lib/x86_64-linux-gnu/liblber.so.2 \
    /usr/lib/x86_64-linux-gnu/libbrotlidec.so.* \
    /usr/lib/x86_64-linux-gnu/libaec.so.0 \
    /usr/lib/x86_64-linux-gnu/libquadmath.so.0 \
    /usr/lib/x86_64-linux-gnu/libsasl2.so.2 \
    /usr/lib/x86_64-linux-gnu/libbrotlicommon.so.1 \
    /usr/lib/x86_64-linux-gnu/libpsm_infinipath.so.1 \
    /usr/lib/x86_64-linux-gnu/libgomp.so.1 \
    /usr/lib/x86_64-linux-gnu/libfftw3.so.3 \
    /usr/lib/x86_64-linux-gnu/libmpi_cxx.so.40 \
    /usr/lib/x86_64-linux-gnu/libmpi.so.40 \
    /usr/lib/x86_64-linux-gnu/libopen-pal.so.40 \
    /usr/lib/x86_64-linux-gnu/libopen-rte.so.40 \
    /usr/lib/x86_64-linux-gnu/libhwloc.so.15 \
    /usr/lib/x86_64-linux-gnu/libevent_core-2.1.so.7 \
    /usr/lib/x86_64-linux-gnu/libevent_pthreads-2.1.so.7 \
    /usr/lib/x86_64-linux-gnu/libnl-3.so.200 \
    /usr/lib/x86_64-linux-gnu/libpmix.so.2 \
    /usr/lib/x86_64-linux-gnu/libmunge.so.2 \
    #-----
    /usr/lib/x86_64-linux-gnu/libmca_common_ofi.so.10 \
    /usr/lib/x86_64-linux-gnu/libibverbs.so.1 \
    /usr/lib/x86_64-linux-gnu/libmca_common_monitoring.so.50 \
    /usr/lib/x86_64-linux-gnu/libucp.so.0 \
    /usr/lib/x86_64-linux-gnu/libpsm2.so.2 \
    /usr/lib/x86_64-linux-gnu/libinfinipath.so.4 \
    /usr/lib/x86_64-linux-gnu/libmca_common_ofi.so.10 \
    /usr/lib/x86_64-linux-gnu/libmca_common_sm.so.40 \
    /usr/lib/x86_64-linux-gnu/libfabric.so.1 \
    /usr/lib/x86_64-linux-gnu/libmca_common_verbs.so.40 \
    /usr/lib/x86_64-linux-gnu/libucs.so.0 \
    /usr/lib/x86_64-linux-gnu/libnuma.so.1 \
    /usr/lib/x86_64-linux-gnu/librdmacm.so.1 \
    /usr/lib/x86_64-linux-gnu/libnl-route-3.so.200 \
    /usr/lib/x86_64-linux-gnu/libmca_common_ucx.so.40 \
    /usr/lib/x86_64-linux-gnu/libefa.so.1 \
    /usr/lib/x86_64-linux-gnu/libuct.so.0 \
    /usr/lib/x86_64-linux-gnu/libatomic.so.1 \
    /usr/lib/x86_64-linux-gnu/libucm.so.0 \
    /usr/lib/x86_64-linux-gnu/

COPY --from=builder --chown=$UID:$UID /usr/lib/x86_64-linux-gnu/openmpi /usr/lib/x86_64-linux-gnu/openmpi
COPY --from=builder --chown=$UID:$UID /usr/lib/x86_64-linux-gnu/pmix2 /usr/lib/x86_64-linux-gnu/pmix2
COPY --from=builder --chown=$UID:$UID /usr/share/openmpi /usr/share/openmpi
COPY --from=builder --chown=$UID:$UID /usr/local/lib/libo2scl.so.* /usr/local/lib/

COPY --from=builder --chown=$UID:$UID /usr/local/share/o2scl /usr/local/share/o2scl

# Copy utk-eos files to the runtime stage
COPY --from=builder --chown=$UID:$UID /opt/e4mma/data /opt/e4mma/data
COPY --from=builder --chown=$UID:$UID /opt/e4mma/src/eos_nuclei /opt/e4mma/src/ 
COPY --chown=$UID:$UID src/makefile \
    src/yaml_validator.py \
    src/point_validator.py \
    src/yaml_generator.py \
    src/point_generator.py \
    src/postprocess.py \
    src/Status.py \
    /opt/e4mma/src/

COPY --chown=$UID:$UID manifest.yaml /opt/e4mma/
COPY --chown=$UID:$UID test /opt/e4mma/test
COPY --chown=$UID:$UID api /opt/e4mma/api
## Create input and output directories
RUN mkdir /opt/e4mma/input/ && chown $UID:$UID /opt/e4mma/input/
RUN mkdir /opt/e4mma/output/ && chown $UID:$UID /opt/e4mma/output/


# Set environment variables
ENV O2SCL_ADDL_LIBS="/usr/lib/gcc/x86_64-linux-gnu/12/libgomp.so"
ENV LD_LIBRARY_PATH="/usr/local/lib" 

# Set working directory
WORKDIR /opt/e4mma/test/
