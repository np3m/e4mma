FROM python:3.10 as o2scl

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -yq \
        g++ \
        make \
        libgsl-dev \
        autoconf \
        automake \
        libtool \
        git \
        libhdf5-dev \
        libncurses-dev \
        libreadline-dev \
        libboost-all-dev \
        libeigen3-dev \
        cmake \
        libopenblas-dev \
        liblapack-dev \
        libarpack2-dev \
        libsuperlu-dev \
        libarmadillo-dev \
        libfftw3-dev \
        curl \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt
RUN git clone https://github.com/awsteiner/o2scl && \
    cd o2scl && \
    git checkout e7b62ca63047ae183b209755a77e51380ecd8780

COPY requirements.o2scl.txt requirements.txt
RUN pip install -r requirements.txt

WORKDIR /opt/o2scl
RUN autoreconf -i
RUN export LDFLAGS="-larmadillo -llapack -lblas -lncurses" && \
    export CXXFLAGS="" && \
    export CXXFLAGS="${CXXFLAGS} -O3 -DO2SCL_UBUNTU_HDF5 -DO2SCL_HDF5_PRE_1_12 -DO2SCL_REGEX -DO2SCL_HDF5_COMP" && \
    export CXXFLAGS="${CXXFLAGS} -I/usr/include -I/usr/lib/x86_64-linux-gnu/hdf5/serial/include" && \
    export CXXFLAGS="${CXXFLAGS} -I/usr/local/lib/python3.10/site-packages/numpy/core/include" && \
    ./configure --enable-eigen --enable-armadillo --enable-openmp --enable-fftw --enable-python
RUN make blank-doc
RUN make 
RUN make install

RUN pip install h5py
ENV LD_LIBRARY_PATH /usr/local/lib
# RUN acol -h ## acol -h executes successfully here
ENV HDF5_DIR=/usr/lib/x86_64-linux-gnu/hdf5/serial
RUN make o2scl-test


FROM python:3.10 as o2sclpy

RUN git clone https://github.com/awsteiner/o2sclpy /opt/o2sclpy && \
    cd /opt/o2sclpy && \
    git checkout 7125dd645f3f0bfd334e8f25a66e325dab571579
WORKDIR /opt/o2sclpy
RUN pip install .


FROM python:3.10

COPY requirements.txt requirements.txt
RUN pip install -r requirements.txt

COPY --from=o2sclpy /usr/local/lib/python3.10/site-packages/o2sclpy /usr/local/lib/python3.10/site-packages/o2sclpy
COPY --from=o2scl /opt/o2scl /opt/o2scl
