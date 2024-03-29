FROM --platform=linux/x86_64 ubuntu:20.04
SHELL ["/bin/bash", "-c"]

# install build tools
RUN apt update \
    && DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC\
    apt install -y sudo wget build-essential \
    gfortran openmpi-bin libopenmpi-dev \
    git software-properties-common curl \
    libblas-dev liblapack-dev

#install python
RUN add-apt-repository ppa:deadsnakes/ppa \
    && apt update \
    && apt install -y python3.10 python3-pip python3.10-distutils \
    && update-alternatives --install /usr/bin/python3 python /usr/bin/python3.10 1 \
    && curl -sS https://bootstrap.pypa.io/get-pip.py | python3.10

RUN pip install ase pytest gitpython

# install fftw3
RUN wget https://www.fftw.org/fftw-3.3.10.tar.gz\
    && tar -zxvf fftw-3.3.10.tar.gz \
    && cd fftw-3.3.10 \
    && ./configure \
    && make -j\
    && sudo make install

# compile state
COPY make-arch .
ARG state_src
COPY ${state_src} .

RUN tar xzf state-5.6.10.tgz \
    && cd ./state-5.6.10/src \
    && ln -fs ../../make-arch make.arch \
    && make


