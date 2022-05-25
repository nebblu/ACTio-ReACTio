FROM ubuntu:20.04

MAINTAINER Ben Bose <ben.bose@ed.ac.uk>

## Make sure building of image doesn't get stuck at any selection stages

ENV DEBIAN_FRONTEND=noninteractive


## Update apt get

RUN apt-get update && apt-get -y install


    #######################################


## GCC, G++ and GSL

RUN apt-get -y install gcc && \
    apt-get -y install g++ && \
    apt-get -y install libgsl-dev

    #######################################


## Python

RUN apt-get -y install python3 && \
    apt-get -y install python3-setuptools


    #######################################


## git and vim

RUN apt-get -y install git && \
    apt-get -y install vim


    #######################################


## Sundials

RUN apt-get -y install cmake && \
    mkdir sundials && cd sundials && \
    git clone -b release-4.1.0 https://github.com/LLNL/sundials.git && \
    mkdir instdir && \
    mkdir builddir && \
    cd builddir && \
    cmake -DCMAKE_INSTALL_PREFIX=/sundials/instdir \
     -DEXAMPLES_INSTALL_PATH=sundials/instdir/examples \
     ../sundials && \
     make && \
     make install && \
     cd


     #######################################

## React

RUN  git clone https://github.com/nebblu/ReACT.git && \
     cd ReACT && \
     sed -i "s/CPPFLAGS +=/CPPFLAGS += -I\/sundials\/instdir\/include/g" pyreact/Makefile && \
     sed -i "s/LDFLAGS +=/LDFLAGS += -L\/sundials\/instdir\/lib/g" pyreact/Makefile && \
     python3 setup.py develop

## set correct library paths in the container environment

ENV LD_LIBRARY_PATH /sundials/instdir/lib:/ReACT/reactions/lib
