FROM ubuntu:24.04

WORKDIR /home/ubuntu

# designed to be used from build_docker.py, if not specified will be not set for `make`
# I'm using below syntax ${PLATFORM:+PLATFORM=$PLATFORM} which will add PLATFORM=... only if PLATFORM was defined
ARG PLATFORM

#it seems docker v 24 works, but v 20 does not because of the rules in dockerignore
COPY . src
RUN  \
    apt-get update -y \
    && apt-get install -y cmake make automake python3 g++ wget time git r-base \
	&& R -e "install.packages(c('data.table', 'glmnet', 'ggplot2', 'gridExtra', 'pheatmap'), repos='https://cloud.r-project.org')" \
    && cp -r src build \ 
	&& cd build \
	&& make ${PLATFORM:+PLATFORM=$PLATFORM} -j \
	&& make install \
	&& cd .. \
	&& rm -rf build \
	&& apt-get remove -y cmake make automake wget g++ git \
	&& apt-get autoremove -y \
	&& apt-get clean \
	&& rm -rf src/.git \
	&& rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
