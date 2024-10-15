FROM ubuntu:23.10

WORKDIR /home/ubuntu

#it seems docker v 24 works, but v 20 does not because of the rules in dockerignore
COPY . src
RUN  \
    apt-get update -y \
    && apt-get install -y cmake make python3 g++ wget time git \
    && cp -r src build \ 
	&& cd build \
	&& make -j \
	&& make install \
	&& cd .. \
	&& rm -rf build \
	&& apt-get remove -y cmake make wget g++ git \
	&& apt-get autoremove -y \
	&& apt-get clean \
	&& rm -rf src/.git \
	&& rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
