FROM ubuntu:23.10

# mkokot_TODO: change repo path and remove credentials management
# add LABEL org.opencontainers.image.source=https://github.com/refresh-bio/splash

# echo -n <github_token> > github_token.txt"
# export DOCKER_BUILDKIT=1
# sudo docker build --no-cache --build-arg GITHUB_USER=<github_username> --secret id=github_token,src=github_token.txt -t splash:1.9.0 .
# rm github_token.txt

# docker tag splash:1.9.0 ghcr.io/refresh-bio/splash:1.9.0

# docker push ghcr.io/refresh-bio/splash:1.9.0

# to login if push does not work:
# 		CR_PAT=github_token
# 		"echo $CR_PAT | docker login ghcr.io -u USERNAME --password-stdin"

# to create tar image: docker save splash:1.9.0 -o splash_1.9.0.tar
# convert to singularity (for example on Sherlock) image: singularity build splash_1.9.0.sif docker-archive://splash_1.9.0.tar
# run example:
# sudo docker run -v $(pwd):/home/ubuntu r-splash:1.9.0 splash input.txt
# example run sudo docker run -v `pwd`:/home/ubuntu splash --gap_len 0 --calculate_stats  --dump_Cjs --n_most_freq_targets 10  --pvals_correction_col_name pval_opt input.txt
LABEL version="1.9.0"

# Declare the build arguments
ARG GITHUB_USER

WORKDIR /home/ubuntu
RUN --mount=type=secret,id=github_token apt-get update -y \
    && apt-get install -y make python3 g++ wget git time gfortran \
    && git clone https://${GITHUB_USER}:$(cat /run/secrets/github_token)@github.com/marekkokot/splash_0.2.2_to_1.9.0/ \
    && cd splash_0.2.2_to_1.9.0 \
    && make -j \
    && cp example/splash.py /usr/bin/splash \
    && cp bin/* /usr/bin \
    && cd .. \
    && rm -rf splash_0.2.2_to_1.9.0 \
    && apt-get remove -y make git wget g++ \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
