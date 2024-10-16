#!/usr/bin/env python3
import build_release
import os
import sys

if __name__ == "__main__":
    if build_release.get_os() != "linux":
        print("Error: this currently works only on linux machine")
        sys.exit(1)

    hardware_platform = build_release.get_hardware()
    if hardware_platform == "x64":
        hardware_platform = "amd64"

    if hardware_platform not in ("amd64", "arm64"):
        print("Error: platform must be amd64 or arm64")
        sys.exit(1)

    build_release.run_cmd("git submodule update --init --recursive")
    # https://stackoverflow.com/questions/15715825/how-do-you-get-the-git-repositorys-name-in-some-git-repository
    repo = build_release.run_cmd_get_stdout("basename -s .git `git config --get remote.origin.url`").strip().lower()
    with open("hash.git", "w") as f:
        f.write(build_release.run_cmd_get_stdout("git rev-parse HEAD"))
    
    ver = build_release.get_ver("src/splash.py")
    #cmd = f"sudo docker build --no-cache -t {repo}:{ver} ."
    cmd = f"sudo docker buildx build --no-cache --platform linux/{hardware_platform} -t {repo}:{ver} --load ."

    print(cmd)
    build_release.run_cmd(cmd)
    
    os.remove("hash.git")

    print("Docker image created!")
    print(f"To export it to file run: docker save {repo}:{ver} -o {repo}_{ver}.tar")
    print(f"To convert to singularity run: singularity build {repo}_{ver}.sif docker-archive://{repo}_{ver}.tar")
    print(f"Example {repo} run from docker: sudo docker run -v $(pwd):/home/ubuntu {repo}:{ver} splash input.txt")

# Pushing to github packages:
# https://docs.github.com/en/packages/learn-github-packages/connecting-a-repository-to-a-package
#  docker tag {repo}:{ver} ghcr.io/refresh-bio/{repo}:{ver}
#  docker push ghcr.io/refresh-bio/{repo}:{ver}
# May require:
#  docker login ghcr.io -u USERNAME
# To use docker image in registry but from singularity:
# singularity run docker://ghcr.io/refresh-bio/{repo}:{ver} splash input.txt
# conver to sif: singularity build {repo}_{ver}.sif docker-archive://{repo}_{ver}.tar
