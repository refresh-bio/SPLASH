#!/usr/bin/env python3
import build_release
import os

if __name__ == "__main__":
    build_release.run_cmd("git submodule init")
    build_release.run_cmd("git submodule update")
    # https://stackoverflow.com/questions/15715825/how-do-you-get-the-git-repositorys-name-in-some-git-repository
    repo = build_release.run_cmd_get_stdout("basename -s .git `git config --get remote.origin.url`").strip().lower()
    with open("hash.git", "w") as f:
        f.write(build_release.run_cmd_get_stdout("git rev-parse HEAD"))
    
    ver = build_release.get_ver("src/splash.py")
    cmd = f"sudo docker build --no-cache -t {repo}:{ver} ."
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
