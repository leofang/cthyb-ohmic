#!/bin/bash
# this script should be placed at the root dirctory of cthyb-ohmic when executing

# TODO: test if this scripts is in the same folder as CMakeLists.txt
# TODO: test if root aceess is granted
# TODO: test if docker experimental feature is turned on (for squash)

# actual build 
docker build --build-arg arg_build_time="$(date +"%b. %d, %T %Z, %Y")"  --squash -t leofang/cthyb-ohmic .

# # push to Docker Hub
# # TODO: test if has access to docker://leofang/cthyb-ohmic
# docker push leofang/cthyb-ohmic
