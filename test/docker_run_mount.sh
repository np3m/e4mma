#!/bin/bash

# Runs the Crust-DFT module in a Docker container
set -euo pipefail
echo -e "\nRunning Crust-DFT module in Docker...\n"

# Default values for Docker image and tag
DOCKER_IMAGE_NAME="nostrad1/utk-eos"
DOCKER_IMAGE_TAG="v1.9.1"

# Get UID and GID of the current user
#UID=$(id -u)
GID=$(id -g)

# Check if command-line arguments are given to overwrite defaults
if [ $# -ge 1 ]; then
    DOCKER_IMAGE_NAME="$1"
fi
if [ $# -ge 2 ]; then
    DOCKER_IMAGE_TAG="$2"
fi

# Run the Crust-DFT Docker container, mapping input, output and data directories,
# mounting eos table as a volume, and executing in test directory.
docker run -it --rm --name utk -u $UID:$GID \
  -v "${PWD}/../input:/opt/e4mma/input:rw" \
  -v "${PWD}/../output:/opt/e4mma/output:rw" \
  -v "${PWD}/../data:/opt/e4mma/data" \
  $DOCKER_IMAGE_NAME:$DOCKER_IMAGE_TAG /bin/bash run_utk_for_lepton.sh

echo -e "\nCrust-DFT Docker run completed\n"
exit 0