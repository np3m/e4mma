#!/bin/bash

set -ex

## Change to parent directory of this script
cd "$(dirname "$(readlink -f "$0")")"

## Temporarily copy the source files into the Docker build context
rsync -vcr --delete ../src/ ./module/

## Build the Docker image
docker build . -t local-module-docs

## Remove temporary copy of source files
rm -rf ./module/

## Run the Docker image to serve the rendered HTML 
## at http://localhost:4000/docs/
docker run --rm -p 4000:80 local-module-docs
