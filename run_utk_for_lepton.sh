#!/bin/bash

# Runs the UTK project locally
# Requires all necessary shared object dynamic libraries to be installed on the local machine, see docs for more information

set -euo pipefail
echo -e "\nRunning UTK module locally...\n"

# Determine the path to the Python executable (preferably python3)
PYTHON="$(command -v python3 2>/dev/null || echo python)"

# Default file paths
USER_CONFIG_YAML_PATH="api/input/config.yaml"
POTENTIAL_DATA_HDF5_PATH="data/fid_3_5_22.o2"

# Check if command-line arguments are given to overwrite defaults
if [ $# -ge 1 ]; then
    USER_CONFIG_YAML_PATH="$1"
fi
if [ $# -ge 2 ]; then
    POTENTIAL_DATA_HDF5_PATH="$2"
fi

# Convert the file paths to absolute paths
USER_CONFIG_YAML_PATH=$(realpath "$USER_CONFIG_YAML_PATH")
POTENTIAL_DATA_HDF5_PATH=$(realpath "$POTENTIAL_DATA_HDF5_PATH")

# Create the 'input' and 'output' directories if they do not already exist
mkdir -p api/input
mkdir -p api/output

# Check if user config file exists
if [ ! -f "$USER_CONFIG_YAML_PATH" ]; then
    echo "YAML configuration file does not exist: $USER_CONFIG_YAML_PATH"
    exit 1
fi

# Check if the user config file is not in the expected location; copy it if needed.
if [ "$USER_CONFIG_YAML_PATH" != "$(realpath "api/input/config.yaml")" ]; then
    cp "$USER_CONFIG_YAML_PATH" api/input/config.yaml
fi

# Check if the EOS file exists
if [ ! -f "$POTENTIAL_DATA_HDF5_PATH" ]; then
    echo "EOS data file does not exist: $POTENTIAL_DATA_HDF5_PATH"
    exit 1
fi

# Check if the EOS file is in the expected location
if [ "$POTENTIAL_DATA_HDF5_PATH" != "$(realpath "data/$(basename "$POTENTIAL_DATA_HDF5_PATH")")" ]; then
    echo "Error: EOS data file is not in data/ directory: $POTENTIAL_DATA_HDF5_PATH"
    exit 1
fi

# ----------------------------------------------------------------
# New EOS parameter sets
# ----------------------------------------------------------------

P_FIDUCIAL="470 738 0.5 13.0 62.4 32.8 0.9"
P_LARGE_MMAX="783 738 0.5 13.0 62.4 32.8 0.9"
P_SMALL_R="214 738 0.5 13.0 62.4 32.8 0.9"
P_SMALLER_R="256 738 0.5 13.0 62.4 32.8 0.9"
P_LARGE_R="0 738 0.5 13.0 62.4 32.8 0.9"
P_SMALL_SL="470 738 0.5 13.0 23.7 29.5 0.9"
P_LARGE_SL="470 738 0.5 13.0 100.0 36.0 0.9"

# Run UTK module for Lepton
./eos_nuclei \
		-select-model $P_FIDUCIAL \
		-set a_virial 10 -set b_virial 10 \
		-set extend_frdm 0 \
		-set fd_A_max 600 -set max_ratio 7.0 \
		-set fixed_dist_alg 1999 \
		-set function_verbose 0 \
		-load data/fid_3_5_22.o2 \
		-utk-for-lepton create

# Check exit status
if [ $? -eq 0 ]; then
  echo -e "\n\tUtk running for Lepton: OK\n"
else
  echo -e "\n\tUtk running for Lepton: Failed\n"
  exit 1
fi

echo -e "\nUtk running for Lepton completed\n"
exit 0
