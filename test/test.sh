#!/bin/bash

# Runs the Crust-DFT project locally
# Requires all necessary shared object dynamic libraries to be installed on the local machine, see docs for more information

set -euo pipefail
echo -e "\nRunning Crust-DFT module locally...\n"

# Determine the path to the Python executable (preferably python3)
PYTHON="$(command -v python3 2>/dev/null || echo python)"

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

# Run Crust-DFT module
../src/eos_nuclei \
	-select-model $P_FIDUCIAL \
	-point-nuclei 0.16 0.5 0.1 

# Check exit status
if [ $? -eq 0 ]; then
  echo -e "\n\tCrust-DFT local test: OK\n"
else
  echo -e "\n\tCrust-DFT local test: Failed\n"
  exit 1
fi

echo -e "\nCrust-DFT local test completed\n"
exit 0
