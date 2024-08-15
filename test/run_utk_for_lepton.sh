#!/bin/bash

# Runs the Crust-DFT locally
# Requires all necessary shared object dynamic libraries to be installed on the local machine, see docs for more information

set -euo pipefail
echo -e "\nRunning Crust-DFT module...\n"

# Determine the path to the Python executable (preferably python3)
PYTHON="$(command -v python3 2>/dev/null || echo python)"

# Default file paths
USER_CONFIG_YAML_PATH="../input/config.yaml"

python3 ../src/Status.py \
    200 "Starting Crust-DFT module" 

# Check if command-line arguments are given to overwrite defaults
if [ $# -ge 1 ]; then
    USER_CONFIG_YAML_PATH="$1"
fi

# Convert the file paths to absolute paths
USER_CONFIG_YAML_PATH=$(realpath "$USER_CONFIG_YAML_PATH")

# Check if user config file exists
if [ ! -f "$USER_CONFIG_YAML_PATH" ]; then
    echo "YAML configuration file does not exist: $USER_CONFIG_YAML_PATH"
    python3 ../src/Status.py \
    400 "Error: config file is not in input/ directory"
fi

# Check if the user config file is not in the expected location; copy it if needed.
if [ "$USER_CONFIG_YAML_PATH" != "$(realpath "../input/config.yaml")" ]; then
    cp "$USER_CONFIG_YAML_PATH" ../input/config.yaml
fi

# ----------------------------------------------------------------
# validate the config.yaml file
$PYTHON ../src/yaml_validator.py
# ----------------------------------------------------------------
# Read validated config.yaml values for Crust-DFT module
read_parameters() {
    while IFS="=" read -r name value; do
        # Remove leading/trailing whitespaces
        name=$(echo "$name" | sed 's/^[ \t]*//;s/[ \t]*$//')
        value=$(echo "$value" | sed 's/^[ \t]*//;s/[ \t]*$//')

        # Replace invalid characters in the name with underscores
        name=$(echo "$name" | tr -cd '[:alnum:]_')

        # Set Bash parameters
        case $name in
            "output_format" | "nB_grid_spec" | "Ye_grid_spec" | "inc_lepton" | "verbose")
                eval "${name}=\$value"
                ;;
        esac
    done < <(awk -F':' '
        {
            # Remove leading/trailing whitespaces
            gsub(/^[ \t]+|[ \t]+$/, "", $1);
            gsub(/^[ \t]+|[ \t]+$/, "", $2);

            # Print parameter name and value
            printf("%s=%s\n", $1, $2);
        }
    ' "$1")
}

read_parameters "../input/validated_config.yaml"
# ----------------------------------------------------------------
# Default EOS table paths
EOS_DATA_HDF5_PATH=../data/EOS_table.o2

# Check if command-line arguments are given to overwrite defaults
if [ $# -ge 2 ]; then
    EOS_DATA_HDF5_PATH="$2"
fi

# Convert the EOS table paths to absolute paths
EOS_DATA_HDF5_PATH=$(realpath "$EOS_DATA_HDF5_PATH")

# Check if the EOS file exists
if [ ! -f "$EOS_DATA_HDF5_PATH" ]; then
    echo "EOS data file does not exist: $EOS_DATA_HDF5_PATH"
    python3 ../src/Status.py \
    400 "Error: EOS data file is not in data/ directory"
fi

# Check if the EOS file is in the expected location
if [ "$EOS_DATA_HDF5_PATH" != "$(realpath "../data/$(basename "$EOS_DATA_HDF5_PATH")")" ]; then
    echo "Error: EOS data file is not in data/ directory: $EOS_DATA_HDF5_PATH"
    python3 ../src/Status.py \
    400 "Error: EOS data file is not in data/ directory" \
    exit 1
fi

# ---------------------------------------------------------------
# Run Crust-DFT module
../src/eos_nuclei \
		-load $EOS_DATA_HDF5_PATH \
        -set nB_grid_spec $nB_grid_spec \
        -set Ye_grid_spec $Ye_grid_spec \
        -set inc_lepton $inc_lepton \
        -set verbose $verbose \
		-muses create 

# Run Postprocess.py
$PYTHON ../src/postprocess.py

# Check exit status
if [ $? -eq 0 ]; then
  echo -e "\n\tCrust-DFT module: Success\n"
  python3 ../src/Status.py \
    200 "Success: Crust-DFT module: Success" 
else
  echo -e "\n\tCrust-DFT module: Failed\n"
  python3 ../src/Status.py \
    400 "Error: Crust-DFT module: Failed"
  exit 1
fi

echo -e "\nCrust-DFT module completed\n"
exit 0
