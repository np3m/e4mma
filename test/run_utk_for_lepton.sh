#!/bin/bash

# Runs the E4MMA locally
# Requires all necessary shared object dynamic libraries to be installed on the local machine, see docs for more information

set -euo pipefail
echo -e "\nRunning E4MMA module...\n"

# Determine the path to the Python executable (preferably python3)
PYTHON="$(command -v python3 2>/dev/null || echo python)"

# Default file paths
USER_CONFIG_YAML_PATH="../input/config.yaml"
USER_STATUS_YAML_PATH="../output/status.yaml"

python3 ../src/Status.py \
    --code 200 \
    --message "Starting E4MMA running w/o Lepton" 

# Check if command-line arguments are given to overwrite defaults
if [ $# -ge 1 ]; then
    USER_CONFIG_YAML_PATH="$1"
fi

# Convert the file paths to absolute paths
USER_CONFIG_YAML_PATH=$(realpath "$USER_CONFIG_YAML_PATH")
USER_STATUS_YAML_PATH=$(realpath "$USER_STATUS_YAML_PATH")

# Create the 'input' and 'output' directories if they do not already exist
mkdir -p ../input
mkdir -p ../output

# Check if user config file exists
if [ ! -f "$USER_CONFIG_YAML_PATH" ]; then
    echo "YAML configuration file does not exist: $USER_CONFIG_YAML_PATH"
    echo "creating default configuration file: $USER_CONFIG_YAML_PATH"
    python3 ../src/yaml_generator.py \
	--verbose 0 \
	--load fid_3_5_22.o2 \
	--output_format HDF5 \
    --nB_grid_spec '301,10^(i*0.04-12)*2.0' \
	--Ye_grid_spec '70,0.01*(i+1)' 
fi

# Check if the user config file is not in the expected location; copy it if needed.
if [ "$USER_CONFIG_YAML_PATH" != "$(realpath "../input/config.yaml")" ]; then
    cp "$USER_CONFIG_YAML_PATH" ../input/config.yaml
fi

# ----------------------------------------------------------------
# validate the config.yaml file
$PYTHON ../src/yaml_validator.py
# ----------------------------------------------------------------
# Read validated config.yaml values for E4MMA module
read_parameters() {
    while IFS="=" read -r name value; do
        # Remove leading/trailing whitespaces
        name=$(echo "$name" | sed 's/^[ \t]*//;s/[ \t]*$//')
        value=$(echo "$value" | sed 's/^[ \t]*//;s/[ \t]*$//')

        # Replace invalid characters in the name with underscores
        name=$(echo "$name" | tr -cd '[:alnum:]_')

        # Set Bash parameters
        case $name in
            "load" | "output_format" | "nB_grid_spec" | "Ye_grid_spec" | "inc_lepton" | "verbose")
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
EOS_DATA_HDF5_PATH=../data/$load

# Check if command-line arguments are given to overwrite defaults
if [ $# -ge 2 ]; then
    EOS_DATA_HDF5_PATH="$2"
fi

# Convert the EOS table paths to absolute paths
EOS_DATA_HDF5_PATH=$(realpath "$EOS_DATA_HDF5_PATH")

# Check if the EOS file exists
if [ ! -f "$EOS_DATA_HDF5_PATH" ]; then
    echo "EOS data file does not exist: $EOS_DATA_HDF5_PATH"
    echo "Downloading EOS table: $EOS_DATA_HDF5_PATH"
    curl https://isospin.roam.utk.edu/public_data/eos_tables/du21/$load --output $EOS_DATA_HDF5_PATH
fi

# Check if the EOS file is in the expected location
if [ "$EOS_DATA_HDF5_PATH" != "$(realpath "../data/$(basename "$EOS_DATA_HDF5_PATH")")" ]; then
    echo "Error: EOS data file is not in data/ directory: $EOS_DATA_HDF5_PATH"
    python3 ../src/Status.py \
    --code 400 \
    --message "Error: EOS data file is not in data/ directory" \
    exit 1
fi

# ---------------------------------------------------------------
# Run E4MMA module
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
  echo -e "\n\tE4MMA running w/o Lepton: OK\n"
  python3 ../src/Status.py \
    --code 200 \
    --message "Error: E4MMA running w/o Lepton: OK" 
else
  echo -e "\n\tE4MMA running w/o Lepton: Failed\n"
  python3 ../src/Status.py \
    --code 400 \
    --message "Error: E4MMA running w/o Lepton: Failed"
  exit 1
fi

echo -e "\nE4MMA running w/o Lepton completed\n"
exit 0
