#!/bin/bash

# Runs the UTK project locally
# Requires all necessary shared object dynamic libraries to be installed on the local machine, see docs for more information

set -euo pipefail
echo -e "\nRunning UTK EoS locally...\n"

# Determine the path to the Python executable (preferably python3)
PYTHON="$(command -v python3 2>/dev/null || echo python)"

# Default file paths
USER_CONFIG_YAML_PATH="../input/point.yaml"
EOS_DATA_HDF5_PATH="../data/fid_3_5_22.o2"

# Check if command-line arguments are given to overwrite defaults
if [ $# -ge 1 ]; then
    USER_CONFIG_YAML_PATH="$1"
fi
if [ $# -ge 2 ]; then
    EOS_DATA_HDF5_PATH="$2"
fi

# Convert the file paths to absolute paths
USER_CONFIG_YAML_PATH=$(realpath "$USER_CONFIG_YAML_PATH")
EOS_DATA_HDF5_PATH=$(realpath "$EOS_DATA_HDF5_PATH")

# Create the 'input' and 'output' directories if they do not already exist
mkdir -p ../input
mkdir -p ../output

# Check if user config file exists
if [ ! -f "$USER_CONFIG_YAML_PATH" ]; then
    echo "YAML configuration file does not exist: $USER_CONFIG_YAML_PATH"
    exit 1
fi

# Check if the user config file is not in the expected location; copy it if needed.
if [ "$USER_CONFIG_YAML_PATH" != "$(realpath "../input/point.yaml")" ]; then
    cp "$USER_CONFIG_YAML_PATH" ../input/point.yaml
fi

# Check if the EOS file exists
if [ ! -f "$EOS_DATA_HDF5_PATH" ]; then
    echo "EOS data file does not exist: $EOS_DATA_HDF5_PATH"
    exit 1
fi

# Check if the EOS file is in the expected location
if [ "$EOS_DATA_HDF5_PATH" != "$(realpath "../data/$(basename "$EOS_DATA_HDF5_PATH")")" ]; then
    echo "Error: EOS data file is not in data/ directory: $EOS_DATA_HDF5_PATH"
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

# ----------------------------------------------------------------
# validate the config.yaml file
$PYTHON ../src/point_validator.py
# Define the validated YAML file
yaml_file="../input/validated_point.yaml"
# ----------------------------------------------------------------
# Read config.yaml values to use to run utk code
read_parameters() {
    while IFS="=" read -r name value; do
        # Remove leading/trailing whitespaces
        name=$(echo "$name" | sed 's/^[ \t]*//;s/[ \t]*$//')
        value=$(echo "$value" | sed 's/^[ \t]*//;s/[ \t]*$//')

        # Replace invalid characters in the name with underscores
        name=$(echo "$name" | tr -cd '[:alnum:]_')

        # Set Bash parameters
        case $name in
            "load" | "select_high_T" | "alt_model" | "select_model" | "a_virial" | "b_virial" | "include_muons" | "max_ratio" | "mh_tol_rel" | "recompute" | "strange_axis" | "use_alt_eos" | "eos_deriv")
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

read_parameters $yaml_file

# Define the command names you want to extract
command_names=("point_nuclei" "point" "get" "random")

# Parse the YAML file and extract the command and its value
command=$(awk '/^ *commands:/{flag=1; next} flag && /^ *[^ ]/{sub(/:$/, "", $1); print $1; exit}' "$yaml_file")
value=$(awk '/^ *commands:/{flag=1; next} flag && /^ *[^ ]/{print substr($0, index($0,$2)); exit}' "$yaml_file")

# ----------------------------------------------------------------
# Run UTK module point
../src/eos_nuclei \
		-select-model $select_model -load $load \
        -select-high-T $select_high_T \
		-set a_virial $a_virial -set b_virial $b_virial \
        -set include_muons $include_muons -set max_ratio $max_ratio \
        -set mh_tol_rel $mh_tol_rel -set recompute $recompute \
        -set strange_axis $strange_axis -set use_alt_eos $use_alt_eos \
        -set verbose 1 \
		-$command $value
        
# ----------------------------------------------------------------
# Run Postprocess.py
#$PYTHON ../src/postprocess.py

# Check exit status
if [ $? -eq 0 ]; then
  echo -e "\n\tUtk running point: OK\n"
else
  echo -e "\n\tUtk running point: Failed\n"
  exit 1
fi

echo -e "\nUtk running point completed\n"
exit 0
