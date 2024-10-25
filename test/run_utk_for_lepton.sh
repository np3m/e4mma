#!/bin/bash

# Runs the Crust-DFT locally
# Requires all necessary shared object dynamic libraries to be installed on the local machine, see docs for more information

set -euo pipefail
echo -e "\nRunning Crust-DFT module...\n"

# Determine the path to the Python executable (preferably python3)
PYTHON="$(command -v python3 2>/dev/null || echo python)"

# Default file paths
USER_CONFIG_YAML_PATH="../input/config.yaml"

$PYTHON ../src/Status.py \
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
    $PYTHON ../src/Status.py \
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
            "generate_table" | "ext_guess" | "output_format" | "nB_grid_spec" | "Ye_grid_spec" | "verbose")
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


# Check if the EOS file is in the expected location
#if [ "$EOS_DATA_HDF5_PATH" != "$(realpath "../data/$(basename "$EOS_DATA_HDF5_PATH")")" ]; then
#    echo "Error: EOS data file is not in data/ directory: $EOS_DATA_HDF5_PATH"
#    python3 ../src/Status.py \
#    400 "Error: EOS data file is not in data/ directory" \
#    exit 1
#fi
# ------------------------------------------------------------
# Run Crust-DFT module
# ------------------------------------------------------------
# ----------------------------------------------------------------
# EOS parameter sets from Du et al. (2022)
# ----------------------------------------------------------------

P_FIDUCIAL="470 738 0.5 13.0 62.4 32.8 0.9"
P_LARGE_MMAX="783 738 0.5 13.0 62.4 32.8 0.9"
P_SMALL_R="214 738 0.5 13.0 62.4 32.8 0.9"
P_SMALLER_R="256 738 0.5 13.0 62.4 32.8 0.9"
P_LARGE_R="0 738 0.5 13.0 62.4 32.8 0.9"
P_SMALL_SL="470 738 0.5 13.0 23.7 29.5 0.9"
P_LARGE_SL="470 738 0.5 13.0 100.0 36.0 0.9"
# ----------------------------------------------------------------
if [ $generate_table ]; then
    echo "Generating new table..."

    if [ $ext_guess ]; then
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
            $PYTHON ../src/Status.py \
            400 "Error: EOS data file is not in data/ directory"
        fi

        echo "With external guess: $EOS_DATA_HDF5_PATH" 
        ../src/eos_nuclei \
            -set data_dir "../data" \
            -select-model $P_FIDUCIAL \
            -set nB_grid_spec $nB_grid_spec \
            -set Ye_grid_spec $Ye_grid_spec \
            -set T_grid_spec "3,0.1+i" \
            -generate-table \
		    "ext_guess=$EOS_DATA_HDF5_PATH" \
		    -eos-deriv \
            -output ../output/Edt.o2 

    else
        echo "Without external guess" 
        ../src/eos_nuclei \
            -set data_dir "../data" \
            -select-model $P_FIDUCIAL \
            -set nB_grid_spec $nB_grid_spec \
            -set Ye_grid_spec $Ye_grid_spec \
            -set T_grid_spec "3,0.1+i" \
            -generate-table \
		    -eos-deriv \
            -output ../output/Edt.o2 
    fi
        ../src/eos_nuclei \
		    -load ../output/Edt.o2 \
            -set nB_grid_spec $nB_grid_spec \
            -set Ye_grid_spec $Ye_grid_spec \
            -set verbose $verbose \
		    -muses "../output/crust_dft.csv"


else
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
        $PYTHON ../src/Status.py \
        400 "Error: EOS data file is not in data/ directory"
    fi

    echo "Using precalculated table..."
    ../src/eos_nuclei \
		-load $EOS_DATA_HDF5_PATH \
        -set nB_grid_spec $nB_grid_spec \
        -set Ye_grid_spec $Ye_grid_spec \
        -set verbose $verbose \
		-muses "../output/crust_dft.csv"
fi
# ---------------------------------------------------------------
# Run Postprocess.py
if [ ! -f "../output/crust_dft.csv" ]; then
    echo "Output file does not exist: ../output/crust_dft.csv"
    $PYTHON ../src/Status.py \
    400 "Error: Output file is not in output/ directory"
fi
$PYTHON ../src/postprocess.py

# Check exit status
if [ $? -eq 0 ]; then
  echo -e "\n\tCrust-DFT module: Success\n"
  $PYTHON ../src/Status.py \
    200 "Success: Crust-DFT module: Success" 
else
  echo -e "\n\tCrust-DFT module: Failed\n"
  $PYTHON ../src/Status.py \
    400 "Error: Crust-DFT module: Failed"
  exit 1
fi

echo -e "\nCrust-DFT module completed\n"
exit 0
