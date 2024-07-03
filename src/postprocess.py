#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import yaml

from muses_porter import Porter

"""
file: postprocess.py
author: Nikolas Cruz Camacho <cnc6@illinois.edu>
purpose: Postprocess UTK output files.

This script provides functionality to postprocess output files generated by the UTK model.
It includes operations such as slicing, copying, and processing specific output files
based on the defined variables and configurations.
"""

LEPTON_VARIABLES = [
    "temperature",
    "mu_B",
    "mu_S",
    "mu_Q",
    "baryon_density",
    "total_strangeness_density",
    "total_charge_density",
    "energy_density",
    "pressure",
    "entropy_density",
]


def load_configuration():
    """
    Load configuration from a config.yaml file.

    Returns:
        dict: A dictionary containing configuration data.
    """
    config_file = "../input/config.yaml"
    with open(config_file, "r") as file:
        yaml_data = yaml.load(file, Loader=yaml.FullLoader)
    return yaml_data


def main():
    """
    Main function for executing postprocess.py.

    This function processes the input folder based on user-provided arguments or default configurations.
    It copies specific UTK output files, creates an output directory if needed, and performs data slicing.

    Args:
        None

    Returns:
        None
    """
    print("\nStarting execution of postprocess.py...")

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Process input folder.")
    parser.add_argument(
        "run_name",
        nargs="?",
        default="default",
        help="Input folder to process. If none is provided, the default folder will be used.",
    )
    args = parser.parse_args()

    # Determine the run_name
    if len(sys.argv) == 1:
        yaml_input = load_configuration()
        run_name = yaml_input.get("run_name", "default")
        print(f"Running with a user-provided folder from config.yaml: {run_name}")
    elif len(sys.argv) == 2:
        run_name = args.run_name
    else:
        raise ValueError(
            "Too many arguments provided. Maximum one argument is allowed."
        )

    output_directory = "../output/"
    output_directory_run = output_directory

    # Create the output directory if it does not exist
    if not os.path.exists(output_directory_run):
        os.makedirs(output_directory_run)

    # Copy UTK output files
    utk_output_files = ["e4mma_wo_lepton.csv"]
    #for file_name in utk_output_files:
    #    shutil.copy(
    #        os.path.join(output_directory_run, file_name),
    #        os.path.join(output_directory, file_name),
    #    )

    # Create a Porter instance
    porter = Porter()

    # Process UTK output for Lepton
    print(f"Processing E4MMA output w/o Lepton: ...")

    porter.import_table(
        os.path.join(output_directory_run, f"e4mma_wo_lepton.csv"),
        extension="CSV",
        filename_schema="../api/OpenAPI_Specifications_UTK.yaml",
        schema="e4mma_wo_lepton",
        dropna=True,
        delimiter=",",
        verbose=False,
    )

    if yaml_input.get("output_format", "") == "HDF5":
        # Export data to an HDF5 file
        porter.export_table(
            os.path.join(
                output_directory, f"e4mma_wo_lepton.h5"
            ),
            filename_schema="../api/OpenAPI_Specifications_UTK.yaml",
            schema="e4mma_wo_lepton",
            extension="HDF5",
            dropna=True,
            verbose=False,
        )
    else:
        # Export data to an CSV file following the OpenAPI specifications UTK_output_for_Lepton
        porter.export_table(
            os.path.join(
                output_directory, f"e4mma_wo_lepton.csv"
            ),
            extension="CSV",
            filename_schema="../api/OpenAPI_Specifications_UTK.yaml",
            schema="e4mma_wo_lepton",
            dropna=True,
            delimiter=",",
            verbose=False,
            float_format="%.12e",
        )

    porter.reset()

    print("\nExecution end of postprocess.py.")


if __name__ == "__main__":
    main()
