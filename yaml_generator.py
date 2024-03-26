#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import yaml
"""
file: yaml_generator.py
This file creates a YAML configuration file for lepton-module based on user input.
Run as:
python3 yaml_generator.py --option1=value --option2=value --option3=value ...
If an option is not included it will be set to default value.
"""

# Constant
DEFAULT_CONFIG_FILE = "config.yaml"

def str2bool(v):
    return v.lower() in ('true', 't', '1', 'yes', 'y')

def create_yaml_file(input, temp_file_path):
    """
    Create a temporary YAML file with the provided command line arguments.

    Args:
        args: Parsed command line arguments.
    """
    # Data dictionary to store parsed arguments
    data = {
            "add_eg": input.add_eg,
            "alias": input.alias,
            "alt_model": input.alt_model,
            "check_virial": input.check_virial,
            "commands": input.commands,
            "compare_tables": input.compare_tables,
            "create_new_table": input.create_new_table,
            "edit_data": input.edit_data,
            "eg_table": input.eg_table,
            "eos_deriv": input.eos_deriv,
            "eos_sn": input.eos_sn,
            "exit": input.exit,
            "fit_frdm": input.fit_frdm,
            "fix_cc": input.fix_cc,
            "generate_table": input.generate_table,
            "get": input.get,
            "help": input.help,
            "hrg_load": input.hrg_load,
            "increase_density": input.increase_density,
            "interp_point": input.interp_point,
            "license": input.license,
            "load": input.load,
            "maxwell": input.maxwell,
            "mcarlo_beta": input.mcarlo_beta,
            "mcarlo_data": input.mcarlo_data,
            "mcarlo_neutron": input.mcarlo_neutron,
            "mcarlo_nuclei": input.mcarlo_nuclei,
            "mcarlo_nuclei2": input.mcarlo_nuclei2,
            "merge_tables": input.merge_tables,
            "muses": input.muses,
            "muses_table": input.muses_table,
            "no_intro": input.no_intro,
            "output": input.output,
            "pns_eos": input.pns_eos,
            "point": input.point,
            "point_nuclei": input.point_nuclei,
            "quit": input.quit,
            "random": input.random,
            "run": input.run,
            "select_high_T": input.select_high_T,
            "select_model": input.select_model,
            "shell": input.shell,
            "stability": input.stability,
            "stats": input.stats,
            "table_Ye": input.table_Ye,
            "table_full": input.table_full,
            "table_nB": input.table_nB,
            "test_deriv": input.test_deriv,
            "test_eg": input.test_eg,
            "test_random": input.test_random,
            "verify": input.verify,
            "vir_comp": input.vir_comp,
            "vir_fit": input.vir_fit,
            "warranty": input.warranty,
            "write_nuclei": input.write_nuclei,
            "xml_to_o2": input.xml_to_o2,
            "set": {
            "S_grid_spec": input.S_grid_spec,
            "T_grid_spec": input.T_grid_spec,
            "Ye_grid_spec": input.Ye_grid_spec,
            "Ye_list": input.Ye_list,
            "a_virial": input.a_virial,
            "alg_mode": input.alg_mode,
            "b_virial": input.b_virial,
            "cs2_verbose": input.cs2_verbose,
            "edge_list": input.edge_list,
            "ext_guess": input.ext_guess,
            "extend_frdm": input.extend_frdm,
            "fd_A_max": input.fd_A_max,
            "file_update_iters": input.file_update_iters,
            "file_update_time": input.file_update_time,
            "fixed_dist_alg": input.fixed_dist_alg,
            "function_verbose": input.function_verbose,
            "inc_hrg": input.inc_hrg,
            "include_detail": input.include_detail,
            "include_muons": input.include_muons,
            "max_ratio": input.max_ratio,
            "max_time": input.max_time,
            "mh_tol_rel": input.mh_tol_rel,
            "nB_grid_spec": input.nB_grid_spec,
            "ns_record": input.ns_record,
            "nucleon_func": input.nucleon_func,
            "old_ns_fit": input.old_ns_fit,
            "propagate_points": input.propagate_points,
            "recompute": input.recompute,
            "rnuc_less_rws": input.rnuc_less_rws,
            "select_cs2_test": input.select_cs2_test,
            "show_all_nuclei": input.show_all_nuclei,
            "six_neighbors": input.six_neighbors,
            "strange_axis": input.strange_axis,
            "survey_eqs": input.survey_eqs,
            "test_ns_cs2": input.test_ns_cs2,
            "use_alt_eos": input.use_alt_eos,
            "verbose": input.verbose,
            "verify_only": input.verify_only,    
        },
            "output_format": input.output_format,
        

    }

    # Write the data to a temporary YAML file
    with open(temp_file_path, "w") as outfile:
        yaml.dump(data, outfile, default_flow_style=False)


def main():
    """
    Entry point for the YAML preprocessing script.
    """

    # Read OpenAPI specification and extract the variables from the schema, their default and description



    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Create a YAML configuration file for lepton module based on user input validated by the OpenAPI specification",
    add_help=False)

    parser.add_argument("--config_path", type=str, default= "api/input/", help="Path to the config file")

    # Computational Parameters
    parser.add_argument("--add_eg",type=None,default=None,help="Add electrons and photons.")
    parser.add_argument("--alias",type=str,default='',help="Create a command alias.")
    parser.add_argument("--alt_model",type=str,default='',help="Use alternate, rather than the Du et al. EOS.")
    parser.add_argument("--check_virial",type=None,default=None,help="Check the virial EOS.")
    parser.add_argument("--commands",type=None,default=None,help="List all available commands.",)
    parser.add_argument("--compare_tables",type=str,default="filename1 filename2",help="Compare two output tables.",)
    parser.add_argument("--create_new_table",type=None,default=None,help="Create new table for Muses at fixed T",)
    parser.add_argument("--edit_data",type=str,default="<select func.> [tensor to modify] [value func.]",help="Edit an EOS table.",)
    parser.add_argument("--eg_table",type=str,default="filename",help="Construct an electrons and photon table.",)
    parser.add_argument("--eos_deriv",type=None,default=None,help="Compute derivatives numerically.",)
    parser.add_argument("--eos_sn",type=str,default="0.0",help="Compare to other EOSs?",)
    parser.add_argument("--exit",type=None,default=None,help="Exit (synonymous with 'quit').",)
    parser.add_argument("--fit_frdm",type=None,default=None,help="Fit the FRDM mass model.",)
    parser.add_argument("--fix_cc", type=str, default="output file", help="Increase nB to optimize the phase transition.",)
    parser.add_argument("--generate_table", type=str, default="output file", help="Generate an EOS table.",)
    parser.add_argument("--get", type=str, default='', help="Get the value of a parameter.",)
    parser.add_argument("--help", type=str, default='', help="Show help information.",)
    parser.add_argument("--hrg_load", type=str, default="filename", help="Load the HRG table.",)
    parser.add_argument("--increase_density", type=str, default="<nB low> <nB high> <Ye low> <Ye high> <T low> <T high> <output file>", help="Use low densities to improve results at high densities.",)
    parser.add_argument("--interp_point", type=None, default=None, help="interpolate a point.")
    
    parser.add_argument(
        "--license", 
        type=str, 
        default="filename", 
        help="Show license information."
    )

    parser.add_argument(
        "--load", 
        type=str, 
        default="data/fid_3_5_22.o2", 
        help="Load an EOS table."
    )

    parser.add_argument(
        "--maxwell", 
        type=str, 
        default="nB Ye T", 
        help="Maxwell construction."
    )
   
    parser.add_argument(
        "--mcarlo_beta",
        type=str,
        default="filename npoint",
        help="Monte Carlo neutrino opacity in beta equilibrium.",
    )

    parser.add_argument(
        "--mcarlo_data",
        type=str,
        default='params',
        help="Compute the data for the Monte Carlo figures.",
    )
    
    parser.add_argument(
        "--mcarlo_neutron",
        type=str,
        default='<filename> [n_point]',
        help="Monte Carlo neutrino opacity in pure neutron matter.",
    )

    parser.add_argument(
        "--mcarlo_nuclei",
        type=str,
        default='params',
        help="Monte Carlo results with nuclei.",
    )
    parser.add_argument(
        "--mcarlo_nuclei2",
        type=str,
        default='<nB> <Ye> <T> <N> <filename>',
        help="Monte Carlo results with nuclei (v2)",
    )
    parser.add_argument(
        "--merge_tables",
        type=str,
        default='<input file 1> <input file 2> <output file>',
        help="Merge two output tables to create a third.",
    )

    parser.add_argument(
        "--muses",
        type=None,
        default=None,
        help="Create a new table for Muses (experimental)",
    )
    parser.add_argument(
        "--muses_table",
        type=None,
        default=None,
        help="Create a new EOS table for Muses (experimental)",
    )
    parser.add_argument(
        "--no_intro",
        type=None,
        default=None,
        help="Do not print introductory text.",
    )

    # electron neutrino parameters
    parser.add_argument(
        "--output",
        type=str,
        default='filename',
        help="Output an EOS table to a file.",
    )
    parser.add_argument(
        "--pns_eos",
        type=str,
        default="<entropy per baryon> <lepton fraction> <output filename>",
        help="Construct the PNS EOS and the M-R curve.",
    )
    parser.add_argument(
        "--point",
        type=str,
        default="params",
        help="Evaluate the EOS at one (nB,Ye,T) point.",
    )

    parser.add_argument(
        "--point_nuclei",
        type=str,
        default="0.16 0.5 30",
        help="Compute and/or show EOS results at one (n_B,Y_e,T) point.",
    )

    parser.add_argument(
        "--quit",
        type=None,
        default=None,
        help="Quit (synonymous with 'exit').",
    )

    parser.add_argument(
        "--random",
        type=None,
        default=None,
        help="Select a random EOS model.",
    )

    parser.add_argument(
        "--run",
        type=str,
        default="filename",
        help="Run a file containing a list of commands.",
    )

    # muon parameters
    parser.add_argument(
        "--select_high_T",
        type=int,
        default=0,
        help="Select the high-temperature Skyrme EOS.",
    )
    parser.add_argument(
        "--select_model",
        type=str,
        default="skyrme model vector of numbers",
        help="Select an EOS model.",
    )
    parser.add_argument(
        "--set",
        type=str,
        default="param value",
        help="Set the value of a parameter.",
    )

    # muon neutrino parameters
    parser.add_argument(
        "--shell",
        type=str,
        default="command",
        help="Run a shell command.",
    )
    parser.add_argument(
        "--stability",
        type=str,
        default="filename",
        help="Compute the second derivatives and the stability matrix.",
    )
    parser.add_argument(
        "--stats",
        type=None,
        default=None,
        help="Output convergence statistics and simple checks.",
    )

    parser.add_argument(
        "--table_Ye",
        type=str,
        default="<filename> <Ye>",
        help="Construct a table at fixed electron fraction.",
    )
    parser.add_argument(
        "--table_full",
        type=str,
        default="filename",
        help="Construct a full 3D EOS table.",
    )
    parser.add_argument(
        "--table_nB",
        type=str,
        default="filename 0.16",
        help="Construct a table at fixed baryon density.",
    )

    parser.add_argument(
        "--test_deriv",
        type=str,
        default="params",
        help="Test the first derivatives of the free energy.",
    )
    parser.add_argument(
        "--test_eg",
        type=str,
        default="filename",
        help="Test the electron and photon EOS.",
    )
    parser.add_argument(
        "--test_random",
        type=int,
        default=5,
        help="Test an EOS at random points in (nB,Ye,T)",
    )
    parser.add_argument(
        "--verify",
        type=str,
        default='random <n_tests> <output file>',
        help="Verify the EOS.",
    )
    parser.add_argument(
        "--vir_comp",
        type=str,
        default="params",
        help="Compare the virial and full EOS.",
    )
    parser.add_argument(
        "--vir_fit",
        type=None,
        default=None,
        help="Fit the virial EOS.",
    )
    parser.add_argument(
        "--warranty",
        type=None,
        default=None,
        help="Show warranty information.",
    )
    parser.add_argument("--write-nuclei",type=str,default="filename",help="Write the nuclear masses to an HDF5 file.")
    parser.add_argument("--xml_to_o2", type=None, default=None, help="Read doxygen XML files to create run-time documentation.")
    parser.add_argument("--S_grid_spec", type=str, default="5,0.1*i", help="The function for S grid specification.")
    parser.add_argument("--T_grid_spec", type=str, default="160,0.1*1.046^i", help="The function for default temperature grid.")
    parser.add_argument("--Ye_grid_spec", type=str, default="70,0.01*(i+1)", help="The function for default electron fraction grid.")
    parser.add_argument("--Ye_list", type=str, default="1-3,5-7,59-60", help="The list of electron fractions for generate-table command.")
    parser.add_argument("--a_virial", type=float, default=10.0, help="Coefficient for modulation of virial EOS.")
    parser.add_argument("--alg_mode", type=int, default=4, help="Algorithm mode.")
    parser.add_argument("--b_virial", type=float, default=10.0, help="Coefficient for modulation of virial EOS.")
    parser.add_argument("--cs2_verbose", type=int, default=0, help="Verbose parameter for cs2.")
    parser.add_argument("--edge_list", type=str2bool,default=False, help="If true, recompute points where Z, A, log_xn, or log_xp reach a maximum or minimum.")
    parser.add_argument("--ext_guess", type=str, default="", help="Filename containing a separate table to use as a guess for the generate-table command.")
    parser.add_argument("--extend_frdm", type=str2bool,default=False, help="If true, attempt to extend FRDM beyond the drip lines.")
    parser.add_argument("--fd_A_max", type=int, default=600, help="The maximum value of A for a fixed distribution when alg_mode is 4.")
    parser.add_argument("--file_update_iters", type=int, default=100000, help="The number of iterations between file updates.")
    parser.add_argument("--file_update_time", type=int, default=1.8e3, help="The time (in seconds) between output file updates for generate_table().")
    parser.add_argument("--fixed_dist_alg", type=int, default=1111, help="Algorithm for eos_fixed_dist().")
    parser.add_argument("--function_verbose", type=int, default=11111, help="A new function verbose parameter.")
    parser.add_argument("--inc_hrg", type=str2bool,default=False, help="Incremental HRG parameter.")
    parser.add_argument("--include_detail", type=str2bool,default=False, help="Include detail parameter.")
    parser.add_argument("--include_muons", type=str2bool,default=False, help="If true, include muons.")
    parser.add_argument("--max_ratio", type=float, default=7.0, help="The maximum value of N/Z or Z/N.")
    parser.add_argument("--max_time", type=float, default=0.0, help="Maximum time, in seconds, for the generate-table command.")
    parser.add_argument("--mh_tol_rel", type=float, default=1.0e-6, help="Relative tolerance for the solver in the eos_fixed_dist() function.")
    parser.add_argument("--nB_grid_spec", type=str, default="301,10^(i*0.04-12)*2.0", help="The function for default baryon density grid.")
    parser.add_argument("--ns_record", type=str2bool,default=False, help="If true, save the results of the neutron star fit to a file, and immediately exit.")
    parser.add_argument("--nucleon_func", type=str, default="(i<100)*10+(i>=100)*sqrt(i)", help="Function for delta Z and delta N in the single nucleus approximation.")
    parser.add_argument("--old_ns_fit", type=str2bool,default=True, help="If true, use the old neutron star fit.")
    parser.add_argument("--propagate_points", type=str2bool,default=True, help="If true, use previously computed points (or guesses) as an initial guess to compute adjacent points.")
    parser.add_argument("--recompute", type=str2bool,default=False, help="If true, recompute all points, irrespective of the value of the convergence flag.")
    parser.add_argument("--rnuc_less_rws", type=str2bool,default=True, help="If true, ensure that the nuclear radius is less than the Wigner-Seitz radius.")
    parser.add_argument("--select_cs2_test", type=str2bool,default=True, help="If true, test cs2 in the select_internal() function.")
    parser.add_argument("--show_all_nuclei", type=str2bool,default=False, help="If true, show all nuclei considered at every point.")
    parser.add_argument("--six_neighbors", type=str2bool,default=False, help="If true, when computing a point, perform the calculation using some of the neighboring points as an initial guess.")
    parser.add_argument("--strange_axis", type=str2bool,default=False, help="If true, strangeness is included as an additional tensor rank")
    parser.add_argument("--survey_eqs", type=str2bool,default=False, help="If true, survey the nB and Ye equations near a failed point.")
    parser.add_argument("--test_ns_cs2", type=str2bool,default=False, help="If true, test the neutron star speed of sound.")
    parser.add_argument("--use_alt_eos", type=str2bool,default=False, help="Use an alternate EOS rather than the Du et al. combined EOS.")
    parser.add_argument("--verbose", type=int, default=0, help="Verbose parameter.")
    parser.add_argument("--verify_only", type=str2bool,default=False, help="If true, verify points only and do not attempt to solve.")
    parser.add_argument("--output_format", type=str, default='CSV', help="Output format for Lepton module")

    args = parser.parse_args()

    create_yaml_file(args, args.config_path+DEFAULT_CONFIG_FILE)

if __name__ == "__main__":
    print("\nStarting execution of yaml_generator.py..")

    main()

    print("\nExecution end of yaml_generator.py..")
