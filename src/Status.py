import argparse
import os
import yaml
from openapi_core import Spec

# Default paths
DEFAULT_API_FILE_PATH = os.path.join(os.path.dirname(__file__), '..', 'api/status_Specifications.yaml')
DEFAULT_STATUS_FILE_PATH = os.path.join(os.path.dirname(__file__), '..', 'output/status.yaml')

def main():
    
    # Create command line argument parser
    parser = argparse.ArgumentParser(
        description="Create a YAML configuration file for E4MMA module based on user input conforming to the OpenAPI specification"
    )

    # Paths
    parser.add_argument(
        "--api_file_path",
        type=str,
        default=DEFAULT_API_FILE_PATH,
        help="Path to the OpenAPI specification file",
    )

    parser.add_argument(
        "--config_file_path",
        type=str,
        default=DEFAULT_STATUS_FILE_PATH,
        help="Path to the created configuration file",
    )

    (args, rest) = parser.parse_known_args()

    # Parse command line arguments based on OpenAPI spec
    args = parse_args_from_spec(parser, spec_file_path=args.api_file_path)

    # Write parsed command line arguments to YAML configuration file based on OpenAPI spec
    write_args_to_config(args, spec_file_path=args.api_file_path, config_file_path=args.config_file_path)


def parse_args_from_spec(argparser, spec_file_path):
    
    # Reading OpenAPI specifications
    with open(spec_file_path, 'r') as fp:
        openapi_specs = Spec.from_file(fp)
    
    # Get schema for config.yaml file
    config_path = openapi_specs / "components" / "schemas" / "Status" / "properties"

    # Recursively iterate through nested schema and create a command line argument for each property
    def add_arguments(node, required=[]):
        for prop, value in node.items():
            if value['type'] == 'object':
                add_arguments(value['properties'], required=value['required'] if 'required' in value else [])
            else:
                argparser.add_argument(
                    "--" + prop,
                    type=str,
                    required=True if prop in required else False,
                    help=value['description']
                )

    with config_path.open() as config:
        add_arguments(config)
    
    # Parse command line arguments
    return argparser.parse_args()


def write_args_to_config(args, spec_file_path, config_file_path):

    # Reading OpenAPI specifications
    with open(spec_file_path, 'r') as fp:
        openapi_specs = Spec.from_file(fp)
    
    # Get schema for config.yaml file
    config_path = openapi_specs / "components" / "schemas" / "Status" / "properties"

    # Recursively append arguments to dictionary based on OpenAPI schema and cast argument types as they are placed
    def append_arguments(node, data):
        for prop, value in node.items():
            if value['type'] == 'object':
                data[prop] = {}
                append_arguments(value['properties'], data[prop])
            else:
                if vars(args)[prop]:
                    data[prop] = type_cast(vars(args)[prop], value['type'])
                
    data = {}
    with config_path.open() as config:
        append_arguments(config, data)
    
    # Write configuration data to config.yaml file
    with open(config_file_path, 'w') as fp:
        yaml.safe_dump(data, fp)


def type_cast(value, type_name):
    if type_name == 'number':
        return float(value)
    elif type_name == 'integer':
        return int(value)
    elif type_name == 'boolean':
        if value.lower() in ("true", "1"):
            return True
        elif value.lower() in ("false", "0"):
            return False
    else:
        return value


if __name__ == "__main__":
    main()
