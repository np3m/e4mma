#!/usr/bin/env python
# -*- coding: utf-8 -*-

## Validate the config.yaml to ensure that it meets the OpenAPI specifications
## If any fields are left empty, the config.yaml will be overwritten with these fields filled by default values

import argparse
import os
import yaml
import json
from openapi_core import Spec
from openapi_core import unmarshal_request
from openapi_core.testing import MockRequest

# Constants
DEFAULT_API_FILE_PATH = os.path.join(os.path.dirname(__file__), '../api/OpenAPI_Specifications_UTK.yaml')
DEFAULT_CONFIG_FILE_PATH = os.path.join(os.path.dirname(__file__), '../input/config.yaml')
DEFAULT_VALID_CONFIG_FILE_PATH = os.path.join(os.path.dirname(__file__), '../input/validated_config.yaml')


def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Validate and unmarshal YAML configuration file for UTK module based on the OpenAPI specification"
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
        default=DEFAULT_CONFIG_FILE_PATH,
        help="Path to the initial configuration file",
    )

    parser.add_argument(
        "--valid_config_file_path",
        type=str,
        default=DEFAULT_VALID_CONFIG_FILE_PATH,
        help="Path to the final validated configuration file",
    )

    args = parser.parse_args()

    # Validate input config file against OpenAPI spec and create validated config
    assert validate_unmarshal_file(spec_file_path=args.api_file_path, input_file_path=args.config_file_path, valid_file_path=args.valid_config_file_path)

def validate_unmarshal_file(spec_file_path, input_file_path, valid_file_path, verbose=True):
    # Reading OpenAPI specifications
    with open(spec_file_path, 'r') as fp:
        openapi_specs = Spec.from_file(fp)
    
    # Reading input file
    with open(input_file_path, 'r') as fp:
        data = yaml.safe_load(fp)
    
    # Define request
    request = MockRequest(
        host_url='/',
        method="PUT",
        path='/input/config.yaml',
        data=json.dumps(data)
    )

    # Validate and unmarshal request
    try:
        valid_result = unmarshal_request(spec=openapi_specs, request=request)
        
    except Exception as err:
        if verbose:
            print(f'''Invalid config file "{input_file_path}" against OpenAPI spec file "{spec_file_path}":''')
            errors = f'''{err}'''
            for error in err.__cause__.schema_errors:
                #~Skip "None for not nullable" error as it always come double with another error, for when no value is given for a property
                if error.message != "None for not nullable":
                    #~For errors associated with a given defined property
                    if len(error.relative_path) != 0:
                        errors += "\n   ['"+ error.relative_path[0] +"'] -> "+ error.message
                    else:
                        errors += "\n  -> "+ error.message
            print(f'''{errors}''')
        return False
    
    # Write unmarshaled request to output file
    with open(valid_file_path, 'w') as fp:
        yaml.safe_dump(valid_result.body, fp)
    return True


if __name__ == "__main__":
    print("\nStarting execution of yaml_validation.py...")
    main()
    print("\nFinished yaml_validation.py...")