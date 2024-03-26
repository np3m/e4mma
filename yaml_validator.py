import sys
import argparse


DEFAULT_API_SCHEMA_PATH='api/OpenAPI_Specifications_UTK.yaml'
DEFAULT_CONFIG_PATH= 'input/config.yaml'
CONFIG_NAME= "config.yaml"

# Set path of api directory containing validation script
sys.path.append("api")
# Import validation script
from OpenAPI_Specifications_validator import OpenAPI_input_validator


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Validate and update YAML configuration.")
    parser.add_argument(
        '--openapi_schema', 
        help="Path to OpenAPI schema (YAML)", 
        default=DEFAULT_API_SCHEMA_PATH
    )

    parser.add_argument(
        "--input_config_path", 
        type=str, 
        default= "api/input/",
        help="Path to the config file"
    )

    args = parser.parse_args()
    # Run validation
    OpenAPI_input_validator(args.openapi_schema, args.input_config_path+CONFIG_NAME, warnings=True)


    print("\nFinished yaml_validator.py...")