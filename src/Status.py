import argparse
import os
import yaml
from openapi_core import Spec

# Default paths
API_FILE_PATH = os.path.join(os.path.dirname(__file__), '..', 'api/OpenAPI_Specifications_UTK.yaml')
STATUS_FILE_PATH = os.path.join(os.path.dirname(__file__), '..', 'output/status.yaml')

def main():

    parser = argparse.ArgumentParser(description='Status Crust DFT EoS module execution')
    parser.add_argument('code', type=int, help='The code to write to the YAML file.')
    parser.add_argument('message', type=str, help='The message to write to the YAML file.')

    args = parser.parse_args()

    with open(STATUS_FILE_PATH, 'a') as statusfile:
        yaml.safe_dump({"code": args.code, "message": args.message}, statusfile)

if __name__ == "__main__":
    main()
