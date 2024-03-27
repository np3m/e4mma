#!/usr/bin/env python
# -*- coding: utf-8 -*-
#-----------------------------------------------------------------------------------------------
#Script writen by J.Jahan (jahan.johannes@gmail.com), last update 2023/12/13
#version: 0.0.1
#-----------------------------------------------------------------------------------------------
# | MUSES |
# ---------
#	
################################################################################################

import sys as sys
import os as os
import yaml as yaml
import json as json
from openapi_core import Spec
from openapi_core import validate_request
from openapi_core.testing import MockRequest

#=========================#
# OPENAPI INPUT VALIDATOR #
#=========================#
def OpenAPI_input_validator(specs_file_path, input_file_path, warnings=False):    
    # Ensuring input files are found at specified locations
    #-------------------------------------------------------
    #Specs file
    if not os.path.exists(specs_file_path):
        print("[OpenAPI_input_validator]> File '" + specs_file_path + "' not found \n\nOPERATION ABORTED")
        sys.exit(3)
    
    #YAML file
    if not os.path.exists(input_file_path):
        print("[OpenAPI_input_validator]> File '" + input_file_path + "' not found \n\nOPERATION ABORTED")
        sys.exit(3)
    
    # Reading OpenAPI specifications
    #--------------------------------
    with open(specs_file_path, 'r') as spec_file:
        openapi_specs = Spec.from_file(spec_file)
    
    #~Extract name of input file
    input_file = input_file_path[input_file_path.rfind('/')+1:]
    
    #~Find corresponding path in specifications
    input_spec_path = ''
    for path in openapi_specs["paths"].keys():
        if input_file in path:
            input_spec_path = path

    #~Stop program if no match found
    if input_spec_path == '':
        print("[OpenAPI_input_validator]> No path found in specifications matching '" + input_file_path + "' \n\nOPERATION ABORTED")
        sys.exit(3)
    
    # Count number of expected parameters in input_file
    #---------------------------------------------------
    input_num_params = 0
    
    #~Extract paths to schemas in input file
    schemas = []
    specs_input_schemas = openapi_specs["paths"][input_spec_path]["put"]["requestBody"]["content"]["application/json"]["schema"]
    #~In case of several schemas
    if "allOf" in specs_input_schemas.keys():
        for schema in specs_input_schemas["allOf"]:
            schemas.append(schema["$ref"].split('/')[1:])
    #~In case of one schema
    else:
        schemas.append(specs_input_schemas["$ref"].split('/')[1:])
    
    #~Count number
    for schema in schemas:
        input_schema = openapi_specs
        #Define
        for path in schema:
            input_schema = input_schema[path]
        input_num_params += len(input_schema["properties"])
    
    # Reading input file
    #--------------------
    with open(input_file_path, 'r') as input_file:
        data = yaml.safe_load(input_file)
    
    # Define request
    #----------------
    request = MockRequest(
        host_url='/',
        method="PUT",
        path=input_spec_path,
        data=json.dumps(data)
    )

    #To store error messages if there are
    errors = "" 
    
    # Validate the request
    #----------------------
    try: #use 'try' to not interrupt the program in case 'validate_request()' finds errors
        #~Warning if extra parameters compared to specifications
        if warnings == True:
            if len(data) > input_num_params:
                print("\n   /!\\ Extra parameter(s) in input files compared to module specifications")
        #~Validate parameters that are properly defined
        validate_request(spec=openapi_specs, request=request)
    
    # Raise error message if problem in validation
    #----------------------------------------------
    except Exception as e:
        for error in e.__cause__.schema_errors:
            #~Skip "None for not nullable" error as it always come double with another error, for when no value is given for a property
            if error.message != "None for not nullable":
                #~For errors associated with a given defined property
                if len(error.relative_path) != 0:
                    errors += "\n   ['"+ error.relative_path[0] +"'] -> "+ error.message
                else:
                    errors += "\n  -> "+ error.message
    
    # Print message for validation 
    #------------------------------
    if errors != "":
        #~Failed validation
        errors = "Failed to validate "+ input_file_path +" with module specifications.\nRaised errors:" + errors
        print(errors)
        sys.exit(1)
    else:
        #~Successful validation
        print("Successful validation of %s"%(input_file_path))
        
    return 0

#======#
# MAIN #
#======#
if __name__ == '__main__':
    # Inputs 
    #--------    
    if len(sys.argv) == 3:
        #~Affect arguments 
        specs_file_path = sys.argv[1] 
        input_file_path = sys.argv[2]
    if len(sys.argv) < 3:
        #~Exit with error message if missing arguments
        print("[OpenAPI-Specifications_validator.py]> Missing arguments: should be 2 (specifications file + input file to validate)")
        sys.exit(2)
    
    # Validation
    #------------
    OpenAPI_input_validator(specs_file_path, input_file_path, warnings=True)

