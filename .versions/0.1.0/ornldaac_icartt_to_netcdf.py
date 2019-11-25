#!/usr/bin/env python3
__version__ = (0, 1, 0)
__updated__ = "2019-11-23"
__contact__ = "mcnelisjj@ornl.gov"
__doc__ = '''
#
#  VERSION:  {} [ {} ]
#  CONTACT:  {}
#
#  DESCRIPTION
#
#     This script provides all the functionality you need to translate
#     ICARTT v2 into netCDF. It parses standard metadata from from an 
#     ICARTT header and provides functions to translate to CF compliant
#     netCDF based on inputs and metadata given in ancillary JSON 
#     reference files.
#
#  USAGE
#
#     ...
#
#  NOTES
#
#  * parser is based on icartt v2.0 spec. most convenient source, imo:
#    https://www-air.larc.nasa.gov/missions/etc/IcarttDataFormat.htm
#
#  * and spec reference on earthdata;s page:
#    https://cdn.earthdata.nasa.gov/conduit/upload/6158/ESDS-RFC-029v2.pdf
#
#  * important Python constructs are named in CAPS.
#
#
'''.format(
    __version__,
    ".".join(list(
    __updated__ )),
    __contact__,
)

import os
import re
import sys
import yaml
import json
import numpy as np
import pandas as pd
import netCDF4 as nc4

# Get the local path of this script.
script_path = os.path.dirname(os.path.abspath(__file__))


def _printer(message: str):
    '''
    Prints a message to inform user about some event in a standardized way.

    Parameters
    ----------
    message (str): A notification that gets sandwiched between two newlines.

    '''
    print("\n - {}\n".format(message))


def _exit_with_error(message: str="(unexpected) notify: mcnelisjj@ornl.gov"):
    '''
    Prints the input error message and exits the script with a failure.

    Parameters
    ----------
    message (str): An error message that gets sandwiched between two newlines.

    '''
    sys.exit(print("ERROR: {} -- Abort.".format(message)))


def _crawl_directory(path: str, extension: str=None):
    """
    Crawl an input directory for a list of ICARTT files.

    Parameters
    ----------
    path (str): full path to an input directory.
    ext (str): An optional extension to limit search.

    Returns:
    -------
    A list of paths to data files (strings).

    """

    selected_files = []

    # Walk directory.
    for root, dirs, files in os.walk(path):

        # Loop over files,
        for f in files:

            # Get the extension.
            fext = os.path.splitext(f)[1]

            # If file matches input extension or if no extension given,
            if extension is None or extension==fext:

                # Join to root for the full path.
                fpath = os.path.join(root, f)

                # Add to list.
                selected_files.append(fpath)

    # Return the complete list.
    return selected_files


def __write_yaml(dictionary: dict, output: str):
    '''
    
    Parameters
    ----------
    
    Returns
    -------

    '''

    # Dump to YAML string.
    dict_yaml = yaml.dump(dictionary)

    # Write output YAML.
    with open(output, "w") as f:
        f.write(dict_yaml)

    
def _write_json(dictionary: dict, output: str, numpy_cls=True):
    '''
    
    Parameters
    ----------
    
    Returns
    -------

    '''

    # Dump to JSON string.
    if numpy_cls:
        dict_json = json.dumps(dictionary, indent=2, cls=NumpyEncoder)
    else:
        dict_json = json.dumps(dictionary, indent=2)

    # And write output JSON.
    with open(output, "w") as f:
        f.write(dict_json)


def __read_netcdf_file_structure(input_structure: str):
    '''
    Description ...

    Parameters
    ----------

    Return
    ------
    s
    '''

    # If it's a dict, return.
    if type(input_structure) is dict:
        return input_structure

    # If type is string and extension is JSON, parse and return.
    if type(input_structure) is str:
        if input_structure.endswith(".json"):
            with open(input_path_json, "r") as f:
                return json.load(f)
        
        # Else if it's a YMAL, parse and return.
        elif input_structure.endswith(".yaml"):
            with open(input_path_yaml, "r") as f:
                return yaml.load(f)

    # Else return None.
    else:
        return None


def _write_variable_jsons(input_structure: str, PATH_OUTPUT: str):
    '''
    Writes one JSON per variable for all variables in a project
    Parameter
    ---------


    '''
    
    # Get dataset_structure as a dictionary if it is a file path.
    netcdf_dataset_structure = _read_netcdf_file_structure(input_structure)

    # Loop over the variables in the structure dictionary and write to json.
    for name, variable in netcdf_dataset_structure['variables'].items():
        _write_json(variable, os.path.join(PATH_OUTPUT, name + ".json"))

    # Call this function on each group.
    for name, group in netcdf_dataset_structure['groups'].items():
        _parse_netcdf_dataset_structure_to_variable_jsons(group)


def parse_icartt_header(icartt_file_path: str):
    '''
    Parse all of the scientifically relevant information from an ICARTT header.

    Parameters:
    -----------
    icartt_file_path (str): path to a standard ICARTT v2.0 file

    Returns
    -------
    A tuple with a variable table in position one (with columns: name, units, 
    scale factor, and fill value) and a dictionary of other relevant metadata
    from the ICARTT header.

    '''

    # Open with context manager for reading,
    with open(icartt_file_path, "r") as f:

        # and get the first line of the file.
        line1 = f.readlines()[0]

    # Select the header length from the index values.
    hdr_length = int(line1.split(",")[0])

    # Get the header index (the line number of header of the table).
    hdr_index = hdr_length - 1

    # Now open again for reading with context manager,
    with open(icartt_file_path, "r") as f:

        # and read entire header's as lines into list.
        hdr_lines = f.readlines()[1:hdr_index]

    # Strip the end of line characters from all header lines.
    hdr_lines = [ln.strip() for ln in hdr_lines]

    #
    #  Parse the rest of the ICARTT header based on the v2.0 specification.
    #  Note that Python sequences are base 0.
    #

	# 2. PI last name, first name/initial.
    pi_last_name, pi_first_name = hdr_lines[0].split(", ")

	# 3. Organization/affiliation of PI.
    pi_organization = hdr_lines[1]

	# 4. Data source description (e.g., instrument, platform, model, etc.).
    data_source_description = hdr_lines[2]

	# 5. Mission name (usually the mission acronym).
    mission_name = hdr_lines[3]

	# 6. File volume number, number of file volumes.
    file_volume_number, file_volume_count = hdr_lines[4].split(", ")

	# 7. UTC date start, UTC date end (yyyy, mm, dd, yyyy, mm, dd).
    utc_start = "-".join(hdr_lines[5].split(", ")[:3])
    utc_stop = "-".join(hdr_lines[5].split(", ")[-3:])

    # 8. Interval (this is complicated; see specification).
    interval = hdr_lines[6].split(",")

	# 9. Description or name of independent variable (see specification).
    independent_var_name, independent_var_units = hdr_lines[7].split(",")

	# 10. Number of dependent variables.
    dependent_var_count = int(hdr_lines[8])

    #
    # The next two lines contain the scale factors and fill value for each
    # dependent variable. We should add an integer 1 for the scale factor
    # and None for the fill value of the independent variable.
    #

	# 11. Scale factors (typically 1; given for every column). Make float.
    scale_factors = [1] + [float(s) for s in hdr_lines[9].split(", ")]

	# 12. Missing data indicators. Make float.
    fill_values = [None] + [float(s) for s in hdr_lines[10].split(", ")]

    #
    #  The variable names and units are listed line by line from here. Loop
    #  over the lines after line 12 until the dictionary is the same length
    #  as the integer 'dependent_var_count' + 1.
    #
    #  Store the variable names and units in a dictionary where keys are 
    #  the integer positions of the columns in left-right order.
    #

    var_names_units = dict()

    # I'm faily confident that independent variable is always first.
    var_names_units[0] = (independent_var_name, independent_var_units)

    # Loop over the lines. 
    for i, ln in enumerate(hdr_lines[11:]):
        
        # Break loop when line with no comma is encountered.
        # if ", " not in ln:
        #    break
        if len(var_names_units)==(dependent_var_count + 1):
            break
        
        # Split line and store a tuple with column indices for dict keys.
        else:
            var_names_units[i + 1] = tuple(ln.split(", "))

    #
    #  The line that immediately follows the variable, unit pairs is an 
    #  integer that indicates the number of 'special comment lines', and it
    #  is followed by the 'special comments'. We'll keep those as a block
    #  of text for now and worry about it later if necessary.
    #  
    #  The only exception is the revision info. We loop over the remaining 
    #  lines and see if one starts with 'REVISION:'. If one does, get the 
    #  revision ID (like 'R#:') and then look for a line describing that 
    #  revision ID, which should start with the ID. 
    #

    revision_id, revision_description = None, None

    # Loop over the lines. 
    for i, ln in enumerate(hdr_lines[11 + dependent_var_count:]):
        
        # If we already found the revision_id, look for it to get description.
        if revision_id is not None and ln.startswith(revision_id):
            revision_description = ln.split(": ")[1]

        # Else if it starts with 'REVISION:', then get the revision number.
        elif ln.startswith("REVISION: "):
            revision_id = ln.split(": ")[1]

        # Else if both are none, keep looping.
        elif revision_id==None and revision_description==None:
            continue

        # Else if neither are none, break.
        elif revision_id!=None and revision_description!=None:
            break

        # Else: probably a bug or there are no revisions. Come back to this.
        else:
            pass

    #
    #  Now combine all of the per-variable details into a table with columns:
    #   1. variable index (column number, base zero),
    #   2. name, 
    #   3. unit, 
    #   4. scale factor,
    #   5. fill value,
    #

    # Loop over the variable names dictionary.
    for i, v in var_names_units.items():

        # Select the scale_factor and fill_values.
        scale = scale_factors[i]
        fill = fill_values[i]

        # Replace tuple in dictionary with a four-part tuple.
        var_names_units[i] = (v[0], v[1], scale, fill)

    # Turn the dictionary's values into a pandas table.
    variable_table = pd.DataFrame(

        # Turn the dictionary values into a list of tuples.
        list(var_names_units.values()),

        # The values of each tuple correspond to one variable's metadata.
        columns=['name', 'units', 'scale', 'fill']
    
    )

    #  Put everything in a dictionary and return it.
    return(variable_table, dict(
        header_length = hdr_length, 
        header_index = hdr_index, 
        pi_last_name = pi_last_name, 
        pi_first_name = pi_first_name,
        pi_organization = pi_organization,
        data_source_description = data_source_description,
        mission_name = mission_name,
        file_volume_number = file_volume_number,
        file_volume_count = file_volume_count,
        utc_start = utc_start,
        utc_stop = utc_stop,
        interval = interval,
        revision_id = revision_id,
        revision_description = revision_description,
    ))


def write_resource_files(DATA: dict):
    '''
    ...

    Parameters
    ----------

    '''

    # Get resources dir from the config.
    resources_dir = DATA['output_ancillary_files']

    # If there is no resources dir, create it.
    if not os.path.isdir(resources_dir):
        os.mkdir(resources_dir)

    # Master variable list.
    variable_list = []

    # Loop over the ICARTT files.
    for f in DATA['icartt_files']:

        # Parse the header of the ICARTT.
        variables, metadata = parse_icartt_header(f)

        # Get the basename of the ICARTT without extension.
        file_name = os.path.splitext(os.path.basename(f))[0]

        # Join the output path.
        file_output = os.path.join(resources_dir, "{}".format(file_name))

        # Write the test variable table to a CSV in tests/
        variables.to_csv(file_output + ".csv")

        # Also dump the header data (excluding table) to a JSON file.
        with open(file_output + ".json", "w") as f:
            f.write(json.dumps(metadata, indent=2))

        # Extend the list of all variables.
        variable_list.extend(variables.name.tolist())

    # Get the list of unique variable names.
    variable_names = sorted(list(set(variable_list)))

    ### Update translation table (key table ict -> nc).
    with open(os.path.join(resources_dir + "translate.csv"), "w") as f:
        f.write("\n".join(variable_names))

    # If there is no netcdf_variables/ subdir, create it.
    netcdf_variables_dir = os.path.join(resources_dir, "netcdf_variables/")
    if not os.path.isdir(netcdf_variables_dir):
        os.mkdir(netcdf_variables_dir)

    ### Write individual variable jsons.
    _write_variable_jsons(, netcdf_variables_dir)

    # If the file already exists read with pandas.
    #if os.path.isfile(translate_csv):
    #    translate_tbl = pd.read_csv(translate_csv)
    #else:
        #pd.DataFrame([(v, , , ) for v in variable_names])


def get_netcdf_dataset_structure(dataset):
    """
    Makes a JSON (Panoply-like) representation of netCDF structure. 
    Groups are also dataset constructs in netCDF4 python.
    
    Parameters
    ----------
    
    Returns
    -------
    
    """
    
    # netCDF structure as a dictionary. 
    structure = {
        'groups': {},
        'variables': {}, 
        'attributes': dataset.__dict__,
    }
    
    # VARIABLES: Add this dataset's variables to structure dictionary.
    for name, variable in dataset.variables.items():
        
        structure['variables'][name] = {
            'dimensions': variable.dimensions,
            'attributes': variable.__dict__,
        }

    # GROUPS: are also datasets. Loop and repeatedly call this function.
    for name, group in dataset.groups.items():
        structure['groups'][name] = get_netcdf_dataset_structure(group)
    
    return structure


def get_netcdf_file_structure(netcdf_file_path: str):
    '''
    
    Parameters
    ----------
    
    Returns
    -------
    
    '''
    
    # Read netCDF and get dictionary structure.
    with nc4.Dataset(netcdf_file_path) as ds:
        structure = get_netcdf_dataset_structure(ds)
    
    # Add dimensions to the dictionary structure.
    for name, dimension in dataset.dimensions.items():

        structure['dimensions'][name] = {
            'UNLIMITED': dimension.isunlimited(),
            'size': dimension.size()
        }
    
    return structure


def construct_PATH_OUTPUT(input_file_path: str, PATH_OUTPUT: str, output_file_extension: str):
    '''
    Simple path constructor.
    
    Parameters
    ----------
    
    Returns
    -------
    
    '''

    # Get output file basename.
    output_file_basename = os.path.splitext(os.path.basename(input_file_path))[0]
    
    # Get output file name.
    output_file_name = "{}.{}".format(output_file_basename, output_file_extension)

    # Use os module to join path.
    output_file_path = os.path.join(PATH_OUTPUT, output_file_name)

    return output_file_path


class NumpyEncoder(json.JSONEncoder):
    '''
    A well known hack for serializing numpy types to JSON.
    '''
    numpy_ints = (  # numpy integer data types.
        np.int_, np.intc, np.intp, np.int8, np.int16, np.int32, 
        np.int64, np.uint8, np.uint16, np.uint32, np.uint64)
    numpy_flts = (  # numpy float data types.
        np.float_, np.float16, np.float32, np.float64)
    def default(self, obj):
        if isinstance(obj, self.numpy_ints):
            return int(obj)
        elif isinstance(obj, self.numpy_flts):
            return float(obj)
        elif isinstance(obj,(np.ndarray,)): #### This is the fix
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def write_netcdf_file_structure(
    input_path_netcdf: str, 
    PATH_OUTPUT_json:str=None, 
    PATH_OUTPUT_yaml:str=None):
    '''
    Get the complete hierarchical structure of a netCDF file as a Python dict.
    
    Parameters
    ----------
    
    Returns
    -------
    
    '''
    
    # Open the nteCDF file.
    netcdf_dataset = nc4.Dataset(input_path_netcdf)
   
    # Get the structure of the dataset root.
    root_structure = get_netcdf_dataset_structure(netcdf_dataset)
    
    # Close the open netCDF dataset.
    netcdf_dataset.close()
       
    # If JSON path was given as arg, write to disk.
    if PATH_OUTPUT_json is not None:
        _write_json(root_structure, PATH_OUTPUT_json)

    # If YAML path was given as arg, write to disk.
    if PATH_OUTPUT_yaml is not None:
        __write_yaml(root_structure, PATH_OUTPUT_yaml)


def handle_configuration():
    '''

    Parameters
    ----------
    
    Returns
    -------

    '''

    ### The configuration file should exist in the same directory as script.
    config_path = os.path.join(script_path, "config.yml")

    #
    #  Parsing configuration YAML file: {}
    #

    # If it doesn't exist, inform the user and exit.
    if not os.path.isfile(config_path):
        _exit_with_error("No config.yml in script directory.")

    # Otherwise read it with the base Python YAML loader.
    with open(config_path, "r") as f:
        
        # Try to parse yaml and exit on failure. Notify on success.
        try:
            DATA = yaml.safe_load(f)
        except yaml.YAMLError as e:
            print(e)
            _exit_with_error("Could not parse YAML config file.")
        else:
            _printer("SUCCESS")

    #
    #  Validate the inputs given in the configuration file:
    #
    #  1. Input ICARTT directory: 
    #     {input_DIR_ICARTT}
    #
    #  2. Output netCDF directory:
    #     {PATH_OUTPUT}
    #
    #  3. Variable name translation reference table (CSV):
    #     {input_CSV_VARIABLES}
    #
    #  4. Variable metadata reference (JSON):
    #     {input_variable_metadata}
    #

    ### Ensure ICARTT dir is valid.

    if not os.path.isdir(DATA['input_DIR_ICARTT']):
        _exit_with_error("Input ICARTT directory is invalid.")

    # If the directory exists, recursively crawl it for ICARTT files.
    DATA['icartt_files'] = _crawl_directory(DATA['input_DIR_ICARTT'])

    # If no icartts were found, exit and notify user.
    if len(DATA['icartt_files'])==0:
        _exit_with_error("No ICARTT files found in the input directory.")
    
    # Else, inform on the number of ICARTT files
    else:
        _printer("Found [ {} ] ICARTTs.".format(len(DATA['icartt_files'])))

    ### Ensure netCDF dir is valid.

    if not os.path.isdir(DATA['PATH_OUTPUT']):
        _exit_with_error("Output netCDF directory is invalid.")

    ### Ensure netCDF file output structure is valid.

    # Read the output netCDF structure.
    output_netcdf_structure = _read_netcdf_file_structure(
        DATA['RESOURCES']['JSON_NETCDF'])

    #
    #  Validating the inputs given in the primary configuration file:
    #
    #  1. Input ICARTT directory: 
    #     {input_DIR_ICARTT}
    #
    #  2. Output netCDF directory:
    #     {PATH_OUTPUT}
    #
    #  3. Variable name translation reference table (CSV):
    #     {input_CSV_VARIABLES}
    #
    #  4. Variable metadata reference (JSON):
    #     {input_variable_metadata}
    #
    # ### Ensure that the input variable names translation table is vlaid.
    # if not os.path.isfile(DATA['variable_names']):
    #     _exit_with_error("Path to variable translation table is invalid.")
    # # Try to read with pandas and exit on failure. Notify user of success.
    # try:
    #     # var_names = pd.read_csv(DATA['variable_names'])
    #     DATA['variables'] = pd.read_csv(DATA['variable_names'])
    # except:
    #     _exit_with_error("Could not parse variable translation table (CSV).")
    # else:
    #     _printer("SUCCESS")
    # ### Ensure that input netcdf header structure is valid.    
    # if not os.path.isfile():
    #     _exit_with_error("Path to variable metadata JSON is invalid.")
    # # Try to parse as a json and exit on failure. Notify user of success.
    # try:
    #     var_metadata = json.loads(DATA['output_netcdf_structure'])
    # except:
    #     _exit_with_error("Could not parse variable metadata file (JSON).")
    # else:
    #     _printer("SUCCESS")
    # ### Consider other settings and notify the user.
    # if DATA['option_multileg_join']:
    #     _printer("Multileg flights WILL be joined to one output netCDF.")
    # else:
    #     _printer("Multileg flights WILL NOT be joined to one output netCDF.")
    #

    return(DATA)



def main(write_resources=True):
    '''
    Description ...

    Parameters
    ----------

    Returns
    -------

    '''

    # Open and dump to JSON for dev.
    DATA = handle_configuration()

    print('''
    #
    #  ## Parsing headers for all input ICARTT files.
    #  
    #  Must intercompare & find all unique variable names, collect their 
    #  netCDF attributes, and write to JSON. 
    # 
    #  ACT-America has multiple aircraft & variables differ between them.
    # 
    #
    #  **STEPS**
    #
    #   1. Parse the metadata from header of every ICARTT file into a dict.
    #
    #
    #''')

    # If write_resources is True, do it ~ a csv and json.
    if write_resources:

        write_
    

        

    
    return


def stat():
    '''
    - duration
    - average meteorology - winds, 
    '''
    return




def test():
    '''
    Description ...

    Parameters
    ----------

    Returns
    -------

    '''

    # Check and dump config to test dir.
    config = handle_configuration()
    with open("tests/config_test.json", "w") as f:
        f.write(json.dumps(config, indent=2))

    # (Make the test ICARTT file configurable).
    test_file = "tests/ACTAMERICA-mrgPFP-c130_merge_20180520_R0.ict"

    # Parse the header of the test ICARTT.
    test_file_variables, test_file_metadata = parse_icartt_header(test_file)

    # Get the full path and the basename of the ICARTT without the extension.
    test_file_path = os.path.dirname(test_file)
    test_file_name = os.path.splitext(os.path.basename(test_file))[0]

    # Join the output path.
    test_file_output = os.path.join(test_file_path, "{}".format(test_file_name))

    # Write the test variable table to a CSV in tests/
    test_file_variables.to_csv(test_file_output + ".csv")

    # Also dump the header data (excluding table) to a JSON file.
    with open(test_file_output + ".json", "w") as f:
        f.write(json.dumps(test_file_metadata, indent=2))



if __name__ == "__main__":
    # test()

    #
    #  Call argument handler. It should return a four item dictionary:
    #
    #   1 - arg1
    #   2 - arg2
    #

    # Pass the main_args to main using dictionary expansion.
    # main(**main_args)

    main()
    #config = handle_configuration()

    #print(config)

    