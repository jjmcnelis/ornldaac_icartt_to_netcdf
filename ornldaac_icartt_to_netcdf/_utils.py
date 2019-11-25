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
    __version__  , ".".join(list(
    __updated__  )),
    __contact__  ,
)

import os
import re
import sys
import yaml
import json
import numpy as np
import pandas as pd
import netCDF4 as nc4
from glob import glob


class NumpyEncoder(json.JSONEncoder):
    '''Well known scratch class for serializing numpy types to JSON.'''
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


def _write_yaml(dictionary: dict, output: str):
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


def read_netcdf_file_structure(input_structure: str):
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
            with open(input_structure, "r") as f:
                return json.load(f)
        
        # Else if it's a YMAL, parse and return.
        elif input_structure.endswith(".yaml"):
            with open(input_structure, "r") as f:
                return yaml.load(f)

    # Else return None.
    else:
        return None


def write_variable_jsons(input_structure: str, output_path: str):
    '''
    Writes one JSON per variable for all variables in a project
    Parameter
    ---------


    '''
    
    # Get dataset_structure as a dictionary if it is a file path.
    netcdf_dataset_structure = read_netcdf_file_structure(input_structure)

    # Loop over the variables in the structure dictionary and write to json.
    for name, variable in netcdf_dataset_structure['variables'].items():
        _write_json(variable, os.path.join(output_path, name + ".json"))

    # Call this function on each group.
    for name, group in netcdf_dataset_structure['groups'].items():
        write_variable_jsons(group)


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

    variable_names = []

    # Loop over the variable names dictionary.
    for i, v in var_names_units.items():

        # Select the scale_factor and fill_values.
        scale = scale_factors[i]
        fill = fill_values[i]

        variable_names.append(v[0])

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
        title = data_source_description,
        mission_name = mission_name,
        file_volume_number = file_volume_number,
        file_volume_count = file_volume_count,
        flight_date = utc_start,
        utc_stop = utc_stop,
        interval = interval,
        revision_id = revision_id,
        revision_description = revision_description,
        variable_list = variable_names,
    ))


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
            'datatype': str(variable.datatype),
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
        structure['dimensions'] = {}
        for name, dimension in ds.dimensions.items():
            structure['dimensions'][name] = {
                'UNLIMITED': dimension.isunlimited(),
                'size': dimension.size
            }
    
    return structure


def write_netcdf_file_structure(
    input_netcdf_path: str, 
    output_json_path: str=None):
    '''
    Writes a netCDF structure JSON to the output path, or the path of the input.
    '''

    # Get the dictionary structure of the input netCDF file.
    netcdf_structure = get_netcdf_file_structure(input_netcdf_path)

    # Open an output JSON for writing.
    if output_json_path is None:
        output_json_path = os.path.splitext(input_netcdf_path)[0]+".json"

    # Dump the netCDF structure dictionary to a JSON string and write to file.
    _write_json(netcdf_structure, output_json_path)


