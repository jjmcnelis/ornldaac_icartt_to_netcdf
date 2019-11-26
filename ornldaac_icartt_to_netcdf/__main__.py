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
from ornldaac_icartt_to_netcdf import *




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




def _construct_path_output(
    output_path: str, 
    input_file_path: str, 
    extension: str):
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
    output_file_name = "{}.{}".format(output_file_basename, extension)

    # Use os module to join path.
    output_file_path = os.path.join(output_path, output_file_name)

    return output_file_path




def _write_resource_files(DATA: dict):
    '''
    ...

    Parameters
    ----------

    '''

    # Get absolute path to resources directory.
    resources_dir = DATA['DIR_RESOURCES']

    # Make resources directory if it doesn't exist.
    if not os.path.isdir(resources_dir):
        os.mkdir(resources_dir)

    def resources_path_constructor(subdirectory: str):
        return os.path.join(resources_dir, subdirectory)

    # Update the paths to the resource subdirectories.
    for name, resource in DATA['RESOURCES'].items():
        DATA['RESOURCES'][name] = resources_path_constructor(resource)

        # Create the subdirectory if it doesn't exist.
        if not os.path.isdir(DATA['RESOURCES'][name]):
            print(" - Creating {} resources folder: {}".format(name, resource))
            os.mkdir(DATA['RESOURCES'][name])
    
    # Now grab the four paths for convenience.
    icartt_header_dir = DATA['RESOURCES']['ICARTT_HEADERS']
    icartt_variable_dir = DATA['RESOURCES']['ICARTT_VARIABLES']
 
    ### DECONSTRUCT ICARTT FILES IN A LOOP.

    print("\n### PARSING ICARTT HEADERS IN A LOOP AND WRITING RESOURCE FILES.")

    DATA['FLIGHT_RESOURCES'] = {}

    # Loop over the ICARTT files.
    for i, ict in enumerate(DATA['ICARTT_FILES']):
        print(" - [ {} / {} ] {} ".format(
            i+1, len(DATA['ICARTT_FILES']), os.path.basename(ict)))

        # Get the basename of the ICARTT without extension.
        file_name = os.path.splitext(os.path.basename(ict))[0]

        # Parse ICARTT header to a dictionary and a table.
        variables, metadata = parse_icartt_header(ict)

        ## Write ICARTT header dictionary to a JSON.
        ict_hdr = _construct_path_output(icartt_header_dir, file_name, "json")
        with open(ict_hdr, "w") as f:
            f.write(json.dumps(metadata, indent=2))

        ## Write ICARTT file variables to a CSV.
        ict_var = _construct_path_output(icartt_variable_dir, file_name, "csv")
        variables.to_csv(ict_var)

        # Add the paths as a dictionary linked to this ICARTT file.
        DATA['FLIGHT_RESOURCES'][ict] = {
            'HEADER': ict_hdr, 
            'VARIABLES': ict_var
        }

        ### BACKFILL NETCDF VARIABLES TO INCLUDE MISSING ICARTT VARS.

        # Loop over ICARTT's variable table as a list of dictionary records:
        for var in variables.to_dict(orient="records"):

            # Get the path to netCDF variable reference JSON.
            json_var = _construct_path_output(
                icartt_variable_dir, var['name'], "json")

            ## Write the current variable a JSON if there isn't one.
            if not os.path.isfile(json_var):
                with open(json_var, "w") as f:
                    f.write(json.dumps(var, indent=2))

    ### RECONCILE MASTER VARIABLES TABLE.

    csv_variables_header = [
        'VARIABLE_MAP',
        'NETCDF_NAME',
        'ICARTT_NAME',
        'ICARTT_UNITS',
        'ICARTT_SCALE',
        'ICARTT_FILL',
        'BOOL_INCLUDE',
        'USERFUNC',
    ]

    # Open a new file for writing at the master variables table path.
    var_csv = _construct_path_output(resources_dir, "VARIABLES", "csv")
    with open(var_csv, "w") as outf:

        # Write a header line.
        outf.write(",".join(csv_variables_header) + "\n")

        # Open the variable reference JSONs in a loop.
        for j in glob(icartt_variable_dir + "*.json"):
            with open(j, "r") as inf:
                
                # Get the input ICARTT variable metadata. Add NULL.
                ict_val = ['NULL']*2+list(json.load(inf).values())+[1, 'NULL']
                
                # Get the values as a list and convert to strings.
                out_val = [str(v) for v in ict_val]

                # Join with commas and write new lines to master variable CSV.
                outf.write(",".join(out_val) + "\n")




def _organize_standard_and_multileg_flights(DATA: dict):
    '''
    
    Parameters
    ----------
    
    Returns
    -------
    
    '''

    # A regular expression catches the multi leg flight suffix.
    multileg_regex = re.compile('_L[0-9].ict')

    # A dictionary stores the output filename and legs as child list.
    flights = {}

    for ict in DATA['ICARTT_FILES']:

        # If regular expression is not matched anywhere in string,
        if re.search(multileg_regex, ict) is None:

            # Add to list of standard flights.
            flights[ict] = ict

        # Else if regular expression is matched in string.
        else:

            # The output file won't have the suffix.
            output_filename = ict[:-7] + ".ict"

            # Add this file to the dict of multi-leg flights.
            if output_filename not in flights:
                flights[output_filename] = [ict]
            else:
                flights[output_filename].append(ict)

    # Return the organized flights as a dictionary.
    return flights




def parse_icartt_table(icartt_file: str):
    '''
    Parses a single ICARTT file to a pandas data frame.

    Parameters
    ----------

    Returns
    -------

    '''

    # Get the header row number from the ICARTT.
    with open(icartt_file, "r") as f:
        header_row = int(f.readlines()[0].split(",")[0])-1

    # Parse the table starting at the header index.
    _flight_table = pd.read_csv(icartt_file, header=header_row, delimiter=",")

    return _flight_table




def _parse_icartt_table_multileg(icartt_files: list):
    '''
    Reorganizes the data for a multileg flight to resemble the standard flight.

    Parameters
    ----------

    Returns
    -------

    '''

    # Sort the list of input ICARTTs.
    icartts = sorted(icartt_files)

    # This is the merged flight data.
    _merged_table = None

    # Loop over the ICARTTs for this flight.
    for ict in icartts:

        # Parse the table starting at the header index.
        _flight_table = parse_icartt_table(ict)

        # If flight_tables slot is None, set. Otherwise append.
        if _merged_table is None:
            _merged_table = _flight_table
        else:
            _merged_table.append(_flight_table)

    # Strip the leading/trailing white spaces from final merged tables.
    _merged_table.columns = [c.strip() for c in list(_merged_table.columns)]

    # Return the merged table.
    return _merged_table




def _parse_flight(icartt):
    '''
    Takes an input ICARTT file path or list of ICARTT file paths and parses
    to a single output data frame.

    Parameters
    ----------

    Returns
    -------

    '''

  

    # If the type(icartt) is a list, must be multileg.
    if type(icartt) is list:

        # Handle multileg flight so it is merged to single table.
        flight_data = _parse_icartt_table_multileg(icartt)

    # Else if the type is a string, parse the file.
    elif type(icartt) is str:

        # Call the parse_icartt_table directly.
        flight_data = parse_icartt_table(icartt)

    # Strip the leading/trailing white spaces from final merged tables.
    flight_data.columns = [c.strip() for c in list(flight_data.columns)]


    ### SPECIAL ACT-AMERICA ROUTINE -------------------------------------------

    # This special routine exists because the science team wants ot merge
    # some custom metadata flags with the merge ICARTT data set.
    
    _temp = icartt[0] if type(icartt) is list else icartt
    _tail = _temp.split("_merge_")[1].split("_R")[0]

    if "b200" in _temp:
        _flag = glob("inputs/ACTAMERICA_metadata_flags/B200/*{}*".format(_tail))
        if len(_flag)==1:
            _flag_table = parse_icartt_table(_flag[0])
        elif len(_flag)==2:
            _flag_table = _parse_icartt_table_multileg(_flag)
        else:
            print("\n\nESCAPED FLAG MERGE! FLAGS NOT WRITTEN TO OUTPUT!\n\n")
            _flag_table = None
    if "c130" in _temp:
        _flag = glob("inputs/ACTAMERICA_metadata_flags/C130/*{}*".format(_tail))
        if len(_flag)==1:
            _flag_table = parse_icartt_table(_flag[0])
        elif len(_flag)==2:
            _flag_table = _parse_icartt_table_multileg(_flag)
        else:
            print("\n\nESCAPED FLAG MERGE! FLAGS NOT WRITTEN TO OUTPUT!\n\n")
            _flag_table = None

    if _flag_table is not None:
        
        # Strip the leading/trailing white spaces from final merged tables.
        _flag_table.columns = [c.strip() for c in list(_flag_table.columns)]

        # If the index sizes match.
        if _flag_table.index.size==flight_data.index.size:

            # Loop over the flag variables.
            for _flag_var in ['Flight_flag', 
                              'Air_flag', 
                              'BL_FT_flag', 
                              'Maneuver_flag', 
                              'Maneuver_flagQC']:

                # Add the flag variables to the larger dataset.
                flight_data[_flag_var] = _flag_table[_flag_var]

    ### SPECIAL ACT-AMERICA ROUTINE -------------------------------------------


    # Return the parsed table.
    return flight_data




def write_netcdf(
    output_file, 
    output_structure,
    flight_data, 
    flight_variables, 
    flight_resources):
    '''
    Write a trajectory style netCDF from input table and references.

    Parameters
    ----------

    Returns
    -------

    '''

    ## Global attributes are templated from reference structure (STRUCTURE).
    global_attributes = output_structure['attributes']

    # Open the ICARTT header resource file (JSON) and read to dictionary.
    with open(flight_resources["HEADER"], "r") as f:

        # The ICARTT header replaces + adds to the global attributes.
        global_attributes.update(json.load(f))

    # Get the variable list from the global attributes.
    variable_list = global_attributes['variable_list']

    # Drop the variable list and a few others.
    del global_attributes['variable_list']
    del global_attributes['header_length']
    del global_attributes['header_index']
    del global_attributes['utc_stop']

    # Get the flight date and reformat the string.
    flight_datestring = 'merge_{}'.format(
        global_attributes['flight_date'].replace("-", ""))

    # Open a netCDF file for writing.
    with nc4.Dataset(output_file, "w") as ds:

        # Set the global attributes for the output netCDF.
        ds.setncatts(global_attributes)

        # Add the dimensions from the reference structure in a loop.
        for name, dimension in output_structure['dimensions'].items():
            #size = dimension['size'] if not dimension['UNLIMITED'] else None

            # If dimension is unlimited, size is None.
            if dimension['UNLIMITED']:
                size = None
            
            # Else if the idmension name is time, set equal to data index size.
            elif name=='time':
                size = flight_data.index.size

            # Else fall back on the size in the reference file.
            else:
                size = dimension['size']
                
            # Write the dimension into the file.
            ds.createDimension(name, size)

        # Loop over the variable list.
        for name in variable_list:
            
            # If this is a duplicate variable, skip it.
            if name in list(ds.variables):
                print(("   WARN: variable '{}' is duplicate variable!"
                       " Skipping.").format(name))
                continue

            # Select the reference variable row for the currernt variable.
            var_row = flight_variables[flight_variables['ICARTT_NAME']==name]

            # Get the variable reference JSON. If doesn't exist, skip variable.
            try:   
                # Get JSON variable reference file.
                var_map = var_row.VARIABLE_MAP.values[0]

                # Read from JSON to Python dicitonary.
                var_ref = read_netcdf_file_structure(var_map)

                 # If the reference is None, skip over this iterated variable.
                if var_ref is None:
                    print(("   WARN: variable '{}' has no JSON reference."
                           " Skipping.").format(name))
                    continue
            
            # If exception is raised, notify user and skip this variable.
            except Exception as e:
                # print(e)
                print(("   WARN: variable '{}' is not in variables table."
                       " Skipping.").format(name))
                continue
            
            # IF no exception is raised, write the variable.
            else:

                # Get the numpy data type from the reference datatype string.
                try:
                    datatype = getattr(np, var_ref['datatype'])
                except:
                    datatype = var_ref['datatype']

                try:
                    # Create a variable will a fill value.
                    x = ds.createVariable(
                        name, 
                        datatype, 
                        var_ref['dimensions'], 
                        fill_value=var_ref['attributes']['_FillValue'],
                        zlib=True)

                    # Delete the fill values from the attributes dict.
                    del var_ref['attributes']['_FillValue']

                except:
                    # Else create a variable without a fill value.
                    x = ds.createVariable(
                        name, 
                        datatype, 
                        var_ref['dimensions'], 
                        zlib=True)
                
                # Set the attirubtes.
                x.setncatts(var_ref['attributes'])
                
                # Add the data.
                if name=='FLIGHT':
                    x[:] = nc4.stringtochar(
                        np.array(list(flight_datestring), dtype='S1'))
                else:
                    x[:] = flight_data[name].to_numpy()


        ### SPECIAL ACT-AMERICA ROUTINE ---------------------------------------

        # This special routine exists because the science team wants ot merge
        # some custom metadata flags with the merge ICARTT data set.

        # Loop over the flag variables.
        for _flagv in [
            'Flight_flag',
            'Air_flag',
            'BL_FT_flag',
            'Maneuver_flag',
            'Maneuver_flagQC',
        ]:

            # Read the flag JSON config.
            flag_config = read_netcdf_file_structure(
                'references/actamerica/special/{}.json'.format(_flagv))

            # Replace serialized lists with tuples.
            for flag_key, flag_value in flag_config['attributes'].items():
                
                # If the type of the value is list, replace with tuple.
                if type(flag_value) is list:
                    flag_config['attributes'][flag_key] = tuple(flag_value)
            
            # If the flag variable name is in parsed ICARTT columns, add it.
            if _flagv in list(flight_data.columns):
                flag_data = flight_data[_flagv]
                flag_data = flag_data.fillna(-999999).astype(int).to_numpy()
                
                # If the variable is already in the output file, replace data.
                if flag_config['name'] in list(ds.variables):
                    ds.variables[flag_config['name']][:] = flag_data
                    continue
                
                # Else, create new variable.
                else:
                    
                    # Make a new variable and add the data.
                    x = ds.createVariable(
                        flag_config['name'], 
                        flag_config['datatype'], 
                        ('time'), 
                        fill_value=-999999, 
                        zlib=True)

                    # Set the attributes.
                    x.setncatts(flag_config['attributes'])

                    # Add the data.
                    x[:] = flag_data

        ### SPECIAL ACT-AMERICA ROUTINE ---------------------------------------




def _write_icartts_to_netcdf_files(DATA: dict):
    '''

    Parameters
    ----------

    Returns
    -------

    '''

    # Open the reference variables table.
    output_variables = pd.read_csv(DATA['VARIABLES'])

    # Open the reference structure json.
    output_structure = read_netcdf_file_structure(DATA['STRUCTURE'])

    # Make groupings of standard and multileg flights.
    DATA['FLIGHTS'] = _organize_standard_and_multileg_flights(DATA)

    print("\n### WRITING NETCDF FILES TO {}".format(DATA['DIR_OUTPUT']))
    counter = 0

    # Loop over the ICARTT files.
    for flight, icartt in DATA['FLIGHTS'].items():
        counter += 1
        print(" - [ {} / {} ] {} ".format(
            counter, len(DATA['FLIGHTS']), os.path.basename(flight)))

        # Parse the flight tables and add the data to the dictionary.
        flight_data = _parse_flight(icartt)

        # Get the full output netCDF path.
        flight_output = _construct_path_output(DATA['DIR_OUTPUT'], flight, "nc")

        # Get the flight resource files. If icartt is a list, use first.
        if type(icartt) is list:
            flight_resources = DATA['FLIGHT_RESOURCES'][icartt[0]]
        else:
            flight_resources = DATA['FLIGHT_RESOURCES'][flight]

        # Pass all of the information to write_netcdf.
        write_netcdf(
            output_file=flight_output, 
            output_structure=output_structure,
            flight_data=flight_data, 
            flight_variables=output_variables,
            flight_resources=flight_resources)





def _handle_input_configuration(input_config_path: str):
    '''
    '''

    print("\n### VALIDATE INPUT CONFIG (YAML): {}\n".format(input_config_path))

    try:
        with open(input_config_path, "r") as f:
            DATA = yaml.load(f, Loader=yaml.BaseLoader)
    except:
        with open("TEST.yml", "r") as f:
            DATA = yaml.load(f, Loader=yaml.BaseLoader)
    
    ### 1. Ensure ICARTT directory is valid.
    print((" ## 1. Input ICARTT directory: {DIR_ICARTT}").format(**DATA))
    if not os.path.isdir(DATA['DIR_ICARTT']):
        _exit_with_error("Input ICARTT directory is invalid.")

    ### 2. Ensure that it has ICARTT files in it.
    DATA['ICARTT_FILES'] = _crawl_directory(DATA['DIR_ICARTT'])
    if len(DATA['ICARTT_FILES'])==0:
        # If no icartts were found, exit and notify user.
        _exit_with_error("No ICARTT files found in the input directory.")
    else:
        # Else, inform on the number of ICARTT files
        print(" - Found [ {} ] ICARTTs.\n".format(len(DATA['ICARTT_FILES'])))

    ### 3. Ensure output directory exists. If not, create it.
    print((" ## 2. Output files directory: {DIR_OUTPUT}").format(**DATA))
    if not os.path.isdir(DATA['DIR_OUTPUT']):
        print(" - Output directory not found. Creating it now.".format(**DATA))
        try:
            os.mkdir(DATA['DIR_OUTPUT'])
        except Exception as e:
            _exit_with_error("Output directory cannot be created.")
            raise(e)

    return DATA




def main(DATA: dict):
    '''
    Description ...

    Parameters
    ----------

    Returns
    -------

    '''

    # Get absolute path to resources directory.
    DATA['DIR_RESOURCES'] = os.path.abspath(DATA['DIR_OUTPUT'])

    print("\n### WRITING RESOURCES IN DIRECTORY: {}".format(
        DATA['DIR_RESOURCES']))

    # Write new resource files.
    _write_resource_files(DATA)

    # Call the netCDF writer.
    _write_icartts_to_netcdf_files(DATA)

    return DATA




if __name__ == "__main__":

    # Handle the input configuration file.
    DATA = _handle_input_configuration(sys.argv[1])
    
    # Pass config to main.
    DATA = main(DATA)

    #print(DATA['variables_table'])
