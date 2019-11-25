#def test(config: str=TESTS_PATH + "CONFIG.yml"):
def test(config: str=TESTS_CONFIG): 
    '''
    Description ...

    Parameters
    ----------

    Returns
    -------

    '''

    from io import StringIO

    CONFIG = yaml.load(StringIO("""
    ### TEST CONFIG
    
    # Relative or absolute path to read directory for ICARTT file(s).
    DIR_ICARTT: /data/actamerica/ACTAMERICA_Merge/data/ict/
    
    # Relative or absolute path to an output directory. 
    PATH_OUTPUT: .tests/test1/
    
    # Path to reference configurations. Will be created if it doesn't exist.
    RESOURCES_PATH: .resources/
    
    # Relative or absolute path to the variable translation table.
    CSV_VARIABLES: .resources/actamerica/translation_table.csv
    
    # netCDF structure for output files.
    JSON_NETCDF: .resources/actamerica/JSON_NETCDF.yml
    
    # Merge multileg flights (have _Lx at tail of filename) end to end.
    JOIN_MULTILEG: true
    """))

    # Check and dump config to test dir.
    RESULT = main(CONFIG)

    print(RESULT)

    return RESULT

    # # (Make the test ICARTT file configurable).
    # test_file = "tests/ACTAMERICA-mrgPFP-c130_merge_20180520_R0.ict"

    # # Parse the header of the test ICARTT.
    # test_file_variables, test_file_metadata = parse_icartt_header(test_file)

    # # Get the full path and the basename of the ICARTT without the extension.
    # test_file_path = os.path.dirname(test_file)
    # test_file_name = os.path.splitext(os.path.basename(test_file))[0]

    # # Join the output path.
    # test_file_output = os.path.join(test_file_path, "{}".format(test_file_name))

    # # Write the test variable table to a CSV in tests/
    # test_file_variables.to_csv(test_file_output + ".csv")

    # # Also dump the header data (excluding table) to a JSON file.
    # with open(test_file_output + ".json", "w") as f:
    #     f.write(json.dumps(test_file_metadata, indent=2))