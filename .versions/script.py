        

    # # Get the list of unique variable names.
    # variable_names = sorted(list(set(list(variable_reference.keys()))))

    # ### Update translation table (key table ict -> nc).
    
    # with open(os.path.join(resources_dir + "translate.csv"), "w") as f:
    #     f.write("\n".join(variable_names))

    # # If there is no netcdf_variables/ subdir, create it.
    # netcdf_variable_dir = os.path.join(resources_dir, "netcdf_variables/")
    # if not os.path.isdir(netcdf_variable_dir):
    #     os.mkdir(netcdf_variable_dir)

    # ### Write individual variable jsons.
    # #_write_variable_jsons(, netcdf_variable_dir)

    # # If the file already exists read with pandas.
    # #if os.path.isfile(translate_csv):
    # #    translate_tbl = pd.read_csv(translate_csv)
    # #else:
    #     #pd.DataFrame([(v, , , ) for v in variable_names])