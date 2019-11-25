# EVS support services: ICARTT to netCDF

## what

This repo makes public the ORNL DAAC workflow to translate ICARTT data files from the EVS missions into archive ready netCDF data sets. 

This is a total rewrite as of 2019-11-22, so expect some bugs. Please email mcnelisjj@ornl.gov if you find any.

## why

The ORNL DAAC requires netCDF data sets to be self-describing according to the CF conventions (as much as is reasonable) because we feel it increases the longevity of the data and it provides conveniences to data users in ways that the people reading this readme will I'm sure understand.

Public code allows investigators in active projects to reformat their data as needed.

## how

The code is described line by line. To summarize:

### 1. read and validate config file

Example: *`ACTAMERICA_B200.yml`*

```yaml
DIR_ICARTT: /data/actamerica/ACTAMERICA_Merge/data/ict/b200/
DIR_OUTPUT: outputs/ACTAMERICA_Merge_B200/
VARIABLES: references/actamerica/VARIABLES_B200.csv
STRUCTURE: references/actamerica/STRUCTURE_B200.json
RESOURCES:
  ICARTT_HEADERS: icartt_headers/
  ICARTT_VARIABLES: icartt_variables/
```

* `DIR_ICARTT`: */data/actamerica/ACTAMERICA_Merge/data/ict/b200/*

Path to the directory containing the input ICARTT files.

* `DIR_OUTPUT`: *outputs/ACTAMERICA_Merge_B200/*

Path to the output directory where resource files (ICARTT header and variable metadata) will be written along with the converted netCDF files.

* `VARIABLES`: *references/actamerica/VARIABLES_B200.csv*

Variable table maps the input ICARTT variable names (`ICARTT_NAME`) to the output netCDF variable names (`NETCDF_NAME`) and, most importantly, to the JSON reference file (`VARIABLE_MAP`). This table needs to be carefully maintained. The script loops over all of the variables in the ICARTT file, compares them against this table, and writes the data to the output netCDF using the `NETCDF_NAME`, and the dimensions/attributes given in the variable reference file (which are located here: *`references/actamerica/variables/[AIRCRAFT]/[VARIABLE].json`*).

* `STRUCTURE`: *references/actamerica/STRUCTURE_B200.json*

Path to the JSON representation of the output netCDF structure (the ORNL DAAC copies are located: *`references/actamerica/STRUCTURE_[AIRCRAFT].json`*). NOTE: The only information sourced from this file are the *dimensions* and *global attributes* that are written to the output netCDF files.

* `RESOURCES`: Don't worry about these. Every time you run the script, the input ICARTTs will be parsed and their headers and variable metadata will be written to "resource" files in these subdirectories in `DIR_OUTPUT`. These are compared against the reference files as the netCDFs are being written and warnings are printed when a match is not found.

### 2. loop over input ICARTT files and write resource files.

The script loops over the ICARTTs in the directory specified by  `DIR_ICARTT`, parses the header and variable metadata, and writes them to JSONs and CSVs in `RESOURCES` subdirectories within the `DIR_OUTPUT` directory.

### 3. validate input variables against reference variables as they are written to netCDF

The script will then attempt to match every unique variable name from the input ICARTT files to a reference JSON file (`VARIABLE_MAP` in `VARIABLES` CSV file) and output variable name (`NETCDF_NAME` in `VARIABLES` CSV file). Any variables that are not represented in the table will not be written to the output netCDF files.

Looping over the ICARTTs, multileg flights are paired and concatenated, and the ICARTT data are translated into netCDF files.

Each variable in the ICARTT is written to the output netCDF with the attributes given in the matching variable reference JSON; e.g. *references/actamerica/variables/B200/AircraftSunAzimuth.json*:

```json
{
  "datatype": "float64",
  "dimensions": [
    "time"
  ],
  "attributes": {
    "_FillValue": -999999.0,
    "units": "degree",
    "long_name": "aircraft sun azimuth",
    "standard_name": "platform_azimuth_angle",
    "coordinates": "time LATITUDE LONGITUDE GPS_ALT"
  }
}
```
The output variable is assigned the type specified in `datatype` field, the dimensions given in the `dimensions` field, and the attributes given in the `attributes` field. 

Any variables that cannot be matched to one of these JSON files WILL NOT be written to the output netCDF. Any variables that are skipped will be printed to stdout with a message to try to give you some clues:

```shell
 - [ 85 / 85 ] ACTAMERICA-mrg05-c130_merge_20160803_R3.ict 
   WARN: 'INDEX' has no reference. Skip
   WARN: 'WNS' not in variables table. Skip
   WARN: 'WND' not in variables table. Skip
   WARN: 'H2O_MixingRatio_Nav' has no reference. Skip
   WARN: 'Altitude-AGL_GoogleMaps' not in variables table. Skip
   WARN: 'GroundHeight-AMSL_CPL' has no reference. Skip
   WARN: 'MLH-AMSL_CPL' has no reference. Skip
   WARN: 'O3_MoleFraction' not in variables table. Skip
   WARN: 'DICL_MoleFraction_PFP' not in variables table. Skip
   WARN: 'C2H6_MoleFraction_PFP' not in variables table. Skip
```

I still need to resolve the missing reference files for the variables listed above for the C130 files.