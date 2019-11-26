# EVS ICARTT to netCDF4

## quick start

**Requirements:** Python 3

* `numpy`: https://numpy.org/
* `pandas`: https://pandas.pydata.org/
* `netCDF4`: https://unidata.github.io/netcdf4-python/netCDF4/index.html

**Create/edit configuration YAML file** (e.g. [`ACTAMERICA_B200.yml`](ACTAMERICA_B200.yml)) 

Point it at the ICARTT input directory ([`DIR_ICARTT`](inputs/ACTAMERICA_Merge/B200/)) and the output directory ([`DIR_OUTPUT`](outputs/ACTAMERICA_Merge_B200/)):

```yaml
# Path to some input ICARTT files. Will crawl recursively.
DIR_ICARTT: inputs/ACTAMERICA_Merge/B200/

# Path to write output resource files and netCDFs.
DIR_OUTPUT: outputs/ACTAMERICA_Merge_B200/

# Path to variable reference table, maps ICARTT<>netCDF<>CF.json names.
VARIABLES: references/actamerica/VARIABLES_B200.csv

# Path to JSON representation of netCDF output structure.
STRUCTURE: references/actamerica/STRUCTURE_B200.json

# Paths within resources subdirectory of DIR_OUTPUT (no need to change).
RESOURCES:

    # Parsed header + helpful info from each ICARTT file are written here.
    ICARTT_HEADERS: icartt_headers/

    # Variable reference parsed from ICARTT headers are written here.
    ICARTT_VARIABLES: icartt_variables/
```

Then pass the configuration file as the only argument to the `ornldaac_icartt_to_netcdf` module (to [`ornldaac_icartt_to_netcdf/__main__.py`](ornldaac_icartt_to_netcdf/__main__.py)):

```python
python -m ornldaac_icartt_to_netcdf [CONFIG].yml
```

As long as you have place the ICARTTs in an accessible folder and the full path is given in the `CONFIG.yml` file, everything should be good. I will refine this soon.

**to do**

* fill value replacement -- `OUTPUT_FILL` column has no affect for now
* user function -- `USER_FUNC` will add soon

## overview

### what

This repo makes public the ORNL DAAC workflow to translate ICARTT data files from the EVS missions into archive ready netCDF data sets. 

This is a total rewrite as of 2019-11-22, so expect some bugs. Please email mcnelisjj@ornl.gov if you find any.

### why

The ORNL DAAC requires netCDF data sets to be self-describing according to the CF conventions (as much as is reasonable) because we feel it increases the longevity of the data and it provides conveniences to data users in ways that the people reading this readme will I'm sure understand.

Public code allows investigators in active projects to reformat their data as needed.

### how

The code is described line by line in the `_utils.py` and `__main__.py` scripts. To summarize:

#### 1. read and validate config file

Example: *`ACTAMERICA_B200.yml`*

```yaml
DIR_ICARTT: inputs/ACTAMERICA_Merge/B200/
DIR_OUTPUT: outputs/ACTAMERICA_Merge_B200/
VARIABLES: references/actamerica/VARIABLES_B200.csv
STRUCTURE: references/actamerica/STRUCTURE_B200.json
RESOURCES:
  ICARTT_HEADERS: icartt_headers/
  ICARTT_VARIABLES: icartt_variables/
```

* `DIR_ICARTT`: *inputs/ACTAMERICA_Merge/B200/*

Path to the directory containing the input ICARTT files.

* `DIR_OUTPUT`: *outputs/ACTAMERICA_Merge_B200/*

Path to the output directory where resource files (ICARTT header and variable metadata) will be written along with the converted netCDF files.

* `VARIABLES`: *references/actamerica/VARIABLES_B200.csv*

Variable table maps the input ICARTT variable names (`ICARTT_NAME`) to the output netCDF variable names (`OUTPUT_NAME`) and, most importantly, to the JSON reference file (`VARIABLE_MAP`). This table needs to be carefully maintained. The script loops over all of the variables in the ICARTT file, compares them against this table, and writes the data to the output netCDF using the `OUTPUT_NAME`, and the dimensions/attributes given in the variable reference file (which are located here: *`references/actamerica/variables/[AIRCRAFT]/[VARIABLE].json`*).

* `STRUCTURE`: *references/actamerica/STRUCTURE_B200.json*

Path to the JSON representation of the output netCDF structure (the ORNL DAAC copies are located: *`references/actamerica/STRUCTURE_[AIRCRAFT].json`*). NOTE: The only information sourced from this file are the *dimensions* and *global attributes* that are written to the output netCDF files.

* `RESOURCES`: Don't worry about these. Every time you run the script, the input ICARTTs will be parsed and their headers and variable metadata will be written to "resource" files in these subdirectories in `DIR_OUTPUT`. These are compared against the reference files as the netCDFs are being written and warnings are printed when a match is not found.

#### 2. loop over input ICARTT files and write resource files.

The script loops over the ICARTTs in the directory specified by  `DIR_ICARTT`, parses the header and variable metadata, and writes them to JSONs and CSVs in `RESOURCES` subdirectories within the `DIR_OUTPUT` directory.

#### 3. validate input variables against reference variables as they are written to netCDF

The script will then attempt to match every unique variable name from the input ICARTT files to a reference JSON file (`VARIABLE_MAP`) and output variable name (`OUTPUT_NAME`). Any variables that are not represented in the table will not be written to the output netCDF files.

Looping over the ICARTTs, multileg flights are paired and concatenated, and the ICARTT data are translated into netCDF files.

Each variable in the ICARTT is written to the output netCDF with the attributes given in the matching variable reference JSON; e.g. *[references/actamerica/variables/B200/AircraftSunAzimuth.json](references/actamerica/variables/B200/AircraftSunAzimuth.json)*:

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

Any variables that cannot be matched to one of these JSON files WILL NOT be written to the output netCDF. Any variables that are skipped will be printed to stdout to clue you in. `INDEX` is deliberately excluded:

```shell
 - [ 85 / 85 ] ACTAMERICA-mrg05-c130_merge_20160803_R3.ict 
   WARN: 'INDEX' has no reference. Skip
```
