# netcdf reference files

Write the JSON representation of the netCDF structure to a file.

```python
from ornldaac_icartt_to_netcdf import *
write_netcdf_file_structure("references/actamerica/netcdf/ACTAMERICA-mrg05-b200_merge_20160711_R0.nc")
write_netcdf_file_structure("references/actamerica/netcdf/ACTAMERICA-mrg05-c130_merge_20160718_R3.nc")
```

Write a reference JSON for every netCDF variable:

```python
from ornldaac_icartt_to_netcdf import *
write_variable_jsons(
    "references/actamerica/netcdf/ACTAMERICA-mrg05-b200_merge_20160711_R0.json", 
    "references/actamerica/variables/B200/")
write_variable_jsons(
    "references/actamerica/netcdf/ACTAMERICA-mrg05-c130_merge_20160718_R3.json", 
    "references/actamerica/variables/C130/")

```