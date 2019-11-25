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
from . import MODULE_PATH
from ._utils import *
print(MODULE_PATH)
