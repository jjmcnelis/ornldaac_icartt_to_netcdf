#!/usr/bin/env python3

import os
import glob
import pandas as pd

files = glob.glob("/data/actamerica/ACTAMERICA_Merge/data/ict/ACTAMERICA-mrg*")
flights = [os.path.basename(f) for f in files]
uflights = list(set(flights))
uflightsd = {}

for flt0 in uflights:
    t = []
    for flt1 in files:
        if flt0 in flt1:
            t.append(flt1)
    t = sorted(t)
    if len(t)>1:
        if (("_L1_" in t[-2]) & ("_L2_" in t[-1])):
            uflightsd[flt0] = t
        else:
            uflightsd[flt0] = t[-1]
    else:
        uflightsd[flt0] = t[0]

print(list(uflightsd.keys())[0])
print(uflightsd[list(uflightsd.keys())[0]])


# uflightsd2 = {}
# for k, v in uflightsd.items():
#     file_ = None
#     for f in ncfiles:
#         fcomp = os.path.basename(f).split("_")
#         aircraft = fcomp[0].split("-")[-1].upper()
#         date = fcomp[-2]  
#         if ((aircraft in k) & (date in k)):
#             file_ = f
#     uflightsd2[k] = (file_, v)

# print(uflightsd2)