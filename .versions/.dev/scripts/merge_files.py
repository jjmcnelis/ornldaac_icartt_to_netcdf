import pandas as pd
import glob
import os

# def read_multileg(icartt):
#     dfs = {}
#     for ict in icartt:
#         with open(ict, "r") as f:
#             headerrow = int(f.readlines()[0].split(",")[0])-1
#         dfs[ict] = pd.read_csv(ict, header=headerrow, delimiter=",")
#     df1, df2 = dfs.values()
#     df1 = df1.append(df2)
#     df1.columns = [c.strip() for c in list(df1.columns)]
#     return df1
    
# def read_singleleg(icartt):
#     with open(icartt, "r") as f:
#         headerrow = int(f.readlines()[0].split(",")[0])-1
#     df = pd.read_csv(icartt, header=headerrow, delimiter=",")
#     df.columns = [c.strip() for c in list(df.columns)]
#     return df

# def read_flight(flight_tuple):
#     icartt = flight_tuple[1]
#     flight_tuple2 = [flight_tuple[0], None]
#     if type(icartt) is list:
#         flight_tuple2[1] = read_multileg(icartt)
#     else:
#         flight_tuple2[1] = read_singleleg(icartt)
#     return flight_tuple2


files = glob.glob("/data/actamerica/ACTAMERICA_Merge/data/ict/ACTAMERICA-mrg*")
flights = [os.path.basename(f)[25:38] for f in files]
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

print(uflightsd)

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
    
# uflightsd3 = {}
# for k, v in uflightsd2.items():
#     try:
#         uflightsd3[k] = read_flight(v)
#     except:
#         uflightsd3[k] = None
#         print(k)
#         pass
    
# print(uflightsd3)