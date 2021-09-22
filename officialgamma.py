import matlab.engine
import nstools
import time
import matplotlib.pyplot as plt

def addOfficialGamma(surfaces):
    if "pottemp" not in surfaces[list(surfaces.keys())[0]]:
        nstools.addOldUnits(surfaces)
    eng = matlab.engine.start_matlab()
    for k in surfaces.keys():
        start = time.time()
        sp = matlab.double(list(surfaces[k]["data"]["psu"]))
        t = matlab.double(list(surfaces[k]["data"]["pottemp"]))
        p = matlab.double(list(surfaces[k]["data"]["pres"]))
        lons = matlab.double(list(surfaces[k]["lons"]))
        lats = matlab.double(list(surfaces[k]["lats"]))
        gamma = eng.eos80_legacy_gamma_n(sp,t,p,lons,lats)
        end = time.time()
        surfaces[k]["data"]["gamma"]=gamma[0]
    return surfaces

