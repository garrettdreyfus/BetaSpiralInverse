from netCDF4 import Dataset
import numpy as np
from functools import partial

def searchBath(bathDataset,lat,lon):
    spacing = bathDataset.variables["spacing"][0]
    startlon = bathDataset.variables["x_range"][0]
    ret = bathDataset.variables["dimension"][:][0]
    z = bathDataset.variables["z"]
    i = int((90-lat)/spacing)
    j = int((lon-startlon)/spacing)
    if np.isnan(z[j+i*ret]):
        return 0
    else:
        return float(z[j+i*ret])

def bathBox(bathDataset,lat,lon,length=28):
    dlat = length/111.0
    dlon = length/(np.cos(np.deg2rad(lat))*111.0)
    dlon = min(dlon,170)
    botleftlat = lat-dlat
    botleftlon = lon-dlon
    toprightlat= min(lat+dlat,90)
    toprightlon = lon+dlon
    spacing = bathDataset.variables["spacing"][0]
    startlon = bathDataset.variables["x_range"][0]
    ret = bathDataset.variables["dimension"][:][0]
    z = bathDataset.variables["z"]
    if botleftlon < -180:
        botleftlon = (botleftlon+360)
        flip = np.arange(botleftlon,180,spacing)
        normal = np.arange(-180,toprightlon,spacing)
        gridlons = np.concatenate((normal,flip))
    elif toprightlon >180:
        toprightlon = toprightlon-360
        normal = np.arange(botleftlon,180,spacing)
        flip = np.arange(-180,toprightlon,spacing)
        gridlons = np.concatenate((normal,flip))
    else:
        gridlons = np.arange(botleftlon,toprightlon,spacing*30)

    gridlats = np.arange(botleftlat,toprightlat,spacing*6)
    idx = np.round(np.linspace(0, len(gridlons) - 1, 10)).astype(int)
    gridlons = gridlons[idx]
    idx = np.round(np.linspace(0, len(gridlats) - 1, 10)).astype(int)
    gridlats = gridlats[idx]
    latindexs = []
    lonindexs = []
    for l in gridlats:
        #i
        latindexs.append(int((90-l)/spacing))
    for l in gridlons:
        #j
        lonindexs.append(int((l-startlon)/spacing))
    depths = []
    for i in latindexs:
        for j in lonindexs:
            if not np.isnan(z[j+i*ret]):
                depths.append(float(z[j+i*ret]))
    return depths



d = Dataset("data/ver1_netcdf_geo.nc")
searchBath =partial(searchBath,d)
bathBox =partial(bathBox,d)
