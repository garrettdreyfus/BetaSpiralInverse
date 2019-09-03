from netCDF4 import Dataset
import numpy as np
from functools import partial


#Return depth at certain lat lon
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

## return the depths within a box
##this has a lot of finagling to deal with the wrap around
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

## given a surfaces object add depth information to each point
def addBathToSurfaces(surfaces):
    dumbcache = {}
    for k in surfaces.keys():
        surfaces[k]["data"]["z"] =  np.full(len(surfaces[k]["lons"]),np.nan)
        for l in range(len(surfaces[k]["lats"])):
            lat = surfaces[k]["lats"][l]
            lon = surfaces[k]["lons"][l]
            if (lat,lon) not in dumbcache.keys():
                dumbcache[(lat,lon)]=searchBath(lat,lon)
            surfaces[k]["data"]["z"][l] = dumbcache[(lat,lon)]
    return surfaces

def addBathToSurface(surface):
    dumbcache = {}
    surface["data"]["z"] =  np.full(len(surface["lons"]),np.nan)
    for l in range(len(surface["lats"])):
        lat = surface["lats"][l]
        lon = surface["lons"][l]
        if (lat,lon) not in dumbcache.keys():
            dumbcache[(lat,lon)]=searchBath(lat,lon)
        surface["data"]["z"][l] = dumbcache[(lat,lon)]
    return surface
            


#simplify functions so that they already have the bathymetry file loaded
d = Dataset("data/ver1_netcdf_geo.nc")
searchBath =partial(searchBath,d)
bathBox =partial(bathBox,d)
