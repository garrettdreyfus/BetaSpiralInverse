from netCDF4 import Dataset
import xarray as xr
import numpy as np
from functools import partial
import scipy
import scipy.io as sio
import pdb
## return the depths within a box
##this has a lot of finagling to deal with the wrap around
def bathBox(lat,lon,region,length=0.125,spacing=0.008333333,boxmethod="deg"):
    if boxmethod == "m":
        dlat = length/111.0
        dlon = length/(np.cos(np.deg2rad(lat))*111.0)
    if boxmethod == "deg":
        dlat = length
        dlon = length
    if lon >0: lon = -(180-lon) + -180
    dlon = min(dlon,170)
    botleftlat = lat-dlat
    botleftlon = lon-dlon
    toprightlat= min(lat+dlat,90)
    toprightlon = lon+dlon
    gridlons = np.linspace(botleftlon,toprightlon,29)
    gridlats = np.linspace(botleftlat,toprightlat,29)
    latindexs = []
    lonindexs = []
    depths = []
    for i in gridlats:
        for j in gridlons:
            z= region["searchBath"](i,j)
            if not np.isnan(z):
                depths.append(z)
    #pdb.set_trace()
    return depths

## given a surfaces object add depth information to each point
def addBathToSurfaces(surfaces,region):
    dumbcache = {}
    for k in surfaces.keys():
        surfaces[k]["data"]["z"] =  np.full(len(surfaces[k]["lons"]),np.nan)
        for l in range(len(surfaces[k]["lats"])):
            lat = surfaces[k]["lats"][l]
            lon = surfaces[k]["lons"][l]
            if (lat,lon) not in dumbcache.keys():
                dumbcache[(lat,lon)]=region["searchBath"](lat,lon)
            surfaces[k]["data"]["z"][l] = dumbcache[(lat,lon)]
    return surfaces

def addBathToSurface(surface,region):
    dumbcache = {}
    surface["data"]["z"] =  np.full(len(surface["lons"]),np.nan)
    for l in range(len(surface["lats"])):
        lat = surface["lats"][l]
        lon = surface["lons"][l]
        if (lat,lon) not in dumbcache.keys():
            dumbcache[(lat,lon)]=region["searchBath"](lat,lon)
        surface["data"]["z"][l] = dumbcache[(lat,lon)]
    return surface

#simplify functions so that they already have the bathymetry file loaded

