from netCDF4 import Dataset
import xarray as xr
import numpy as np
from functools import partial
import scipy
import scipy.io as sio
import pdb
## given a surfaces object add depth information to each point
def addBathToSurfaces(surfaces,region):
    dumbcache = {}
    for k in surfaces.keys():
        surfaces[k]["data"]["z"] =  np.full(len(surfaces[k]["lons"]),np.nan)
        for l in range(len(surfaces[k]["lats"])):
            lat = surfaces[k]["maplats"][l]
            lon = surfaces[k]["lons"][l]
            if (lat,lon) not in dumbcache.keys():
                dumbcache[(lat,lon)]=searchBath(lat,lon)
            surfaces[k]["data"]["z"][l] = dumbcache[(lat,lon)]
    return surfaces

def addBathToSurface(surface,region):
    dumbcache = {}
    surface["data"]["z"] =  np.full(len(surface["lons"]),np.nan)
    for l in range(len(surface["lats"])):
        lat = surface["maplats"][l]
        lon = surface["lons"][l]
        if (lat,lon) not in dumbcache.keys():
            dumbcache[(lat,lon)]=searchBath(lat,lon)
        surface["data"]["z"][l] = dumbcache[(lat,lon)]
    return surface

#Return depth at certain lat lon
def searchBath(bathDataset,lat,lon):
    #lons = bathDataset.variables["lon"]
    #lats = bathDataset.variables["lat"]
    #loni = int(86400*(lon+180)/360.0)
    #lati = int(43200*(lat+90)/180.0)
    if np.isnan(lon) or np.isnan(lat):
        return np.nan
    f = bathDataset.sel(lon=lon,lat=lat,method="nearest")
    if ~np.isnan(f):
        return float(f.z.values)
    else:
        return f

def bathBox(bathDataset,lat,lon):
    f = bathDataset.sel(lon=slice(lon-0.1,lon+0.1),lat=slice(lat-0.1,lat+0.1))
    if ~np.isnan(f):
        f=f.stack(g=("lon","lat"))
        return np.ndarray.flatten(f.lon.values),np.ndarray.flatten(f.lat.values),np.ndarray.flatten(f.z.values)
    else:
        return f

searchBath = partial(searchBath,xr.open_dataset("data/SRTM.nc"))
bathBox = partial(bathBox,xr.open_dataset("data/SRTM.nc"))


#simplify functions so that they already have the bathymetry file loaded

