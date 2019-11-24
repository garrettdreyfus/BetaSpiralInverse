from netCDF4 import Dataset
import xarray as xr
import numpy as np
from functools import partial
import scipy
import scipy.io as sio
import pdb

#Return depth at certain lat lon
def arcticSearchBath(bathDataset,lat,lon):
    spacing = bathDataset.variables["spacing"][0]
    startlon = bathDataset.variables["x_range"][0]
    ret = bathDataset.variables["dimension"][:][0]
    z = bathDataset.variables["z"]
    if np.isnan(lat) or np.isnan(lon):
        return np.nan
    i = int((90-lat)/spacing)
    j = int((lon-startlon)/spacing)
    if np.isnan(z[j+i*ret]):
        return 0
    else:
        return float(z[j+i*ret])

#Return depth at certain lat lon
def nepbSearchBath(bathDataset,lat,lon):
    if np.isnan(lon) or np.isnan(lat):
        return np.nan
    f = bathDataset.sel(lon=lon, lat=lat, method='nearest')
    if ~np.isnan(f):
        return float(-f.elevation.values)
    else:
        return f

def nepbMatlabSearchBath(bathDataset,lat,lon,includeindex=False):
    lati = np.abs(bathDataset["lat_SRTM"]-lat).argmin()
    loni = np.abs(bathDataset["lon_SRTM"]-lon).argmin()
    #print("---")
    #print(lati,loni)
    if includeindex:
        return bathDataset["z_SRTM"][lati][loni], lati, loni
    else:
        return bathDataset["z_SRTM"][lati][loni]



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
            z= searchBath(i,j,region)
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
                dumbcache[(lat,lon)]=searchBath(lat,lon,region)
            surfaces[k]["data"]["z"][l] = dumbcache[(lat,lon)]
    return surfaces

def addBathToSurface(surface,region):
    dumbcache = {}
    surface["data"]["z"] =  np.full(len(surface["lons"]),np.nan)
    for l in range(len(surface["lats"])):
        lat = surface["lats"][l]
        lon = surface["lons"][l]
        if (lat,lon) not in dumbcache.keys():
            dumbcache[(lat,lon)]=searchBath(lat,lon,region)
        surface["data"]["z"][l] = dumbcache[(lat,lon)]
    return surface

def searchBath(lat,lon,region):
    regionToFunction = {"arctic":arcticSearchBath,\
            "nepb":nepbSearchBath,"nepbmatlab":nepbMatlabSearchBath}
    return regionToFunction[region](lat,lon)


#simplify functions so that they already have the bathymetry file loaded
arcticSearchBath =partial(arcticSearchBath,Dataset("data/ver1_netcdf_geo.nc"))
nepbSearchBath =partial(nepbSearchBath,xr.open_dataset('data/nepbbath.nc',decode_times=False))
nepbMatlabSearchBath =partial(nepbMatlabSearchBath,sio.loadmat("data/Run0.new.mat"))

