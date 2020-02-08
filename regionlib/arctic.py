from netCDF4 import Dataset
from functools import partial
import xarray as xr

def createMesh(n,xvals,yvals,coord="xy"):
    xmin= -1793163
    xmax = 971927
    ymin = -1455096
    ymax = 1200385
    return np.meshgrid(np.linspace(xmin,xmax,n), np.linspace(ymin,ymax,n),indexing="xy")

def geoFilter(lon,lat): return True


#Return depth at certain lat lon
def searchBath(bathDataset,lat,lon):
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


arcticSearchBath =partial(searchBath,Dataset("data/ver1_netcdf_geo.nc"))

mapbounds = {"lllon":-136,"urlon":78,"lllat":55,"urlat":63,\
        "lat_0":90,"lon_0":-60}


