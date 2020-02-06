from functools import partial
import aleutianline as al
import xarray as xr
import scipy.io as sio

def hautalaGrid():
    grd = np.meshgrid(\
            np.concatenate((np.linspace(170,178,5),np.linspace(-180,-122,30))),

            np.linspace(18,58,21))
    for i in range(grd[0].shape[0]):
        for j in range(grd[0].shape[1]):
            x,y = singleXY((grd[0][i][j],grd[1][i][j]))
            grd[0][i][j] = x
            grd[1][i][j] = y
    return grd


def createMesh(n,xvals,yvals,coord="xy"):
    if coord=="latlon":
        return hautalaGrid()
    elif coord == "xy":
        x1,y1 = singleXY((-144,3))
        x2,y2 = singleXY((155,63))
        xmin = min(x1,x2)
        xmax = max(x1,x2)
        ymin = min(y1,y2)
        ymax = max(y1,y2)
        return np.meshgrid(np.linspace(xmin,xmax,n), np.linspace(ymin,ymax,n),indexing="xy")
    else:
        print("This grid type is not supported", sys.exc_info()[0])

def geoFilter(lon,lat):
    if lon > 0: lon=-180-(abs(lon)-180)
    if lon>=np.min(al.lon) and lon<=np.max(al.lon):
        if lat-np.interp(lon,al.lon,al.lat)>-2:
            #IN THE ALEUTIANS
            return False
    return True

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


nepbSearchBath =partial(nepbSearchBath,xr.open_dataset('data/nepbbath.nc',decode_times=False))
nepbMatlabSearchBath =partial(nepbMatlabSearchBath,sio.loadmat("data/Run0.new.mat"))


