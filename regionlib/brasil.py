import glob 
from netCDF4 import Dataset 
import numpy as np
from profile import Profile
from functools import partial
import sys
import xarray as xr
import pdb

def geoFilter(lon,lat):
    latinrange = (lat<0 and lat >-31)
    loninrange = (lon>-69 and lon < -12)
    return (latinrange and loninrange)

#generate a unique id
def idgenerator():
    if not hasattr(idgenerator,"id"):
        idgenerator.id=0
    idgenerator.id = idgenerator.id+1
    return idgenerator.id


##lat lon to x y
def singleXY(coord):
    theta = np.deg2rad(coord[0])
    r = ((90-coord[1]) *111*1000)
    x = (r*np.cos(theta))
    y = (r*np.sin(theta))
    return x,y


def extractProfiles(ncfolder): 
    profiles= []
    for f in glob.glob(ncfolder+"/*.nc"):
        ncdf = Dataset(f, 'r')  # Dataset is the class behavior to open the file
        profiledata = {}
        profiledata["lon"] = ncdf.variables["longitude"][0]
        profiledata["lat"] = ncdf.variables["latitude"][0]
        profiledata["cruise"] = ncdf.WOCE_ID
        profiledata["station"] = ncdf.STATION_NUMBER
        profiledata["time"] = ncdf.variables["time"][0]
        profiledata["cruise"] = ncdf.WOCE_ID
        profiledata["station"] = ncdf.STATION_NUMBER
        if "temperature" in ncdf.variables.keys():
            profiledata["temp"] = np.asarray(ncdf.variables["temperature"][:])
        if "salinity" in ncdf.variables.keys():
            profiledata["sal"] = np.asarray(ncdf.variables["salinity"][:])
        if "pressure" in ncdf.variables.keys():
            profiledata["pres"] = np.asarray(ncdf.variables["pressure"][:])
        if "CTDTMP" in ncdf.variables.keys():
            profiledata["temp"] = np.asarray(ncdf.variables["CTDTMP"][:])
        if "CTDSAL" in ncdf.variables.keys():
            profiledata["sal"] = np.asarray(ncdf.variables["CTDSAL"][:])
        if len(profiledata["pres"])>4 and max(profiledata["pres"])>1500\
                and geoFilter(profiledata["lon"],profiledata["lat"]):
                eyed=idgenerator()
                prof=Profile(eyed,profiledata)
                profiles.append(prof)
    return profiles

def createMesh(n,xvals,yvals,coord="xy"):
    if coord == "xy":
        x1,y1 = singleXY((-68,-54))
        x2,y2 = singleXY((-21,-7))
        xmin = min(x1,x2)
        xmax = max(x1,x2)
        ymin = min(y1,y2)
        ymax = max(y1,y2)
        return np.meshgrid(np.linspace(xmin,xmax,n), np.linspace(ymin,ymax,n),indexing="xy")
    if coord =="latlon":
        grd = np.meshgrid(\
                np.linspace(-71,-12,45),np.linspace(-31,-2,30))
        for i in range(grd[0].shape[0]):
            for j in range(grd[0].shape[1]):
                x,y = singleXY((grd[0][i][j],grd[1][i][j]))
                grd[0][i][j] = x
                grd[1][i][j] = y
        return grd
    else:
        print("This grid type is not supported", sys.exc_info()[0])






