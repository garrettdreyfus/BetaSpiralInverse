import glob 
from netCDF4 import Dataset 
import numpy as np
from profile import Profile
from functools import partial
import sys
import xarray as xr
import pdb
from progress.bar import Bar


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

def createMesh(xvals,yvals,coord="xy",spacingscale=1):
    if coord =="latlon":
        grd = np.meshgrid(\
                np.linspace(-75,-45,35),np.linspace(20,44,30))
        for i in range(grd[0].shape[0]):
            for j in range(grd[0].shape[1]):
                x,y = singleXY((grd[0][i][j],grd[1][i][j]))
                grd[0][i][j] = x
                grd[1][i][j] = y
        return grd
    else:
        print("This grid type is not supported", sys.exc_info()[0])



def extractArgoProfiles(ncfolder): 
    profiles= []
    for f in Bar("file:" ).iter(glob.glob(ncfolder+"/**/*.nc",recursive=True)):
        ncdf = Dataset(f, 'r')  # Dataset is the class behavior to open the file
        #pdb.set_trace()
        #print("yup")
        #print(ncdf.dimensions["N_PROF"].size)
        for prof in range(ncdf.dimensions["N_PROF"].size):
            profiledata = {}
            profiledata["lon"] = ncdf.variables["LONGITUDE"][prof]
            profiledata["lat"] = ncdf.variables["LATITUDE"][prof]
            if geoFilter(profiledata["lon"],profiledata["lat"]) and prof%5 ==0:
                #profiledata["cruise"] = ncdf.WOCE_ID
                #profiledata["station"] = ncdf.STATION_NUMBER
                profiledata["time"] = ncdf.variables["JULD"][prof]
                #profiledata["cruise"] = ncdf.WOCE_ID
                #profiledata["station"] = ncdf.STATION_NUMBER
                if "TEMP_ADJUSTED" in ncdf.variables.keys():
                    profiledata["temp"] = np.asarray(ncdf.variables["TEMP_ADJUSTED"][prof][:])
                if "PSAL_ADJUSTED" in ncdf.variables.keys():
                    profiledata["sal"] = np.asarray(ncdf.variables["PSAL_ADJUSTED"][prof][:])
                if "PRES_ADJUSTED" in ncdf.variables.keys():
                    profiledata["pres"] = np.asarray(ncdf.variables["PRES_ADJUSTED"][prof][:])
                if len(profiledata["pres"])>4 and max(profiledata["pres"])>1500\
                        and geoFilter(profiledata["lon"],profiledata["lat"]):
                        if {"sal","temp","pres","lat","lon"}.issubset(profiledata.keys())\
                            and abs(max(profiledata["pres"])-min(profiledata["pres"])) > 100\
                            and 99999 not in profiledata["pres"] and 99999 not in profiledata["temp"]\
                            and 99999 not in profiledata["sal"]:
                            
                            eyed=idgenerator()
                            prof=Profile(eyed,profiledata)
                            profiles.append(prof)
                            #print(profiledata)
        del ncdf
    return profiles

def geoFilter(lon,lat):
    latinrange = (lat>20 and lat <44)
    loninrange = (lon>-75 and lon < -35)
    return (latinrange and loninrange)

mapbounds = {"lllon":-81,"urlon":-35,"lllat":20,"urlat":44,\
        "lat_0":31,"lon_0":-50}



