from netCDF4 import Dataset 
import numpy as np
from profile import Profile
from functools import partial
import sys,pdb,julian,glob,csv
import xarray as xr
from progress.bar import Bar
import matplotlib.pyplot as plt
import datetime
import json

def geoFilter(lon,lat):
    latinrange = (lat<-10 and lat >-46)
    loninrange = (lon>-79 and lon < 0)
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
    return coord


def extractWoceProfiles(ncfolder): 
    profiles= []
    for f in Bar("WOCE").iter(glob.glob(ncfolder+"/*.nc") + glob.glob(ncfolder+"/**/*.nc",recursive=True)):
        ncdf = Dataset(f, 'r')  # Dataset is the class behavior to open the file
        profiledata = {}
        profiledata["lon"] = ncdf.variables["longitude"][0]
        profiledata["lat"] = ncdf.variables["latitude"][0]
        profiledata["cruise"] = ncdf.WOCE_ID
        profiledata["station"] = ncdf.STATION_NUMBER
        profiledata["time"] = julian.to_jd(datetime.datetime(1980,1,1,0) + datetime.timedelta(minutes=int(ncdf.variables["time"][0])) )
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
                prof=Profile(eyed,profiledata,"insitu","practical")
                profiles.append(prof)
    return profiles


def extractBodcProfiles(ncfolder): 
    profiles= []
    maxs= []
    fieldspass = 0
    locationpass = 0
    depthpass = 0
    for f in Bar(" cdf file: ").iter(glob.glob(ncfolder+"/*.nc")):
        ncdf = Dataset(f, 'r')  # Dataset is the class behavior to open the file
        profiledata = {}
        profiledata["lon"] = ncdf.variables["LONGITUDE"][0]
        profiledata["lat"] = ncdf.variables["LATITUDE"][0]
        #profiledata["cruise"] = ncdf.SDN_CRUISE
        #profiledata["station"] = ncdf.SDN_STATION
        profiledata["time"] = ncdf.variables["TIME"][0]
        #profiledata["station"] = ncdf.STATION_NUMBER
        if "POTMCV01" in ncdf.variables.keys():
            profiledata["temp"] = np.asarray(ncdf.variables["POTMCV01"][:][0])
        if "PSALCC01" in ncdf.variables.keys():
            profiledata["sal"] = np.asarray(ncdf.variables["PSALCC01"][:][0])
        if "PRES" in ncdf.variables.keys():
            profiledata["pres"] = np.asarray(ncdf.variables["PRES"][:][0])
        maxs.append(np.max(profiledata["pres"]))
        if np.max(profiledata["pres"]) >1000:
            depthpass+=1
        if len(profiledata["pres"])>4 and np.max(profiledata["pres"])>1000\
                and geoFilter(profiledata["lon"],profiledata["lat"]):
                locationpass+=1
                eyed=idgenerator()
                try:
                    prof=Profile(eyed,profiledata,"potential","practical")
                    profiles.append(prof)
                    fieldspass+=1
                except ValueError as e:
                    print(profiledata.keys())
                    for k in ncdf.variables.keys():
                        if "SAL" in k:
                            print(k)
                    print(e)
        #else:
            #print((profiledata["lon"],profiledata["lat"]),geoFilter(profiledata["lon"],profiledata["lat"]))
            #print(profiledata["pres"])
            #plt.plot(profiledata["pres"])
            #plt.show()
    print("fieldspass: ",fieldspass)
    print("depthpass: ",depthpass)
    print("locationpass: ",locationpass)
    print("above 1000: ",np.count_nonzero(np.asarray(maxs)>1000))
    return profiles


def extractArgoProfiles(ncfolder): 
    profiles= []
    for f in Bar("file:" ).iter(glob.glob(ncfolder+"/**/*.nc",recursive=True)):
        ncdf = Dataset(f, 'r')  # Dataset is the class behavior to open the file
        #pdb.set_trace()
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
                            prof=Profile(eyed,profiledata,"insitu","practical")
                            profiles.append(prof)
                            #print(profiledata)
        del ncdf
    return profiles


def extractDeepArgoProfiles(jsonfile): 
    profiles= []
    with open(jsonfile) as f:
        data = json.load(f)
        for prof in Bar("deep profiles").iter(data):
            profiledata = {}
            profiledata["lon"] = prof["lon"]
            profiledata["lat"] = prof["lat"]
            if geoFilter(profiledata["lon"],profiledata["lat"]):
                profiledata["cruise"] = prof["platform_number"]
                profiledata["station"] = prof["cycle_number"]
                profiledata["time"] = prof["date"]
                profiledata["sal"]=[]
                profiledata["temp"]=[]
                profiledata["pres"]=[]
                for m in prof["measurements"]:
                    if "psal" in m.keys() and "temp" in m.keys():
                        profiledata["sal"].append(m["psal"])
                        profiledata["temp"].append(m["temp"])
                        profiledata["pres"].append(m["pres"])
                s = np.argsort(profiledata["pres"])
                profiledata["sal"] = np.asarray(profiledata["sal"])[s]
                profiledata["temp"] = np.asarray(profiledata["temp"])[s]
                profiledata["pres"] = np.asarray(profiledata["pres"])[s]
                #profiledata["cruise"] = ncdf.WOCE_ID
                #profiledata["station"] = ncdf.STATION_NUMBER

                if len(profiledata["pres"])>4 and max(profiledata["pres"])>1500:

                        if {"sal","temp","pres","lat","lon"}.issubset(profiledata.keys())\
                            and abs(max(profiledata["pres"])-min(profiledata["pres"])) > 100:
                            eyed=idgenerator()
                            prof=Profile(eyed,profiledata,"insitu","practical")
                            profiles.append(prof)
    return profiles


def createMesh(xvals,yvals,coord="xy",spacingscale=1):
    if coord == "xy":
        n=25
        x1,y1 = singleXY((-68,-54))
        x2,y2 = singleXY((-21,-7))
        xmin = min(x1,x2)
        xmax = max(x1,x2)
        ymin = min(y1,y2)
        ymax = max(y1,y2)
        return np.meshgrid(np.linspace(xmin,xmax,n), np.linspace(ymin,ymax,n),indexing="xy")
    if coord =="latlon":
        grd = np.meshgrid(\
                np.linspace(-50,-8,42*spacingscale),np.linspace(42,12,30*spacingscale))
        for i in range(grd[0].shape[0]):
            for j in range(grd[0].shape[1]):
                x,y = singleXY((grd[0][i][j],grd[1][i][j]))
                grd[0][i][j] = x
                grd[1][i][j] = y
        return grd
    else:
        print("This grid type is not supported", sys.exc_info()[0])


mapbounds = {"lllon":-53,"urlon":-10,"lllat":-45,"urlat":-8,\
        "lat_0":-17,"lon_0":-33}


def getWOCEExpoCodes(ncfolder): 
    expocodes= {}
    for f in glob.glob(ncfolder+"/*.nc"):
        ncdf = Dataset(f, 'r')  # Dataset is the class behavior to open the file
        expocodes.setdefault(ncdf.EXPOCODE,0)
        expocodes[ncdf.EXPOCODE] = expocodes[ncdf.EXPOCODE] + 1

    a_file = open("wocecruises.csv", "w")
    writer = csv.writer(a_file)
    for key, value in expocodes.items():
        writer.writerow([key, value])
    a_file.close()
 




