from netCDF4 import Dataset
from profile import Profile
from progress.bar import Bar
import numpy as np
import pickle
import xarray as xr
from geopy.distance import great_circle

##return index of certain coord
def getLatLonIndex(latin,lonin):
    thetaset = xr.open_dataset('ecco/THETA.2015.nc',decode_times=False)
    lat = thetaset["lat"]
    lon = thetaset["lon"]
    latindex = int(((round(latin*10)-3)-((round(latin*10)-3)%5))/5) + 180
    lonindex = int(((round(lonin*10)-3)-((round(lonin*10)-3)%5))/5) + 360
    dist = great_circle((latin,lonin),(lat[latindex][lonindex],lon[latindex][lonindex])).km
    return latindex,lonindex,dist

#return the ssh at a coord
def getSSHAt(latin,lonin):
    latindex,lonindex,dist = getLatLonIndex(latin,lonin)
    if not hasattr(getSSHAt,"sshset"):
        getSSHAt.sshset = xr.open_dataset('ecco/SSH.2015.nc',decode_times=False)
    return getSSHAt.sshset["SSH"][0][latindex][lonindex]

#return the vel at a coord
def getVelAt(latin,lonin,d):
    latindex,lonindex,dist = getLatLonIndex(latin,lonin)
    if not hasattr(getVelAt,"nvelset"):
        getVelAt.nvelset = xr.open_dataset('ecco/NVEL.2015.nc',decode_times=False)
        getVelAt.evelset = xr.open_dataset('ecco/EVEL.2015.nc',decode_times=False)
        getVelAt.depths = getVelAt.nvelset["dep"].values
    if d >= getVelAt.depths[0] and d <= getVelAt.depths[-1]:
        before = np.argwhere(getVelAt.depths<=d)[-1][0]
        after = np.argwhere(getVelAt.depths>=d)[0][0]
        depthset = [getVelAt.depths[before],getVelAt.depths[after]]
        nveldepthset= [getVelAt.nvelset["NVEL"][0][before][latindex][lonindex],getVelAt.nvelset["NVEL"][0][after][latindex][lonindex]]
        eveldepthset= [getVelAt.evelset["EVEL"][0][before][latindex][lonindex],getVelAt.evelset["EVEL"][0][after][latindex][lonindex]]
        nvel = np.interp(d,depthset,nveldepthset)
        evel = np.interp(d,depthset,eveldepthset)
        return evel,nvel,dist
    else:
        return 0,0,0


##add ssh values to a surface
def addSSHToSurface(surfaces):
    dumbcache = {}
    for k in surfaces.keys():
        surfaces[k]["data"]["ssh"] =  np.full(len(surfaces[k]["lons"]),np.nan)
        for l in range(len(surfaces[k]["lats"])):
            lat = surfaces[k]["lats"][l]
            lon = surfaces[k]["lons"][l]
            if (lat,lon) not in dumbcache.keys():
                dumbcache[(lat,lon)]=getSSHAt(lat,lon)
            surfaces[k]["data"]["ssh"][l] = dumbcache[(lat,lon)]
    return surfaces

def addModelEccoUV(surfaces):
    distances = []
    for k in Bar("Adding model uv").iter(surfaces.keys()):
        surfaces[k]["data"]["knownu"] =  np.full(len(surfaces[k]["lons"]),np.nan)
        surfaces[k]["data"]["knownv"] =  np.full(len(surfaces[k]["lons"]),np.nan)
        for l in range(len(surfaces[k]["lats"])):
            lat = surfaces[k]["lats"][l]
            lon = surfaces[k]["lons"][l]
            d = surfaces[k]["data"]["pres"][l]
            u,v,dist = getVelAt(lat,lon,d)
            distances.append(dist)
            surfaces[k]["data"]["knownu"][l] = u
            surfaces[k]["data"]["knownv"][l] = v
    print("distances mean: ",np.mean(distances))
    print("distances max: ",np.max(distances))
    print("distances stdev: ",np.std(distances))
    return surfaces

 

#generate a unique id
def idgenerator():
    if not hasattr(idgenerator,"id"):
        idgenerator.id=0
    idgenerator.id = idgenerator.id+1
    return idgenerator.id

#read nc files, load into profiles and save into pickle
def generateProfiles(savepath='data/eccoprofiles.pickle'):
    thetaset = xr.open_dataset('ecco/THETA.2015.nc',decode_times=False)
    lat = thetaset["lat"]
    lon = thetaset["lon"]
    depths = thetaset["dep"]
    theta= xr.open_dataset('ecco/THETA.2015.nc',decode_times=False)["THETA"]

    salt = xr.open_dataset('ecco/SALT.2015.nc',decode_times=False)["SALT"]
    u = xr.open_dataset('ecco/EVEL.2015.nc',decode_times=False)["EVEL"]
    v = xr.open_dataset('ecco/NVEL.2015.nc',decode_times=False)["NVEL"]

    profiles = []
    print("LETS RIP")
    for latindex in Bar("lats").iter(range(len(lat))):
        if lat[latindex][0] >= 68:
            for lonindex in Bar("lons").iter(range(len(lon[latindex]))):
                if lonindex%5==0 and  not (lat[latindex][0]<81 and -93 < lon[latindex][lonindex] < 20):
                    data = {}
                    data["lat"]=lat.values[latindex][lonindex]
                    data["lon"]=lon.values[latindex][lonindex]
                    data["temp"]=[]
                    data["sal"]=[]
                    data["pres"]=[]
                    data["knownu"]=[]
                    data["knownv"]=[]
                    for depthindex in range(len(depths)):
                        if ~np.isnan(theta.values[0][depthindex][latindex][lonindex]) and ~np.isnan(salt.values[0][depthindex][latindex][lonindex]):
                            data["pres"].append(float(depths.values[depthindex]))
                            data["temp"].append(float(theta.values[0][depthindex][latindex][lonindex]))
                            data["sal"].append(float(salt.values[0][depthindex][latindex][lonindex]))
                            data["knownu"].append(float(u.values[0][depthindex][latindex][lonindex]))
                            data["knownv"].append(float(v.values[0][depthindex][latindex][lonindex]))

                    if len(data["pres"])>4 and max(data["pres"])>1500:
                        eyed=idgenerator()
                        p=Profile(eyed,data)
                        profiles.append(p)
                

    print(len(profiles))
    with open(savepath, 'wb') as outfile:
        pickle.dump(profiles, outfile)

getLatLonIndex(87.9,-56.9)


