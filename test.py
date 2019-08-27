from netCDF4 import Dataset
from profile import Profile
from progress.bar import Bar
import numpy as np
import pickle
import xarray as xr


def idgenerator():
    if not hasattr(idgenerator,"id"):
        idgenerator.id=0
    idgenerator.id = idgenerator.id+1
    return idgenerator.id

print(Profile)
thetaset = xr.open_dataset('ecco/THETA.0001.nc',decode_times=False)
lat = thetaset["lat"]
lon = thetaset["lon"]
depths = thetaset["dep"]
theta= xr.open_dataset('ecco/THETA.0001.nc',decode_times=False)["THETA"]

salt = xr.open_dataset('ecco/SALT.0001.nc',decode_times=False)["SALT"]
u = xr.open_dataset('ecco/EVEL.0001.nc',decode_times=False)["EVEL"]
v = xr.open_dataset('ecco/NVEL.0001.nc',decode_times=False)["NVEL"]

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
with open('data/eccoprofiles.pickle', 'wb') as outfile:
    pickle.dump(profiles, outfile)



