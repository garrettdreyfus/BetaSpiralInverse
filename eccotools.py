from netCDF4 import Dataset
from profile import Profile
from progress.bar import Bar
import numpy as np
import pickle
import xarray as xr
from geopy.distance import great_circle
import pdb

##return index of certain coord
def closestGridPoint(x,y):
    if not hasattr(closestGridPoint,"grid"):
        closestGridPoint.grid = xr.open_dataset('ecco/TILEDATA/ARCTICGRID.nc',decode_times=False)
        lons = closestGridPoint.grid.XC
        lats = closestGridPoint.grid.YC
        print("lons shape",lons.shape)
        theta = np.deg2rad(lons)
        r = ((90-lats)*111.0*1000.0)
        closestGridPoint.x = (r*np.cos(theta))
        closestGridPoint.y = (r*np.sin(theta))
        print("xs",closestGridPoint.x.shape)
    dists = (closestGridPoint.x-x)**2 + (closestGridPoint.y-y)**2
    loc = np.unravel_index(np.argmin(dists, axis=None), dists.shape)
    return loc
    
    #return the ssh at a coord
def getSSHAt(latin,lonin):
    latindex,lonindex,dist = getLatLonIndex(latin,lonin)
    if not hasattr(getSSHAt,"sshset"):
        getSSHAt.sshset = xr.open_dataset('ecco/SSH.2015.nc',decode_times=False)
    return getSSHAt.sshset["SSH"][0][latindex][lonindex]

#return the vel at a coord
def getVelAt(x,y,d):
    loc = closestGridPoint(x,y)
    if not hasattr(getVelAt,"nvelset"):
        getVelAt.nvelset = xr.open_dataset('ecco/TILEDATA/NVELARCTIC.nc',decode_times=False)
        getVelAt.evelset = xr.open_dataset('ecco/TILEDATA/EVELARCTIC.nc',decode_times=False)
        getVelAt.depths = getVelAt.nvelset["dep"].values
    if d >= getVelAt.depths[0] and d <= getVelAt.depths[-1]:
        before = np.argwhere(getVelAt.depths<=d)[-1][0]
        after = np.argwhere(getVelAt.depths>=d)[0][0]
        depthset = [getVelAt.depths[before],getVelAt.depths[after]]
        nveldepthset= [getVelAt.nvelset["NVEL"][0][before][loc],getVelAt.nvelset["NVEL"][0][after][loc]]
        eveldepthset= [getVelAt.evelset["EVEL"][0][before][loc],getVelAt.evelset["EVEL"][0][after][loc]]
        nvel = np.interp(d,depthset,nveldepthset)
        evel = np.interp(d,depthset,eveldepthset)
        return evel,nvel
    else:
        return 0,0,0

def extractArctic(arr,depth):
    arr = arr[105300*depth:105300*(depth+1)]
    arr = arr.byteswap().reshape([105300,1])

    arr = xr.DataArray(arr[48600:56700,0].reshape([90,90]),coords=[np.arange(0,90,1),np.arange(0,90,1)],
                            dims=['j','i'])
    return arr

def formatMixData(var):
    dat = []
    for depth in Bar("depth").iter(range(50)):
        geoflx = np.fromfile('ecco/mixingdata/'+var+'.bin', dtype=np.float32)
        geoflx06 = extractArctic(geoflx,depth)
        geoflx = np.fromfile('ecco/mixingdata/'+var+'.data', dtype=np.float32)
        geoflx06 = extractArctic(geoflx,depth)+geoflx06
        dat.append(geoflx06)
    return dat



def getMixAt(x,y,d):
    loc = closestGridPoint(x,y)
    if not hasattr(getMixAt,"nvelset"):
        getMixAt.nvelset = xr.open_dataset('ecco/TILEDATA/NVELARCTIC.nc',decode_times=False)
        getMixAt.depths = getMixAt.nvelset["dep"].values
        getMixAt.diffkr = formatMixData("diffkr")
        getMixAt.kapredi = formatMixData("kapredi")
        getMixAt.kapgm = formatMixData("kapgm")
    if d >= getMixAt.depths[0] and d <= getMixAt.depths[-1]:
        before = np.argwhere(getMixAt.depths<=d)[-1][0]
        after = np.argwhere(getMixAt.depths>=d)[0][0]
        depthset = [getMixAt.depths[before],getMixAt.depths[after]]
        diffkrset= [getMixAt.diffkr[before][loc],getMixAt.diffkr[after][loc]]
        kaprediset= [getMixAt.kapredi[before][loc],getMixAt.kapredi[after][loc]]
        kapgmset= [getMixAt.kapgm[before][loc],getMixAt.kapgm[after][loc]]
        diffkr = np.interp(d,depthset,diffkrset)
        kapredi = np.interp(d,depthset,kaprediset)
        kapgm = np.interp(d,depthset,kapgmset)
        return diffkr,kapredi,kapgm
    else:
        return 0,0,0

##add ssh values to a surface
#def addSSHToSurface(surfaces):
    #dumbcache = {}
    #for k in surfaces.keys():
        #surfaces[k]["data"]["ssh"] =  np.full(len(surfaces[k]["lons"]),np.nan)
        #for l in range(len(surfaces[k]["lats"])):
            #lat = surfaces[k]["lats"][l]
            #lon = surfaces[k]["lons"][l]
            #if (lat,lon) not in dumbcache.keys():
                #dumbcache[(lat,lon)]=getSSHAt(lat,lon)
            #surfaces[k]["data"]["ssh"][l] = dumbcache[(lat,lon)]
            #surfaces[k]["data"]["psi"][l] = surfaces[k]["data"]["psi"][l]-9.8*dumbcache[(lat,lon)]
    #return surfaces


def addModelEccoUV(surfaces):
    for k in Bar("Adding model uv").iter(surfaces.keys()):
        surfaces[k]["data"]["knownu"] =  np.full(len(surfaces[k]["lons"]),np.nan)
        surfaces[k]["data"]["knownv"] =  np.full(len(surfaces[k]["lons"]),np.nan)
        for l in range(len(surfaces[k]["lats"])):
            x = surfaces[k]["x"][l]
            y = surfaces[k]["y"][l]
            d = surfaces[k]["data"]["pres"][l]
            u,v = getVelAt(x,y,d)
            surfaces[k]["data"]["knownu"][l] = u
            surfaces[k]["data"]["knownv"][l] = v
    return surfaces

def addModelEccoMix(surfaces):
    for k in Bar("Adding model mix").iter(surfaces.keys()):
        surfaces[k]["data"]["diffkr"] =  np.full(len(surfaces[k]["lons"]),np.nan)
        surfaces[k]["data"]["kapredi"] =  np.full(len(surfaces[k]["lons"]),np.nan)
        surfaces[k]["data"]["kapgm"] =  np.full(len(surfaces[k]["lons"]),np.nan)
        for l in range(len(surfaces[k]["lats"])):
            x = surfaces[k]["x"][l]
            y = surfaces[k]["y"][l]
            d = surfaces[k]["data"]["pres"][l]
            diffkr,kapredi,kapgm = getMixAt(x,y,d)
            surfaces[k]["data"]["diffkr"][l] = diffkr
            surfaces[k]["data"]["kapgm"][l] = kapgm
            surfaces[k]["data"]["kapredi"][l] = kapredi
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
                if not (lat[latindex][0]<81 and -93 < lon[latindex][lonindex] < 20):
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

#read nc files, load into profiles and save into pickle
def generateProfilesNative(savepath='data/eccoprofiles.pickle'):
    thetaset= xr.open_dataset('ecco/TILEDATA/THETAARCTIC.nc',decode_times=False)
    saltset = xr.open_dataset('ecco/TILEDATA/SALTARCTIC.nc',decode_times=False)
    uset = xr.open_dataset('ecco/TILEDATA/EVELARCTIC.nc',decode_times=False)
    vset = xr.open_dataset('ecco/TILEDATA/NVELARCTIC.nc',decode_times=False)
    depths = thetaset.dep
    ecco_grid = xr.open_dataset('ecco/TILEDATA/ARCTICGRID.nc')
    lons = ecco_grid.XC
    lats = ecco_grid.YC

    profiles = []
    print("LETS RIP")
    for i in Bar("Row: ").iter(range(90)):
        for j in range(90):
            lon = lons[i][j]
            lat = lats[i][j]
            if lat >= 68 and not (lat<81 and -93 < lon < 20):
                data = {}
                data["lat"]=lat
                data["lon"]=lon
                data["temp"]=[]
                data["sal"]=[]
                data["pres"]=[]
                for depthindex in range(len(depths)):
                    if ~np.isnan(thetaset["THETA"].values[0][depthindex][i][j]) and ~np.isnan(saltset["SALT"].values[0][depthindex][i][j]):
                        data["pres"].append(float(depths.values[depthindex]))
                        data["temp"].append(float(thetaset["THETA"].values[0][depthindex][i][j]))
                        data["sal"].append(float(saltset["SALT"].values[0][depthindex][i][j]))

                if len(data["pres"])>4 and max(data["pres"])>1500:
                    eyed=idgenerator()
                    p=Profile(eyed,data)
                    profiles.append(p)
                

    print(len(profiles))
    graph.plotProfiles(profiles,"profs")
    with open(savepath, 'wb') as outfile:
        pickle.dump(profiles, outfile)

