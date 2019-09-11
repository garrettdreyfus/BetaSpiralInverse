from netCDF4 import Dataset
from profile import Profile
from progress.bar import Bar
import numpy as np
import pickle
import xarray as xr
from geopy.distance import great_circle
import pdb
import graph

##return index of certain coord
def closestGridPoint(x,y,prefix):
    if not hasattr(closestGridPoint,"grid"):
        closestGridPoint.grid = xr.open_dataset('ecco/TILEDATA/'+prefix+'GRID.nc',decode_times=False)
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
def getVelAt(x,y,d,prefix):
    loc = closestGridPoint(x,y,prefix)
    if not hasattr(getVelAt,"nvelset"):
        getVelAt.nvelset = xr.open_dataset('ecco/TILEDATA/'+prefix+'NVEL.nc',decode_times=False)
        getVelAt.evelset = xr.open_dataset('ecco/TILEDATA/'+prefix+'EVEL.nc',decode_times=False)
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

def extractArctic(arr,depth,tilenumber):
    arr = arr[105300*depth:105300*(depth+1)]
    arr = arr.byteswap().reshape([105300,1])

    arr = xr.DataArray(arr[(tilenumber-1)*8100:(tilenumber)*8100,0].reshape([90,90]),coords=[np.arange(0,90,1),np.arange(0,90,1)],
                            dims=['j','i'])
    return arr

def formatMixData(var,prefix):
    prefixToTile = {"ARCTIC":7,"NEPB":8}
    tilenumber = prefixToTile[prefix]
    dat = []
    for depth in Bar("depth").iter(range(50)):
        geoflx = np.fromfile('ecco/mixingdata/'+var+'.bin', dtype=np.float32)
        geoflx06 = extractArctic(geoflx,depth,tilenumber)
        geoflx = np.fromfile('ecco/mixingdata/'+var+'.data', dtype=np.float32)
        geoflx06 = extractArctic(geoflx,depth,tilenumber)+geoflx06
        dat.append(geoflx06)
    return dat



def getMixAt(x,y,d,prefix):
    loc = closestGridPoint(x,y,prefix)
    if not hasattr(getMixAt,"nvelset"):
        getMixAt.nvelset = xr.open_dataset('ecco/TILEDATA/'+prefix+'NVEL.nc',decode_times=False)
        getMixAt.depths = getMixAt.nvelset["dep"].values
        getMixAt.diffkr = formatMixData("diffkr",prefix)
        getMixAt.kapredi = formatMixData("kapredi",prefix)
        getMixAt.kapgm = formatMixData("kapgm",prefix)
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


def addModelEccoUV(surfaces,prefix):
    for k in Bar("Adding model uv").iter(surfaces.keys()):
        surfaces[k]["data"]["knownu"] =  np.full(len(surfaces[k]["lons"]),np.nan)
        surfaces[k]["data"]["knownv"] =  np.full(len(surfaces[k]["lons"]),np.nan)
        for l in range(len(surfaces[k]["lats"])):
            x = surfaces[k]["x"][l]
            y = surfaces[k]["y"][l]
            d = surfaces[k]["data"]["pres"][l]
            u,v = getVelAt(x,y,d,prefix)
            surfaces[k]["data"]["knownu"][l] = u
            surfaces[k]["data"]["knownv"][l] = v
    return surfaces

def addModelEccoMix(surfaces,prefix):
    for k in Bar("Adding model mix").iter(surfaces.keys()):
        surfaces[k]["data"]["diffkr"] =  np.full(len(surfaces[k]["lons"]),np.nan)
        surfaces[k]["data"]["kapredi"] =  np.full(len(surfaces[k]["lons"]),np.nan)
        surfaces[k]["data"]["kapgm"] =  np.full(len(surfaces[k]["lons"]),np.nan)
        for l in range(len(surfaces[k]["lats"])):
            x = surfaces[k]["x"][l]
            y = surfaces[k]["y"][l]
            d = surfaces[k]["data"]["pres"][l]
            diffkr,kapredi,kapgm = getMixAt(x,y,d,prefix)
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

def arcticRestrict(lat,lon):
    return (lat >= 68 and not (lat<81 and -93 < lon < 20))

def nepbRestrict(lat,lon):
    return 60>lat> 20 and 0>lon > -170

#read nc files, load into profiles and save into pickle
def generateProfilesNative(prefix,coordFilter,savepath='data/eccoprofiles.pickle'):
    thetaset= xr.open_dataset('ecco/TILEDATA/'+prefix+'THETA.nc',decode_times=False)
    saltset = xr.open_dataset('ecco/TILEDATA/'+prefix+'SALT.nc',decode_times=False)
    uset = xr.open_dataset('ecco/TILEDATA/'+prefix+'EVEL.nc',decode_times=False)
    vset = xr.open_dataset('ecco/TILEDATA/'+prefix+'NVEL.nc',decode_times=False)
    depths = thetaset.dep
    ecco_grid = xr.open_dataset('ecco/TILEDATA/'+prefix+'GRID.nc')
    lons = ecco_grid.XC
    lats = ecco_grid.YC

    profiles = []
    print("LETS RIP")
    for i in Bar("Row: ").iter(range(90)):
        for j in range(90):
            lon = lons[i][j]
            lat = lats[i][j]
            if coordFilter(lat,lon):
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
    graph.plotProfiles(profiles,"profs",region="nepb")
    with open(savepath, 'wb') as outfile:
        pickle.dump(profiles, outfile)


#generateProfilesNative("ARCTIC",arcticRestrict,"data/ecconprofiles.pickle")
#generateProfilesNative("NEPB",nepbRestrict,"data/ecconepbprofiles.pickle")
