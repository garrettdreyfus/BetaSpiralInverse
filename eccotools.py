from netCDF4 import Dataset
from profile import Profile
from progress.bar import Bar
import numpy as np
import pickle
import xarray as xr
from regionlib import brasil
from geopy.distance import great_circle
import pdb
import graph
import pdb

##return index of certain coord
def closestGridPoint(x,y,prefix):
    if not hasattr(closestGridPoint,"grid"):
        closestGridPoint.grid = xr.open_dataset('ecco/TILEDATA/'+prefix+'GRID.nc',decode_times=False)
        lons = closestGridPoint.grid.XC
        print(lons)
        lats = closestGridPoint.grid.YC
        print("lons shape",lons.shape)
        theta = np.deg2rad(lons)
        r = ((90-lats)*111.0*1000.0)
        closestGridPoint.x = (r*np.cos(theta))
        closestGridPoint.y = (r*np.sin(theta))
        print("xs",closestGridPoint.x.shape)
    dists = (closestGridPoint.x-x)**2 + (closestGridPoint.y-y)**2
    if ~np.isnan(dists).all():
        loc = np.unravel_index(np.nanargmin(dists, axis=None), dists.shape)
        return loc
    else:
        return np.nan
    
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

def mixDataArray(arr,depth,tilenumber,lons,lats):
    arr = arr[105300*depth:105300*(depth+1)]
    arr = arr.byteswap().reshape([105300,1])
    arr = xr.DataArray(arr[(tilenumber-1)*8100:(tilenumber)*8100,0].reshape([90,90]),coords=[np.ravel(lons),np.ravel(lats)],
                            dims=['lon','lat'])
    #ds = xr.Dataset(data_vars={"m":(["i","j"],arr[(tilenumber-1)*8100:(tilenumber)*8100,0].reshape([90,90])), 
           #"lat":(["i","j"],lats),
           #"lon":(["i","j"],lons)}, 
        #coords={"x": (["i"],range(90)), 
            #"y": (["j"], range(90))})

    #ds = ds.set_coords(["lat","lon"])

    return arr

def formatMixData(var,prefix):
    prefixToTile = {"ARCTIC":7,"NEPB":8,"BRASILEAST":2,"BRASILWEST":12}
    dat = []
    tilenumber = prefixToTile[prefix]
    for depth in Bar("depth").iter(range(50)):
        grd = xr.open_dataset('ecco/TILEDATA/'+prefix+'GRID.nc',decode_times=False)
        lons = grd.XC
        lats = grd.YC
        geoflx = np.fromfile('ecco/mixingdata/'+var+'.bin', dtype=np.float32)
        geoflx06 = mixDataArray(geoflx,depth,tilenumber,lons,lats)
        geoflx = np.fromfile('ecco/mixingdata/'+var+'.data', dtype=np.float32)
        geoflx06 = mixDataArray(geoflx,depth,tilenumber,lons,lats)+geoflx06
        dat.append(geoflx06)
    return dat



def getMixAt(lon,lat,d,prefix):
    if not hasattr(getMixAt,"nvelset"):
        getMixAt.nvelset = xr.open_dataset('ecco/TILEDATA/'+prefix+'NVEL.nc',decode_times=False)
        getMixAt.depths = getMixAt.nvelset["dep"].values
        getMixAt.diffkr = formatMixData("diffkr")
        getMixAt.kapredi = formatMixData("kapredi")
        getMixAt.kapgm = formatMixData("kapgm")
    if d >= getMixAt.depths[0] and d <= getMixAt.depths[-1]:
        before = np.argwhere(getMixAt.depths<=d)[-1][0]
        after = np.argwhere(getMixAt.depths>=d)[0][0]
        for k in range(len(getMixAt.diffkr)):
            lonlow = np.nanmin(getMixAt.diffkr[k][before]["lon"])
            lonhigh = np.nanmax(getMixAt.diffkr[k][before]["lon"])
            latlow = np.nanmin(getMixAt.diffkr[k][before]["lat"])
            lathigh = np.nanmax(getMixAt.diffkr[k][before]["lat"])
            print(k,lonlow,lon,latlow)
            print(k,lonhigh,lat,lathigh)
            if lonlow< lon < lonhigh and latlow < lat < lathigh: 
                depthset = [getMixAt.depths[before],getMixAt.depths[after]]
                closest = (getMixAt.diffkr[k][before].lon-lon)**2 + (getMixAt.diffkr[k][before].lat-lat)**2
                print(np.argmin(closest))
                print("-"*5)
                print(getMixAt.diffkr[k][before][np.argmin(closest)])
                diffkrset= [getMixAt.diffkr[k][before][np.argmin(closest)]\
                        ,getMixAt.diffkr[k][after].sel(lon=lon,lat=lat,method="nearest")]
                kaprediset= [getMixAt.kapredi[k][before].sel(lon=lon,lat=lat,method="nearest")\
                        ,getMixAt.kapredi[k][after].sel(lon=lon,lat=lat,method="nearest")]
                kapgmset= [getMixAt.kapgm[k][before].sel(lon=lon,lat=lat,method="nearest")\
                        ,getMixAt.kapgm[k][after].sel(lon=lon,lat=lat,method="nearest")]
                diffkr = np.interp(d,depthset,diffkrset)
                kapredi = np.interp(d,depthset,kaprediset)
                kapgm = np.interp(d,depthset,kapgmset)
                #print(diffkr,kapredi,kapgm)
                return diffkr,kapredi,kapgm
        return np.nan,np.nan,np.nan
    else:
        return np.nan,np.nan,np.nan


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
            lon = surfaces[k]["lons"][l]
            lat = surfaces[k]["lats"][l]
            d = surfaces[k]["data"]["pres"][l]
            #print(getMixAt(lons,lats,d,prefix))
            diffkr,kapredi,kapgm = getMixAt(lon,lat,d,prefix)
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

def brasilRestrict(lat,lon):
    latinrange = (lat<0 and lat >-80)
    loninrange = (lon>-69 and lon < -12)
    return (latinrange and loninrange)



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
    diffkr,kapredi,kapgm = formatMixData("diffkr"),formatMixData("kapredi"),formatMixData("kapgm")
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
                data["knownu"]=[]
                data["knownv"]=[]
                for depthindex in range(len(depths)):
                    if ~np.isnan(thetaset["THETA"].values[0][depthindex][i][j]) and ~np.isnan(saltset["SALT"].values[0][depthindex][i][j]) \
                            and ~np.isnan(thetaset["land"].values[0][i][j]):

                        data["pres"].append(float(depths.values[depthindex]))
                        t,s,u,v = [],[],[],[]
                        for month in range(12):
                            t.append(float(thetaset["THETA"].values[month][depthindex][i][j]))
                            s.append(float(saltset["SALT"].values[month][depthindex][i][j]))
                            u.append(float(uset["EVEL"].values[month][depthindex][i][j]))
                            v.append(float(vset["NVEL"].values[month][depthindex][i][j]))
                        data["temp"].append(np.mean(t))
                        data["sal"].append(np.mean(s))
                        data["knownu"].append(np.mean(u))
                        data["knownv"].append(np.mean(v))

                        data["kapredi"].append(kapredi[depthindex][i][j])
                        data["diffkr"].append(diffkr[depthindex][i][j])
                        data["kapgm"].append(kapgm[depthindex][i][j])

                if len(data["pres"])>4 and max(data["pres"])>1500:
                    eyed=idgenerator()
                    p=Profile(eyed,data,tempunit="potential",salunit="practical")
                    profiles.append(p)
                

    graph.plotProfiles(brasil,profiles,"profs")
    return profiles


#p1 = generateProfilesNative("BRASILWEST",brasilRestrict,"data/eccobrasilprofiles.pickle")
#p2 = generateProfilesNative("BRASILEAST",brasilRestrict,"data/eccobrasilprofiles.pickle")
#with open("data/eccobrasilprofiles.pickle", 'wb') as outfile:
    #pickle.dump(p1+p2, outfile)
#generateProfilesNative("ARCTIC",arcticRestrict,"data/ecconprofiles.pickle")
#generateProfilesNative("NEPB",nepbRestrict,"data/ecconepbprofiles.pickle")
