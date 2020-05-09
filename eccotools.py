from netCDF4 import Dataset
from profile import Profile
from progress.bar import Bar
import numpy as np
import pickle
import xarray as xr
from regionlib import brasil
from regionlib import nepb
import matplotlib.pyplot as plt
from geopy.distance import great_circle
import pdb
import graph
import gsw
import xmitgcm
from scipy.io import loadmat
from random import randint

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


#generate a unique id
def idgenerator():
    if not hasattr(idgenerator,"id"):
        idgenerator.id=0
    idgenerator.id = idgenerator.id+1
    return idgenerator.id

def arcticRestrict(lat,lon):
    return (lat >= 68 and not (lat<81 and -93 < lon < 20))

def nepbRestrict(lat,lon):
    return 60>lat> 20 and (0>lon > -181 or lon>170)

def brasilRestrict(lat,lon):
    latinrange = (lat<0 and lat >-80)
    loninrange = (lon>-69 and lon < -12)
    return (latinrange and loninrange)

def genPsi(uvel,vvel,dx,dy,xc,yc):
     uvel = np.asarray(uvel)
     vvel = np.asarray(vvel)
     uvel = np.nanmean(uvel,axis=0)
     vvel = np.nanmean(vvel,axis=0)
     psi = np.full_like(uvel,-99999)
     queue = [np.asarray((0,89))]
     for l in range(psi.shape[0]):
         psi[l][queue[0][0]][queue[0][1]] = 0
     tuplesteps = ((1,0),(0,1),(-1,0),(0,-1))
     steps=[]
     dx = np.asarray(dx)
     dy = np.asarray(dy)
     for s in tuplesteps:
         steps.append(np.asarray(s))

     count = 0
     figcount=0
     while len(queue)>0:
         curr = queue.pop(randint(0,len(queue)-1))
         # curr = queue.pop(0)
         for step in steps:
             n = step+curr
             if np.max(n) <90 and np.min(n)>=0 and psi[0][n[0]][n[1]] < -90000 and ~np.isnan(psi[0][n[0]][n[1]]):
                queue.append(step+curr)
                for d in range(uvel.shape[0]):
                    if ~np.isnan(psi[d][curr[0]][curr[1]]):
                        vavg = np.mean((vvel[d][n[0]][n[1]],vvel[d][curr[0]][curr[1]]))
                        uavg = np.mean((uvel[d][n[0]][n[1]],uvel[d][curr[0]][curr[1]]))
                        psid = np.array(uavg*dy[curr[0]][curr[1]],-vavg*dx[curr[0]][curr[1]])
                        newpsi = gsw.f(yc[curr[0]][curr[1]])*np.sum(psid.dot(step)) + psi[d][curr[0]][curr[1]]
                        if ~np.isnan(newpsi):
                            psi[d][n[0]][n[1]] = newpsi
                        else:
                            psi[d][n[0]][n[1]] = np.nan
         if count>1000:
             figcount+=1
             count=0
             g = np.full_like(psi[10],np.nan)
             g[psi[10] > -90000] = psi[10][psi[10] > -90000]
             print(np.nanmin(g))
             plt.imshow(g)
             plt.colorbar()
             plt.savefig("../arcticcirc-pics/surfaces/psisearch/"+str(figcount))
             plt.close()
         count+=1
     fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)
     ax1.imshow(uvel[10])
     ax2.imshow(np.multiply(np.diff(psi[10],axis=1),1.0/gsw.f(yc[:90,:89]))/dy[:90,:89])
     ax3.imshow(vvel[10])
     ax4.imshow(np.multiply(np.diff(psi[10],axis=0),1.0/gsw.f(yc[:89,:90]))/dx[:89,:90])
     plt.show()
     return psi


#read nc files, load into profiles and save into pickle
def generateProfilesNative(prefix,coordFilter,savepath='data/eccoprofiles.pickle'):
    prefixToMixCoord={"BRASILEAST":(90,0),"BRASILWEST":(90,270),"NEPBWEST":(180,180),"NEPBEAST":(180,270)}
    thetaset= xr.open_dataset('ecco/TILEDATA/'+prefix+'THETA.nc',decode_times=False)
    saltset = xr.open_dataset('ecco/TILEDATA/'+prefix+'SALT.nc',decode_times=False)
    uset = xr.open_dataset('ecco/TILEDATA/'+prefix+'EVEL.nc',decode_times=False)
    vset = xr.open_dataset('ecco/TILEDATA/'+prefix+'NVEL.nc',decode_times=False)
    depths = thetaset.dep
    ecco_grid = xr.open_dataset('ecco/TILEDATA/'+prefix+'GRID.nc')
    lons = ecco_grid.XC
    lats = ecco_grid.YC
    dy = ecco_grid.DYC
    dx = ecco_grid.DXC
    psi = genPsi(uset["EVEL"],vset["NVEL"],dx,dy,lons,lats)
    #diffkr,kapredi,kapgm = formatMixData("diffkr",prefix),formatMixData("kapredi",prefix),formatMixData("kapgm",prefix)
    profiles = []
    llc90_extra_metadata = xmitgcm.utils.get_extra_metadata(domain='llc', nx=90)
    grid = xmitgcm.utils.get_grid_from_input('./ecco/mixingdata/nctiles_grid/tile<NFACET>.mitgrid',geometry='llc',extra_metadata=llc90_extra_metadata)
    print("LETS RIP")
    diffkrField  = loadmat('./ecco/mixingdata/diffkr.mat')["finalField"]
    kaprediField  = loadmat('./ecco/mixingdata/kapredi.mat')["finalField"]
    kapgmField  = loadmat('./ecco/mixingdata/kapgm.mat')["finalField"]
    #plt.imshow(diffkrField[:,:,3])

    #plt.colorbar()
    #plt.show()
    latgraph,longraph,diffkrgraph = [], [],[]
    print(np.min(lats))
    for i in Bar("Row: ").iter(range(90)):
        for j in range(90):
            lon = lons[i][j]
            lat = lats[i][j]
            #latgraph.append(lat)
            #if lon<0:lon+=360
            #longraph.append(lon)
            #diffkrgraph.append(diffkrField[prefixToMixCoord[prefix][0]+89-j,prefixToMixCoord[prefix][1]+i,1])
    #plt.scatter(longraph,latgraph,c=diffkrgraph)
    #plt.show()
            if coordFilter(lat,lon):
                data = {}
                data["lat"]=lat
                data["relcoord"] = []
                data["lon"]=lon
                data["temp"]=[]
                data["sal"]=[]
                data["pres"]=[]
                data["knownu"]=[]
                data["knownv"]=[]
                data["kapredi"]=[]
                data["kapgm"]=[]
                data["diffkr"]=[]
                for depthindex in range(len(depths)):
                    if (~np.isnan(thetaset["THETA"].values[0][depthindex][i][j]) and ~np.isnan(saltset["SALT"].values[0][depthindex][i][j]) \
                            and ~np.isnan(thetaset["land"].values[0][i][j])):

                        data["pres"].append(float(depths.values[depthindex]))
                        t,s,u,v = [],[],[],[]
                        for month in range(12):
                            t.append(float(thetaset["THETA"].values[month][depthindex][i][j]))
                            s.append(float(saltset["SALT"].values[month][depthindex][i][j]))
                            u.append(float(uset["EVEL"].values[month][depthindex][i][j]))
                            v.append(float(vset["NVEL"].values[month][depthindex][i][j]))
                        data["temp"].append(np.mean(t))
                        data["relcoord"] = (x[i][j],y[i][j])
                        data["sal"].append(np.mean(s))
                        data["knownu"].append(np.mean(u))
                        data["knownv"].append(np.mean(v))

                        #data["diffkr"].append(diffkr[depthindex][i][j])
                        if prefix=="BRASILEAST":
                            data["diffkr"].append(diffkrField[prefixToMixCoord[prefix][0]+i,prefixToMixCoord[prefix][1]+j,depthindex])
                            data["kapgm"].append(kapgmField[prefixToMixCoord[prefix][0]+i,prefixToMixCoord[prefix][1]+j,depthindex])
                            data["kapredi"].append(kaprediField[prefixToMixCoord[prefix][0]+i,prefixToMixCoord[prefix][1]+j,depthindex])
                        if prefix=="BRASILWEST":
                            data["diffkr"].append(diffkrField[prefixToMixCoord[prefix][0]+89-j,prefixToMixCoord[prefix][1]+i,depthindex])
                            data["kapgm"].append(kapgmField[prefixToMixCoord[prefix][0]+89-j,prefixToMixCoord[prefix][1]+i,depthindex])
                            data["kapredi"].append(kaprediField[prefixToMixCoord[prefix][0]+89-j,prefixToMixCoord[prefix][1]+i,depthindex])
                        if prefix in ["NEPBEAST","NEPBWEST"]:
                            data["diffkr"].append(diffkrField[prefixToMixCoord[prefix][0]+89-j,prefixToMixCoord[prefix][1]+i,depthindex])
                            data["kapgm"].append(kapgmField[prefixToMixCoord[prefix][0]+89-j,prefixToMixCoord[prefix][1]+i,depthindex])
                            data["kapredi"].append(kaprediField[prefixToMixCoord[prefix][0]+89-j,prefixToMixCoord[prefix][1]+i,depthindex])

                if len(data["pres"])>4 and max(data["pres"])>1500:
                    eyed=idgenerator()
                    p=Profile(eyed,data,tempunit="potential",salunit="practical")
                    profiles.append(p)
                

    return profiles


#p2 = generateProfilesNative("BRASILEAST",brasilRestrict,"data/eccobrasilprofiles.pickle")
#p1 = generateProfilesNative("BRASILWEST",brasilRestrict,"data/eccobrasilprofiles.pickle")
#with open("data/eccobrasilprofiles.pickle", 'wb') as outfile:
    #pickle.dump(p1+p2, outfile)
#generateProfilesNative("ARCTIC",arcticRestrict,"data/ecconprofiles.pickle")
p1  = generateProfilesNative("NEPBWEST",nepbRestrict,"data/ecconepbprofiles.pickle")
p2  = generateProfilesNative("NEPBEAST",nepbRestrict,"data/ecconepbprofiles.pickle")
with open("data/ecconepbprofiles.pickle", 'wb') as outfile:
    pickle.dump(p1+p2, outfile)
