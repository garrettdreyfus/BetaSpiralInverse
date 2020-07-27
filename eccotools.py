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
import itertools
import inverttools
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

def genPsi(uvel,vvel,dx,dy,xc,yc,d):
     uvel = np.asarray(uvel)
     vvel = np.asarray(vvel)
     uvel = np.nanmean(uvel,axis=0)
     vvel = np.nanmean(vvel,axis=0)
     for j in range(uvel.shape[0]):
         uvel[j] = uvel[j].T
         vvel[j] = vvel[j].T

     dx = dx.T
     dy = dy.T
     xc = xc.T
     yc = yc.T

     psi = np.full_like(uvel,-99999)
     queue = [(89,89)]
     for l in range(psi.shape[0]):
         psi[l][queue[0][0]][queue[0][1]] = 0
     tuplesteps = ((1,0),(-1,0),(0,1),(0,-1))
     steps=[]
     dx = np.asarray(dx)
     dy = np.asarray(dy)
     for s in tuplesteps:
         steps.append(np.asarray(s))
     count = 0
     figcount=0
     while len(queue)>0:
         curr = queue.pop(0)
         for step in steps:
             n = step+curr
             if np.max(n) <90 and np.min(n)>=0 and psi[d][n[0]][n[1]] < -90000:
                if ~np.isnan(uvel[d][curr[0]][curr[1]]):
                    vavg = np.mean((vvel[d][n[0]][n[1]],vvel[d][curr[0]][curr[1]]))
                    uavg = np.mean((uvel[d][n[0]][n[1]],uvel[d][curr[0]][curr[1]]))
                    psid = np.array(vavg*dx[curr[0]][curr[1]],-uavg*dy[curr[0]][curr[1]])
                    newpsi = gsw.f(yc[curr[0]][curr[1]])*np.sum(psid.dot(step)) + psi[d][curr[0]][curr[1]]
                    if ~np.isnan(newpsi):
                        psi[d][n[0]][n[1]] = newpsi
                        queue.append(step+curr)
         if count>100:
             figcount+=1
             count=0
             g = np.full_like(psi[d],np.nan)
             g[psi[d] > -90000] = psi[d][psi[d] > -90000]
             plt.imshow(g,interpolation=None)
             plt.colorbar()
             plt.savefig("../arcticcirc-pics/surfaces/psisearch/"+str(figcount))
             plt.close()
         count+=1
     fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)
     g = np.full_like(psi[d],np.nan)
     g[psi[d] > -90000] = psi[d][psi[d] > -90000]
     boxIntegral(uvel[d],vvel[d],g,(80,80),(85,85),dx,dy)
     ax1.imshow(uvel[d],interpolation=None)
     ax2.imshow(np.multiply(np.diff(g,axis=1),1.0/gsw.f(yc[:90,:89]))/dy[:90,:89],interpolation=None)
     ax3.imshow(vvel[d],interpolation=None)
     ax4.imshow(np.multiply(np.diff(g,axis=0),1.0/gsw.f(yc[:89,:90]))/dx[:89,:90],interpolation=None)
     plt.show()
     return psi

def generateDivMatrix(uvel,vvel,dx,dy,yc):
    dudx = np.diff(uvel,axis=1)[:uvel.shape[0]-1,:uvel.shape[1]-1]/dx[:uvel.shape[0]-1,:uvel.shape[1]-1]
    dvdy = np.diff(vvel,axis=0)[:uvel.shape[0]-1,:uvel.shape[1]-1]/dy[:uvel.shape[0]-1,:uvel.shape[1]-1]
    return np.multiply(dudx + dvdy,gsw.f(yc[:uvel.shape[0]-1,:uvel.shape[1]-1]))

def generateCurlMatrix(uvel,vvel,dx,dy,yc):
    dudy = np.diff(uvel,axis=0)[:uvel.shape[0]-1,:uvel.shape[1]-1]/dy[:uvel.shape[0]-1,:uvel.shape[1]-1]
    dvdx = np.diff(vvel,axis=1)[:uvel.shape[0]-1,:uvel.shape[1]-1]/dx[:uvel.shape[0]-1,:uvel.shape[1]-1]
    return np.multiply(-dudy + dvdx,gsw.f(yc[:uvel.shape[0]-1,:uvel.shape[1]-1]))

def uvFromPhi(phi,dx,dy,lats):
    u = np.diff(phi,axis=0)[:phi.shape[0]-1,:phi.shape[1]-1]/dy[:phi.shape[0]-1,:phi.shape[1]-1]
    v = np.diff(phi,axis=1)[:phi.shape[0]-1,:phi.shape[1]-1]/dx[:phi.shape[0]-1,:phi.shape[1]-1]
    return np.multiply(np.asarray(u),1/gsw.f(lats[:phi.shape[0]-1,:phi.shape[1]-1])),np.multiply(np.asarray(v),1/gsw.f(lats[:phi.shape[0]-1,:phi.shape[1]-1]))


def uvFromPsi(psi,dx,dy,lats):
    u = -np.diff(psi,axis=0)[:psi.shape[0]-1,:psi.shape[1]-1]/dy[:psi.shape[0]-1,:psi.shape[1]-1]
    v = np.diff(psi,axis=1)[:psi.shape[0]-1,:psi.shape[1]-1]/dx[:psi.shape[0]-1,:psi.shape[1]-1]
    return np.multiply(np.asarray(u),1/gsw.f(lats[:psi.shape[0]-1,:psi.shape[1]-1])),np.multiply(np.asarray(v),1/gsw.f(lats[:psi.shape[0]-1,:psi.shape[1]-1]))


def outlierCheck(div,i,j):
    if j!=div.shape[1]-1:
        a = np.isnan(div[i][j+1])
    else:
        a =True
    if j !=0:
        b = np.isnan(div[i][j-1])
    else:
        b=True
    if i!=div.shape[0]-1:
        c = np.isnan(div[i+1][j])
    else:
        c=True
    if i!=0:
        d = np.isnan(div[i-1][j])
    else:
        d= True
    return a and b and c and d

def retrieveColRef(c,coord):
    if coord not in c:
        c[coord] = c["m"]
        c["m"] = c["m"]+1
    return c[coord]




def doubleLaplacian(div,dx,dy):
    print(div.shape)
    div = np.asarray(div)
    phi = np.full_like(div,np.nan)
    sol = []
    colref = {"m":0}
    result = []
    once = True
    for i in range(div.shape[0]):
        for j in range(div.shape[1]):
            if not outlierCheck(div,i,j) and ~np.isnan(div[i][j]):
                result.append(div[i][j])
                center = retrieveColRef(colref,(i,j))
                if i != div.shape[0]-1 and ~np.isnan(div[i+1][j]):
                    iplus1 = retrieveColRef(colref,(i+1,j))
                elif i != div.shape[0]-1 and np.isnan(div[i+1][j]):
                    iplus1 = retrieveColRef(colref,(i+1,j))
                    result.append(0)
                    row = np.zeros(colref["m"])
                    row[iplus1] = 1/(dy[i][j])
                    row[center] = -1/(dy[i][j])
                    sol.append(np.asarray(row))
                else:
                    iplus1 = None
                if i!= 0 and ~np.isnan(div[i-1][j]):
                    iminus1 = retrieveColRef(colref,(i-1,j))
                elif i!= 0 and np.isnan(div[i-1][j]):
                    iminus1 = retrieveColRef(colref,(i-1,j))
                    result.append(0)
                    row = np.zeros(colref["m"])
                    row[iminus1] = -1/(dy[i-1][j])
                    row[center] = 1/(dy[i-1][j])
                    sol.append(np.asarray(row))
                else:
                    iminus1 = None
                if j!=div.shape[1]-1 and ~np.isnan(div[i][j+1]):
                    jplus1 = retrieveColRef(colref,(i,j+1))
                elif j!=div.shape[1]-1 and  np.isnan(div[i][j+1]):
                    jplus1 = retrieveColRef(colref,(i,j+1))
                    result.append(0)
                    row = np.zeros(colref["m"])
                    row[jplus1] = 1/(dx[i][j])
                    row[center] = -1/(dx[i][j])
                    sol.append(np.asarray(row))
                else:
                    jplus1 = None
                if j!=0 and ~np.isnan(div[i][j-1]):
                    jminus1 = retrieveColRef(colref,(i,j-1))
                elif j!=0 and np.isnan(div[i][j-1]):
                    jminus1 = retrieveColRef(colref,(i,j-1))
                    result.append(0)
                    row = np.zeros(colref["m"])
                    row[jminus1] = -1/(dx[i][j-1])
                    row[center] = 1/(dx[i][j-1])
                    sol.append(np.asarray(row))
                else:
                    jminus1 = None
                row = np.zeros(colref["m"])
                row[center] = 0
                if jminus1:
                    row[jminus1] = 1/((dx[i][j-1])*((dx[i][j-1]+dx[i][j])/2))
                    row[center] = row[center]  - 1/((dx[i][j-1])*((dx[i][j-1]+dx[i][j])/2))
                if jplus1:
                    row[jplus1] = 1/((dx[i][j])*((dx[i][j-1]+dx[i][j])/2))
                    row[center] = row[center] -1/((dx[i][j])*((dx[i][j-1]+dx[i][j])/2)) 
                if iminus1:
                    row[iminus1] = 1/((dy[i-1][j])*((dy[i-1][j]+dy[i][j])/2))
                    row[center] = row[center] - 1/((dy[i-1][j])*((dy[i-1][j]+dy[i][j])/2))
                if iplus1:
                    row[iplus1] = 1/((dy[i][j])*((dy[i-1][j]+dy[i][j])/2))
                    row[center] = row[center] - 1/((dy[i][j])*((dy[i-1][j]+dy[i][j])/2)) 

                sol.append(np.asarray(row))
    start = np.zeros(colref["m"])
    start[int(len(start)/2)]=1
    sol.append(start)
    result.append(0)
    sol, misc = inverttools.combineAs([colref["m"]],[0],np.asarray(sol))
    j,[VT, D, U] = inverttools.SVDdecomp(sol,n_elements=int(sol.shape[1]))
    print(inverttools.condition(D))
    c = np.asarray(result).T
    prime = np.matmul(j,c)
    print(prime)
    for k in colref.keys():
        if k !="m":
            phi[k[0]][k[1]] = prime[colref[k]]
    return phi


def boxIntegral(uvel,vvel,psi,botleft,topright,dx,dy):
    line1 = psi[botleft[0]:topright[0],botleft[1]].dot(dx[botleft[0]:topright[0],botleft[1]])
    line2 = psi[topright[0],botleft[1]:topright[1]].dot(dy[topright[0],botleft[1]:topright[1]])
    line3 = psi[botleft[0]:topright[0],topright[1]].dot(dx[botleft[0]:topright[0],topright[1]])
    line4 = psi[botleft[0],botleft[1]:topright[1]].dot(dy[botleft[0],botleft[1]:topright[1]])
    print(line1 +line2 + line3 +line4)

def helmholtzStitch(prefix1,prefix2):
    uset1 = np.asarray(xr.open_dataset('ecco/TILEDATA/'+prefix1+'EVEL.nc',decode_times=False)["EVEL"])
    vset1 = np.asarray(xr.open_dataset('ecco/TILEDATA/'+prefix1+'NVEL.nc',decode_times=False)["NVEL"])
    uset2 = np.asarray(xr.open_dataset('ecco/TILEDATA/'+prefix2+'EVEL.nc',decode_times=False)["EVEL"])
    vset2 = np.asarray(xr.open_dataset('ecco/TILEDATA/'+prefix2+'NVEL.nc',decode_times=False)["NVEL"])


    depths = xr.open_dataset('ecco/TILEDATA/'+prefix1+'NVEL.nc',decode_times=False)["dep"].values
    print(depths[42:])
    ecco_grid1 = xr.open_dataset('ecco/TILEDATA/'+prefix1+'GRID.nc')
    lons1 = ecco_grid1.XC
    lats1 = ecco_grid1.YC
    dy1 = ecco_grid1.DYC
    dx1 = ecco_grid1.DXC

    ecco_grid2 = xr.open_dataset('ecco/TILEDATA/'+prefix2+'GRID.nc')
    lons2 = ecco_grid2.XC
    lats2 = ecco_grid2.YC
    dy2 = ecco_grid2.DYC
    dx2 = ecco_grid2.DXC

    uset1 = np.nanmean(uset1,axis=0)
    vset1 = np.nanmean(vset1,axis=0)
    uset2 = np.nanmean(uset2,axis=0)
    vset2 = np.nanmean(vset2,axis=0)
    dx = np.concatenate((dx1,dx2))
    dy = np.concatenate((dy1,dy2))
    lats = np.concatenate((lats1,lats2))
    lons = np.concatenate((lons1,lons2)).T
    ud1,vd1 = [],[]
    ud2,vd2 = [],[]

    for j in Bar("helmholtz decomp").iter(range(uset1.shape[0])):
        uset = np.concatenate((uset1[j],uset2[j]))
        vset = np.concatenate((vset1[j],vset2[j]))


        div = generateDivMatrix(uset.T,vset.T,dx.T,dy.T,lats.T)
        # fig, (ax1,ax2) = plt.subplots(1,2)
        # ax1.imshow(div)
        # ax2.imshow(divb)
        # plt.show()

        div[np.logical_and(lons[:lons.shape[0]-1,:lons.shape[1]-1]>-110,lons[:lons.shape[0]-1,:lons.shape[1]-1]<150)] = np.nan
        phi = doubleLaplacian(div,dx.T,dy.T)
        ud,vd = uvFromPhi(phi,dx.T,dy.T,np.asarray(lats).T)

        ud1.append(ud[0:90,0:89].T)
        vd1.append(vd[0:90,0:89].T)
        ud2.append(ud[0:89,90:179].T)
        vd2.append(vd[0:89,90:179].T)
    return (ud1,ud2),(vd1,vd2)

#read nc files, load into profiles and save into pickle
def generateProfilesNative(prefix,coordFilter,savepath='data/eccoprofiles.pickle',knownu=None,knownv=None):
    prefixToMixCoord={"BRASILEAST":(90,0),"BRASILWEST":(90,270),"NEPBWEST":(180,180),"NEPBEAST":(180,270)}
    saltset = xr.open_dataset('ecco/TILEDATA/'+prefix+'SALT.nc',decode_times=False)
    thetaset = xr.open_dataset('ecco/TILEDATA/'+prefix+'THETA.nc',decode_times=False)
    uset = xr.open_dataset('ecco/TILEDATA/'+prefix+'EVEL.nc',decode_times=False)
    vset = xr.open_dataset('ecco/TILEDATA/'+prefix+'NVEL.nc',decode_times=False)
    depths = thetaset.dep
    ecco_grid = xr.open_dataset('ecco/TILEDATA/'+prefix+'GRID.nc')
    lons = ecco_grid.XC
    lats = ecco_grid.YC
    dy = ecco_grid.DYC
    dx = ecco_grid.DXC
    # psi = genPsi(uset["EVEL"],vset["NVEL"],dx,dy,lons,lats,10)
    # psi = doubleLaplacian(curl,dx,dy)

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
                        data["sal"].append(np.mean(s))
                        if not knownu:
                            data["knownu"].append(u[i][j])
                            data["knownv"].append(v[i][j])
                        else:
                            data["knownu"].append(knownu[depthindex][i][j])
                            data["knownv"].append(knownv[depthindex][i][j])

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
# ud,vd = helmholtzStitch("NEPBWEST","NEPBEAST")
# with open("data/udvd.pickle", 'wb') as outfile:
    #pickle.dump([ud,vd],outfile)
with open("data/udvd.pickle", 'rb') as outfile:
    [ud,vd] = pickle.load( outfile)
ud = np.asarray(ud)
vd = np.asarray(vd)
udprefixwest = []
udprefixeast = []
vdprefixwest = []
vdprefixeast = []
print(vd.shape)
print(ud.shape)
for j in range(len(ud[1])):
    ud1 = ud[0][j]
    vd1 = vd[0][j]
    ud2 = ud[1][j]
    vd2 = vd[1][j]
    ud1 = np.vstack([ud1,ud1[-1]])
    ud1 = np.hstack((ud1, np.tile(ud1[:, [-1]], 2)))
    ud2 = np.vstack([ud2,ud2[-1]])
    ud2 = np.vstack([ud2,ud2[-1]])
    ud2 = np.hstack((ud2, np.tile(ud2[:, [-1]], 2)))
    vd1 = np.vstack([vd1,vd1[-1]])
    vd1 = np.hstack((vd1, np.tile(vd1[:, [-1]], 2)))
    vd2 = np.vstack([vd2,vd2[-1]])
    vd2 = np.vstack([vd2,vd2[-1]])
    vd2 = np.hstack((vd2, np.tile(vd2[:, [-1]], 2)))
    print("vd1",vd1.shape)
    print("vd2",vd2.shape)
    print("ud1",ud1.shape)
    print("ud2",ud2.shape)
    udprefixwest.append(ud1)
    udprefixeast.append(ud2)
    vdprefixwest.append(vd1)
    vdprefixeast.append(vd2)

p1  = generateProfilesNative("NEPBWEST",nepbRestrict,"data/ecconepbprofiles.pickle",knownu=udprefixwest,knownv=vdprefixwest)
p2  = generateProfilesNative("NEPBEAST",nepbRestrict,"data/ecconepbprofiles.pickle",knownu=udprefixeast,knownv=vdprefixeast)
with open("data/ecconepbprofiles.pickle", 'wb') as outfile:
    pickle.dump(p1+p2, outfile)
