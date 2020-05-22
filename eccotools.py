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
                if i==0:
                    i+=1
                if j==0:
                    j+=1
                if j==div.shape[1]-1:
                    j+=-1
                if i==div.shape[0]-1:
                    i+=-1
                if ~np.isnan(div[i+1][j]):
                    iplus1 = retrieveColRef(colref,(i+1,j))
                else:
                    iplus1 = retrieveColRef(colref,(i+1,j))
                    result.append(0)
                    row = np.zeros(colref["m"])
                    row[iplus1] = 1/(dy[i][j])
                    row[center] = -1/(dy[i][j])
                    sol.append(np.asarray(row))
                    iplus1 = None
                if ~np.isnan(div[i-1][j]):
                    iminus1 = retrieveColRef(colref,(i-1,j))
                else:
                    iminus1 = retrieveColRef(colref,(i-1,j))
                    result.append(0)
                    row = np.zeros(colref["m"])
                    row[iminus1] = -1/(dy[i-1][j])
                    row[center] = 1/(dy[i-1][j])
                    sol.append(np.asarray(row))
                    iminus1 = None
                if ~np.isnan(div[i][j+1]):
                    jplus1 = retrieveColRef(colref,(i,j+1))
                else:
                    jplus1 = retrieveColRef(colref,(i,j+1))
                    result.append(0)
                    row = np.zeros(colref["m"])
                    row[jplus1] = 1/(dx[i][j])
                    row[center] = -1/(dx[i][j])
                    sol.append(np.asarray(row))
                    jplus1 = None
                if ~np.isnan(div[i][j-1]):
                    jminus1 = retrieveColRef(colref,(i,j-1))
                else:
                    jminus1 = retrieveColRef(colref,(i,j-1))
                    result.append(0)
                    row = np.zeros(colref["m"])
                    row[jminus1] = -1/(dx[i][j-1])
                    row[center] = 1/(dx[i][j-1])
                    sol.append(np.asarray(row))
                    jminus1 = None
                row = np.zeros(colref["m"])
                if jminus1:
                    row[jminus1] = 1/((dx[i][j-1])*((dx[i][j-1]+dx[i][j])/2))
                if jplus1:
                    row[jplus1] = 1/((dx[i][j])*((dx[i][j-1]+dx[i][j])/2))
                if iminus1:
                    row[iminus1] = 1/((dy[i-1][j])*((dy[i-1][j]+dy[i][j])/2))
                if iplus1:
                    row[iplus1] = 1/((dy[i][j])*((dy[i-1][j]+dy[i][j])/2))
                row[center] = -1/((dy[i][j])*((dy[i-1][j]+dy[i][j])/2)) + -1/((dy[i-1][j])*((dy[i-1][j]+dy[i][j])/2))
                row[center] = row[center] + -1/((dx[i][j])*((dx[i][j-1]+dx[i][j])/2)) - 1/((dx[i][j-1])*((dx[i][j-1]+dx[i][j])/2))
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
    # psi = genPsi(uset["EVEL"],vset["NVEL"],dx,dy,lons,lats,10)
    print(depths[30])
    print(uset["EVEL"].shape)
    div = generateDivMatrix(np.nanmean(uset["EVEL"],axis=0)[10].T,np.nanmean(vset["NVEL"],axis=0)[10].T,dx.T,dy.T,lats.T)
    curl = generateCurlMatrix(np.nanmean(uset["EVEL"],axis=0)[10].T,np.nanmean(vset["NVEL"],axis=0)[10].T,dx.T,dy.T,lats.T)
    fig, (ax1,ax2) = plt.subplots(1,2)
    ax1.imshow(div)
    ax2.imshow(curl)
    plt.show()
    phi = doubleLaplacian(div,dx,dy)
    # psi = doubleLaplacian(curl,dx,dy)
    fig, (ax1,ax2,ax3) = plt.subplots(1,3)
    psi = phi
    ax1.imshow(phi)
    ax2.imshow(psi)
    ax3.imshow(phi+psi)
    plt.show()
    fig, ((ax11,ax12,ax13,ax14),(ax21,ax22,ax23,ax24)) = plt.subplots(2,4)

    ud,vd = uvFromPhi(phi,dx,dy,np.asarray(lats))
    ur,vr = uvFromPsi(psi,dx,dy,np.asarray(lats))
    curl2 = generateCurlMatrix(ud,vd,dx.T,dy.T,lats.T)
    div2 = generateDivMatrix(ud,vd,dx.T,dy.T,lats.T)
    div = generateDivMatrix(np.nanmean(uset["EVEL"],axis=0)[10].T,np.nanmean(vset["NVEL"],axis=0)[10].T,dx.T,dy.T,lats.T)
    curl = generateCurlMatrix(np.nanmean(uset["EVEL"],axis=0)[10].T,np.nanmean(vset["NVEL"],axis=0)[10].T,dx.T,dy.T,lats.T)

    c11 = ax11.imshow(ud)
    ax11.title.set_text("U divergent")
    c12 = ax12.imshow(vd)
    ax12.title.set_text("V divergent")
    c21 = ax21.imshow(np.nanmean(uset["EVEL"],axis=0)[10].T)
    ax21.title.set_text("U")
    c22 = ax22.imshow(np.nanmean(vset["NVEL"],axis=0)[10].T)
    ax22.title.set_text("V")
    c13 = ax13.imshow(curl2)
    ax13.title.set_text("Curl of Ud,Vd")
    c23 = ax23.imshow(div2)
    ax23.title.set_text("Divergence of Ud,Vd")
    c14 = ax14.imshow(curl)
    ax14.title.set_text("Curl of U,V")
    c24 = ax24.imshow(div)
    ax24.title.set_text("Divergence of U,V")
    fig.colorbar(c11,ax=ax11)
    fig.colorbar(c12,ax=ax12)
    fig.colorbar(c21,ax=ax21)
    fig.colorbar(c22,ax=ax22)
    fig.colorbar(c13,ax=ax13)
    fig.colorbar(c23,ax=ax23)
    fig.colorbar(c14,ax=ax14)
    fig.colorbar(c24,ax=ax24)
    plt.show()
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
                        data["knownu"].append(ud[i][j])
                        data["knownv"].append(vd[i][j])

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
