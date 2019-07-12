import csv
from mpl_toolkits.basemap import Basemap
from numpy import linspace
import matplotlib.pyplot as plt
from numpy import meshgrid
import numpy as np
from netCDF4 import Dataset
import json
from geopy.distance import geodesic
from geopy.distance import great_circle
from profile import Profile
import random
from scipy.interpolate import Rbf
import pyproj
import bathtools
import copy
import gsw
import pygam
import pickle
from gswmatlab.pyinterface import geo_strf_isopycnal
import scipy.io as sio

def extractProfiles(fnames):
    ##Load JSON data into profile objects
    ##and return index of deepest one
    profiles = []
    deepestindex =[] 
    deepestdepth = 0 
    print(fnames)
    for fname in fnames:
        json_file = open(fname,'r') 
        data = json.load(json_file)
        for p in data.keys():
            profile = Profile(p,data[p])
            if len(profile.ipres)>0:
                profiles.append(profile)
                if data[p]["pres"][-1] > deepestdepth:
                    deepestindex = len(profiles)-1
                    deepestdepth=data[p]["pres"][-1]
    return profiles, deepestindex



def deepestProfile(profiles):
    deepestindex=-1
    deepestdepth=-1
    for p in range(len(profiles)):
        if profiles[p].ipres[-1] >deepestdepth:
            deepestdepth=profiles[p].ipres[-1]
            deepestindex=p
    return deepestindex


def cruiseCount(profiles):
    cruises ={"null"}
    for profile in profiles:
        if profile.cruise not in cruises:
            cruises.add(profile.cruise)
    return cruises

def cruiseSearch(profiles,cruisename,year=None):
    results =[]
    for p in profiles:
        if p.cruise == cruisename:
            if year:
                if year == p.time.year:
                    results.append(p)
            else:
                results.append(p)
    return results


def transArcticSearch(profiles):
    results = []
    for name in cruiseCount(profiles):
        westernprofile = False
        easternprofile = False
        for profile in cruiseSearch(profiles,name):
            if profile.lat > 70 and abs(profile.lon)<20:
                easternprofile =True
            if profile.lat > 70 and abs(profile.lon)>160:
                westernprofile = True
        if westernprofile and easternprofile:
            results.append(name)
    return results

        


def extractProfilesMonths(fnames,months):
    ##Load JSON data into profile objects
    ##and return index of deepest one
    profiles = []
    deepestindex =[] 
    deepestdepth = 0 
    for fname in fnames:
        json_file = open(fname) 
        data = json.load(json_file)
        for p in data.keys():
            profile = Profile(p,data[p])
            if profile.time.month in months :
                if len(profile.ipres)>0:
                    profiles.append(profile)
                    if data[p]["pres"][-1] > deepestdepth:
                        deepestindex = len(profiles)-1
                        deepestdepth=data[p]["pres"][-1]
    return profiles, deepestindex

def extractProfilesBox(fnames,lonleft,lonright,latbot,lattop):
    ##Load JSON data into profile objects
    ##and return index of deepest one
    profiles = []
    deepestindex =[] 
    deepestdepth = 0 
    for fname in fnames:
        json_file = open(fname) 
        data = json.load(json_file)
        for p in data.keys():
            profile = Profile(p,data[p])
            if latbot <= profile.lat <= lattop and lonleft <= profile.lon <= lonright:
                if len(profile.ipres)>0:
                    profiles.append(profile)
                    if data[p]["pres"][-1] > deepestdepth:
                        deepestindex = len(profiles)-1
                        deepestdepth=data[p]["pres"][-1]
    return profiles, deepestindex

def removeNorwegianSea(profiles):
    finalprofiles = []
    deepestindex = -1
    deepestdepth = 0
    for pindex in range(len(profiles)):
        p = profiles[pindex]
        if not (p.lat<81 and -23 < p.lon < 20):
            finalprofiles.append(p)
            if p.ipres[-1] > deepestdepth:
                deepestindex = pindex
                deepestdepth = p.ipres[-1]
    return finalprofiles,deepestindex 

def closestIdentifiedNS(profiles,queryprofile,depth,radius):
    minimumdistance = radius
    minprofile = None
    for p in profiles:
        if depth in p.neutraldepth.keys():
            dist = geodesic((queryprofile.lat,queryprofile.lon),(p.lat,p.lon)).km
            if dist<minimumdistance:
                minimumdistance = dist
                minprofile = p
    #print(minprofile.neutraldepth)
    return minprofile

def profileInBox(profiles,lonleft,lonright,latbot,lattop):
    results=[]
    for p in profiles:
        if lonleft< p.lon <lonright and latbot < p.lat < lattop:
            results.append(p)
    return results

def emptySurface():
    return {"lats":[],"lons":[],"ids":[],"data":{"pres":[],"t":[],"s":[]}}

def peerSearch(profiles,deepestindex,depth,profilechoice,radius=500):
    surfaces = {}
    surfaces[depth]= emptySurface()
    profilechoice.neutraldepth[depth] = depth
    print(profilechoice.eyed)
    references=[]
    references.append(profilechoice)
    print(len(profiles))
    while len(profiles)>0:
        foundcounter = 0
        for p in profiles.copy():
            closest = closestIdentifiedNS(references,p,depth,radius)
            if closest:
                #print(closest.neutraldepth)
                ns = closest.neutralDepth(p,closest.neutraldepth[depth],depthname=depth,searchrange=100) 
                if ns != None:
                    surfaces[depth]["lons"].append(p.lon)
                    surfaces[depth]["lats"].append(p.lat)
                    surfaces[depth]["data"]["pres"].append(ns)
                    surfaces[depth]["ids"].append(p.eyed)
                    profiles.remove(p)
                    foundcounter +=1
                    references.append(p)
        if foundcounter ==0:
            break

        print("found: ,",foundcounter,"left: ", len(profiles))
        #plotProfiles(references,"ITS SPREADING",specialprofile = profilechoice)
    return surfaces

def search(profiles,deepestindex):
    #Lets look for neutral surfaces every 200 dbar below 1000 dbar
    deeprange = range(1000,max(profiles[deepestindex].ipres),200)
    #A dictionary mapping neutral surface pressures to pressures at lat,lon
    surfaces = {}
    for r in deeprange:
        #note: an annoying format to store this information but easily graphable
        surfaces[r]=emptySurface()
    #Ever profile
    for j in range(len(profiles)):
        #every neutral surface
        for r in deeprange:
            ns = profiles[deepestindex].neutralDepth(profiles[j],r) 
            if ns != None:
                surfaces[r]["lons"].append(profiles[j].lon)
                surfaces[r]["lats"].append(profiles[j].lat)
                surfaces[r]["data"]["pres"].append(ns)
    return surfaces

def addDataToSurfaces(profiles,surfaces,stdevs,debug=True):
    if debug:
        print("adding data to surfaces")
    tempSurfs = {}
    negativecount = 0
    for k in surfaces.keys():
        tempSurf = emptySurface()
        for l in range(len(surfaces[k]["lons"])):
            p = getProfileById(profiles,surfaces[k]["ids"][l])
            p.interpolate()
            t,s = p.atPres(surfaces[k]["data"]["pres"])
            pv = p.potentialVorticity(surfaces[k]["data"]["pres"][l],debug=False)
            if pv and pv<0:
                negativecount +=1 
            if t and s and pv and p and pv != np.Inf and pv != np.nan and not np.isnan(t):
                tempSurf["lons"].append(surfaces[k]["lons"][l])
                tempSurf["lats"].append(surfaces[k]["lats"][l])
                tempSurf["data"]["pres"].append(-surfaces[k]["data"]["pres"][l])
                tempSurf["data"]["t"].append(t)
                tempSurf["data"]["s"].append(s)
                tempSurf["data"]["pv"].append(pv)
                tempSurf["ids"].append(surfaces[k][3][l])
        #for j in tempSurf["data"].keys():
            #m = np.median(tempSurf["data"][j])
            #s = np.std(tempSurf["data"][j])
            #a = np.asarray(np.where(abs(tempSurf["data"][j] -m)>stdevs*s))
            ##np.asarray(tempSurf[2][j])[a[0]] == np.nan
            

        if len(tempSurf[0])>5:
            tempSurfs[k] = tempSurf

        print("ns: ",k," negative count: ",negativecount)
    return tempSurfs


def getProfileById(profiles,eyed):
    for p in profiles:
        if p.eyed == eyed:
            return p

def filterCruises(profiles,cruisenames):
    finalprofiles = []
    for p in profiles:
        if p.cruise in cruisenames:
            finalprofiles.append(p)
    return finalprofiles

def filterSurfacesByLine(surfaces,lon,radius=20):
    print("filtering by crosssection")
    for k in surfaces.keys():
        distFilter = np.zeros(len(surfaces[k]["lons"]), dtype=bool)
        for index in range(len(surfaces[k]["lons"])):
            p = (surfaces[k]["lats"][index],surfaces[k]["lons"][index])
            if ((np.cos(p[0]*(np.pi/180)))*abs(((p[1]+180)%180)-lon)*111)<radius:
                #print(((np.cos(lpoint[0]*(np.pi/180)))*abs(p[1]-lpoint[1])*111))
                distFilter[index] = True
        surfaces[k]["lons"] = np.asarray(surfaces[k]["lons"])[distFilter]
        surfaces[k]["lats"] = np.asarray(surfaces[k]["lats"])[distFilter]
        for j in surfaces[k]["data"].keys():
            surfaces[k]["data"][j] = np.asarray(surfaces[k]["data"][j])[distFilter]
        if len(surfaces[k]["ids"])>0:
            surfaces[k]["ids"] = np.asarray(surfaces[k]["ids"])[distFilter]
    return surfaces
        
def plotProfile(p):
    fig, (ax1,ax2) = plt.subplots(1,2)
    ax1.plot(p.itemps,p.ipres)
    ax2.plot(p.isals,p.ipres)
    ax1.invert_yaxis()
    ax2.invert_yaxis()
    fig.suptitle(str(p.eyed)+" lat: "+str(p.lat)+" lon: "+ str(p.lon))
    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
    plt.show()


def addXYToSurfaces(surfaces,debug=True):
    if debug:
        print("converting surfaces to xyz")
    newsurfaces = {}
    for k in surfaces.keys():
        x,y = homemadeXY(surfaces[k]["lons"],surfaces[k]["lats"])
        surfaces[k]["x"]=x
        surfaces[k]["y"]=y
    return surfaces

def deduplicateXYZ(x,y,z):
    return zip(*set(list(zip(x,y,z))))

def removeDiscontinuities(surface,radius=10,debug=True):
    if debug:
        print("removing discontinuities")
    x=np.asarray(surface["x"])
    y=np.asarray(surface["y"])
    z=np.asarray(surface["data"]["pres"])
    final = np.zeros(x.shape)
    #print(x)
    for i in range(len(x)):
        if final[i] == False:
            r = np.sqrt((x- x[i])**2 + (y - y[i])**2)
            inside = r<radius*1000
            inside[i]=False
            s = 0
            counter = 0
            for t in range(len(inside)):
                if t:
                    s+=z[i]
                    counter+=1
            if counter >0:
                z[i] = s/counter
                    
            #print(r)
            if np.count_nonzero(final) == 0  :
                final =inside
            else:
                final = np.logical_or(final,inside)
    final = np.invert(final)
    for k in surface.keys():
        if k == "data":
            for d in surface[k]:
                surface[k][d] = np.asarray(surface[k][d])[final]
        else:
            surface[k] = np.asarray(surface[k])[final]
    return surface

def createMesh(n,xvals,yvals):
    return np.meshgrid(np.linspace(np.min(xvals),np.max(xvals),n), np.linspace(np.min(yvals),np.max(yvals),n),indexing="xy")

def generateMaskedMesh(x,y,radius=100):
    xi,yi = createMesh(100,x,y)
    final = np.zeros(xi.shape)
    for i in range(len(x)):
        r = np.sqrt((xi- x[i])**2 + (yi - y[i])**2)
        inside = r<radius*1000
        if np.count_nonzero(final) == 0  :
            final =inside
        else:
            final = final+inside
    for i in range(len(final[0])):
        if final[0][i] > 2:
            final[0][i]=True
        else:
            final[0][i] = False
    return xi[final],yi[final]

def removeOutlierSurfaces(surfaces,stdevs=2):
    for k in surfaces.keys():
        m = np.median(surfaces[k]["data"]["pres"])
        s = np.std(surfaces[k]["data"]["pres"])
        filt = np.where(np.abs(surfaces[k]["data"]["pres"]-m)<s*stdevs)
        for field in surfaces[k].keys():
            if field=="data":
                for datafield in surfaces[k][field].keys():
                    surfaces[k][field][datafield]=surfaces[k][field][datafield][filt]
                surfaces[k][field]=surfaces[k][field][filt]
            else:
                surfaces[k][field]=surfaces[k][field][filt]
    return surfaces


def interpolateSurface(surface,debug=True):
    #print("######")
    interpsurf={}
    X = np.zeros((len(surface["x"]),2))
    X[:,0]=surface["x"]
    X[:,1]=surface["y"]
    xi,yi = generateMaskedMesh(surface["x"],surface["y"])
    interpdata={}
    interpsurf["x"] =xi
    interpsurf["y"] =yi
    for k in surface["data"].keys():
        gam = pygam.LinearGAM(pygam.te(0,1)).fit(X,surface["data"][k])
        Xgrid = np.zeros((yi.shape[0],2))
        Xgrid[:,0] = xi
        Xgrid[:,1] = yi
        interpdata[k] = gam.predict(Xgrid)
    interpsurf["data"] = interpdata
    interpsurf = addLatLonToSurface(interpsurf)
    return interpsurf

def homemadeXY(lon,lat):
    x=[]
    y=[]
    for i in range(len(lat)):
        theta = np.deg2rad(lon[i])
        r = ((90-lat[i]) *111*1000)
        x.append(r*np.cos(theta))
        y.append(r*np.sin(theta))
    return np.asarray(x),np.asarray(y)

def addLatLonToSurface(surface,debug = True):
    if debug:
        print("converting xyz back to lat lon a")
    lat = 90-(np.sqrt((surface["x"]**2+surface["y"]**2))/111000.0)
    lon = np.degrees(np.arctan2(surface["y"],surface["x"]))
    #print("lat: ",lat, " lon: ",lon)
    surface["lons"]=lon
    surface["lats"]=lat
    return surface
     
        
def findNeighboringPoints(profiles,lat,lon,radius=30):
    lats =[]
    lons = []
    eyeds= []
    for p in profiles:
        if geodesic((p.lat,p.lon),(lat,lon)).km<radius:
            lats.append(p.lat)
            lons.append(p.lon)
            eyeds.append(p.eyed)
    plt.scatter(lons,lats)
    plt.show()

def plotASpiral(profiles,center=None,x=None,y=None):
    fig,ax = plt.subplots(1,1) 
    center = getProfileById(profiles,"120389")
    x = getProfileById(profiles,"120387")
    y = getProfileById(profiles,"114688")
    dx = geodesic((x.lat,x.lon),(center.lat,center.lon)).m
    dy = geodesic((x.lat,x.lon),(center.lat,center.lon)).m
    u = 0
    v = 0
    f = center.f
    us = [0]
    vs = [0]
    for i in range(100,1700,100):
        dpdx = (x.densityAtPres(i)-center.densityAtPres(i))/dx
        dpdy = (y.densityAtPres(i)-center.densityAtPres(i))/dy
        v = (9.8/-f)*(dpdx)
        u = (9.8/f)*(dpdy)
        us.append(u)
        vs.append(v)
        plt.plot([0,u],[0,v])
        ax.annotate(i, (u, v))
    #plt.scatter(us,vs)
    #plt.plot(us,vs,c="r")
    plt.show()

def generateNeighborsList(x,y):
    print("calculating neighbors")
    ref = {}
    for center in range(len(x)):
        dist = ((x-x[center])**2 + (y-y[center])**2)
        ref[center]=np.argsort(dist)[1:5]
    return ref
                     
def componentDistance(surfaces,k,i1,i2):
    #x = surfaces[k][0][i1] - surfaces[k][0][i2]
    #if abs(x) > abs(360 - (x+360)%360):
        #x = np.sign(x)*(360-(x+360)%360)
    if  (surfaces[k]["lons"][i1]+180 ) > (surfaces[k]["lons"][i2]+180):
        x = surfaces[k]["lons"][i1]+180 - (surfaces[k]["lons"][i2]+180)
        if x>180:
            x = -(360-x)
    else:
        x = surfaces[k]["lons"][i1]+180 - (surfaces[k]["lons"][i2]+180)
        if x < -180:
            x= -(-360-x)
    x=x*np.cos(np.deg2rad(surfaces[k]["lats"][i2]))*111000.0
    #print(surfaces[k][0][i1],surfaces[k][0][i2],x)

    y = (surfaces[k]["lats"][i1]-surfaces[k]["lats"][i2])*111000.0
    return x,y


def addPrimeToSurfaces(surfaces,neighbors,debug=False):
    for k in surfaces.keys():
        surfaces[k]["data"]["uz"] = np.zeros(len(surfaces[k]["lons"]))
        surfaces[k]["data"]["vz"] = np.zeros(len(surfaces[k]["lons"]))
    alldxs = []
    for k in neighbors.keys():
        print("adding primes to: ",k)
        for r in neighbors[k].keys():
            #alright so here k is our NS
            #r is the index of the point for which we are calculating these prime values
            #adjacent is the list of adjacent points
            adjacent = neighbors[k][r]
            dxsum = 0
            dysum = 0
            dxs = []
            dys = []
            dhs = []
            dists = 0
            for i in adjacent:
                dx,dy = componentDistance(surfaces,k,i,r)
                dxs.append(dx)
                dys.append(dy)
            dxindexs = [np.argmin(dxs),np.argmax(dxs)]
            dyindexs = [np.argmin(dys),np.argmax(dys)]
            dxfinal,b = componentDistance(surfaces,k,adjacent[dxindexs[1]],adjacent[dxindexs[0]])
            b,dyfinal = componentDistance(surfaces,k,adjacent[dyindexs[1]],adjacent[dyindexs[0]])
            #print(r,adjacent)
            dhx = surfaces[k]["data"]["pres"][adjacent[dxindexs[1]]] - surfaces[k]["data"]["pres"][adjacent[dxindexs[0]]]
            #dhx = surfaces[k][2][0][adjacent[dxindexs[1]]]
            dhy = surfaces[k]["data"]["pres"][adjacent[dyindexs[1]]]-surfaces[k]["data"]["pres"][adjacent[dyindexs[0]]]
            #dhy = surfaces[k][2][0][adjacent[dyindexs[1]]] - surfaces[k][2][0][adjacent[dyindexs[0]]]
            dhdtheta = dhx/dxfinal
            dhdr = dhy/dyfinal
            #surfaces[k][2][4][r] = dhdtheta *(1/((90-surfaces[k][2][0][r])*111000))*(1/np.tan(np.deg2rad(surfaces[k][2][0][r])))
            surfaces[k]["data"]["uz"][r] = dhdtheta
            surfaces[k]["data"]["vz"][r] = dhdr 
    return surfaces

def addStreamFunc(surfaces,profiles):
    neutraldepths={}
    for k in surfaces.keys():
        for i in range(len(surfaces[k]["ids"])):
            if surfaces[k]["ids"][i] not in neutraldepths:
                neutraldepths[surfaces[k]["ids"][i]] =[[],[]]
            neutraldepths[surfaces[k]["ids"][i]][0].append(k)
            neutraldepths[surfaces[k]["ids"][i]][1].append(surfaces[k]["data"]["pres"][i])
    s = []
    t = []
    ip = []
    p_ref = []
    ns = []
    ns_p = []
    ks = []
    count =0
    for k in neutraldepths.keys():
        count+=1
        if count > 2000:
            break
        p = getProfileById(profiles,k)
        s.append(p.sals)
        t.append(p.temps)
        ip.append(p.pres)
        #p_ref.append([10.1235]*len(p.ipres))
        if len(np.unique(p.ipres))!=len(p.ipres):
            print("NOOOOOOOO")
        p_ref.append([0]*len(p.ipres))
        ns.append(neutraldepths[k][0])
        ns_p.append(neutraldepths[k][1])
        ks.append(k)
    
    results = geo_strf_isopycnal(s,t,ip,p_ref,ns,ns_p,ks)
    
    #with open('data/geoisopycnal.pickle', 'wb') as outfile:
        #pickle.dump(results, outfile)
    #print(results)

        
    

def addStreamFuncFromFile(surfaces,profiles,isopycnalfile,referencefile):
    psi = np.asarray(sio.loadmat(isopycnalfile)["geoisopycnal"]).transpose()
    ids = np.asarray(sio.loadmat(referencefile)["ks"])
    tags = np.asarray(sio.loadmat(referencefile)["ns"]).transpose()
    for k in surfaces.keys():
        surfaces[k][2]["psi"] = np.full_like(surfaces[k]["data"]["pres"]
        for i in range(len(surfaces[k]["lons"])):
            if surfaces[k]["ids"][i] in ids:
                col = np.where(ids == surfaces[k]["ids"][i])[0][0]
                row = np.where(k == tags[col])[0]
                if len(row) >0:
                    print(psi[col][row])
                    surfaces[k]["data"]["psi"].append(psi[col][row[0]])
                else:
                    surfaces[k]["data"]["psi"].append(np.nan)
                
            else:
                surfaces[k]["data"]["psi"].append(np.nan)
    return surfaces

def convertOldSurfaces(surfaces):
    newsurfaces = {}
    for k in surfaces.keys():
        surface = {}
        surface["lons"] = surfaces[k][0]
        surface["lats"] = surfaces[k][1]
        surface["data"] = {"pres":surfaces[k][2][0],"t":surfaces[k][2][1],"s":surfaces[k][2][0],"pv":surfaces[k][2][3]}
        surface["ids"] = surfaces[k][3]
        newsurfaces[k] = surface
    return newsurfaces
